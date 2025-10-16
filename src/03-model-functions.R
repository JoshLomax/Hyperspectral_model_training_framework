# ===================================================================
# 03-model-functions.R
# Individual model fitting functions
# ===================================================================

# BayesA model ------------
fit_bayes_model <- function(df,
                            seed,
                            method = "BayesA",
                            initial_nIter = 30000,
                            initial_burnIn = 10000,
                            initial_thin = 100,
                            max_nIter = 20000,
                            n_folds = 5,
                            convergence_threshold = 1.05,
                            max_attempts = 5) {
  set.seed(seed)
  assess_convergence <- function(fit, thin = 5) {
    if (is.null(fit$ETA[[1]]$b)) {
      return(list(
        effective_size = 0,
        gelman_rubin = Inf,
        geweke = Inf,
        autocorr = 1
      ))
    }
    
    post_burnin_samples <- fit$ETA[[1]]$b
    chain_length <- length(post_burnin_samples)
    n_chains <- 3
    samples_per_chain <- floor(chain_length / n_chains)
    
    # Create chains with equal lengths
    chains <- lapply(1:n_chains, function(i) {
      start_idx <- (i - 1) * samples_per_chain + 1
      end_idx <- i * samples_per_chain
      # Ensure we don't exceed array bounds
      end_idx <- min(end_idx, chain_length)
      as.mcmc(post_burnin_samples[start_idx:end_idx])
    })
    
    mcmc_list <- as.mcmc.list(chains[2:3])
    
    # Calculate diagnostics with default returns on failure
    gelman <- tryCatch({
      gelman.diag(mcmc_list)$psrf[1]
    }, error = function(e) {
      Inf
    })
    
    eff_size <- tryCatch({
      effectiveSize(mcmc_list)
    }, error = function(e) {
      0
    })
    
    geweke <- tryCatch({
      mean(abs(geweke.diag(mcmc_list[[1]])$z), na.rm = TRUE)
    }, error = function(e) {
      Inf
    })
    
    autocorr <- tryCatch({
      mean(autocorr(mcmc_list[[1]], lag = 1), na.rm = TRUE)
    }, error = function(e) {
      1
    })
    
    list(
      effective_size = if (is.na(mean(eff_size)))
        0
      else
        mean(eff_size),
      gelman_rubin = if (is.na(gelman))
        Inf
      else
        gelman,
      geweke = if (is.na(geweke))
        Inf
      else
        geweke,
      autocorr = if (is.na(autocorr))
        1
      else
        autocorr
    )
  }
  
  # parameter adjustment
  adjust_parameters <- function(diagnostics,
                                current_nIter,
                                current_burnIn,
                                current_thin,
                                attempt) {
    # Base adjustment factors
    nIter_factor <- 1.2
    burnIn_factor <- 1.2
    thin_factor <- 1.5
    
    # Analyze convergence issues
    eff_size_issue <- diagnostics$effective_size <= 10
    geweke_issue <- abs(diagnostics$geweke) >= 3
    autocorr_issue <- diagnostics$autocorr >= 0.9
    
    # Adjust parameters based on specific issues
    if (eff_size_issue) {
      # Increase iterations more aggressively if effective sample size is too small
      nIter_factor <- 1.5
      burnIn_factor <- 1.3
    }
    
    if (geweke_issue) {
      # If Geweke statistic indicates poor convergence, increase burn-in
      burnIn_factor <- 1.4
    }
    
    if (autocorr_issue) {
      # If autocorrelation is high, increase thinning interval
      thin_factor <- 2.0
    }
    
    # Apply adjustments with limits
    new_nIter <- min(max_nIter, ceiling(current_nIter * nIter_factor))
    new_burnIn <- min(floor(new_nIter * 0.3),
                      ceiling(current_burnIn * burnIn_factor))
    new_thin <- if (autocorr_issue)
      ceiling(current_thin * thin_factor)
    else
      current_thin
    
    # Ensure minimum values
    new_nIter <- max(new_nIter, 5000)
    new_burnIn <- max(new_burnIn, 1000)
    new_thin <- max(new_thin, 50)
    
    list(
      nIter = new_nIter,
      burnIn = new_burnIn,
      thin = new_thin,
      needs_adjustment = any(eff_size_issue, geweke_issue, autocorr_issue)
    )
  }
  
  model_data <- prepare_model_data(df, seed = seed)
  all_results <- tibble()  # Store all attempts
  best_result <- NULL      # Store best performing iteration
  best_rmse <- Inf         # Track best RMSE
  
  # Initialize parameters
  current_nIter <- initial_nIter
  current_burnIn <- initial_burnIn
  current_thin <- initial_thin
  attempt <- 1
  
  start_time <- Sys.time()
  
  run_model_attempts <- function() {
    while (attempt <= max_attempts) {
      cli::cli_inform(
        "Attempt {attempt}: nIter={current_nIter}, burnIn={current_burnIn}, thin={current_thin}"
      )
      
      # Create y with NA for test set
      y_with_NA <- model_data$y
      y_with_NA[model_data$test_idx] <- NA
      
      # Fit model
      fit <- BGLR(
        y = y_with_NA,
        ETA = list(X = list(X = model_data$X, model = method)),
        nIter = current_nIter,
        burnIn = current_burnIn,
        thin = current_thin,
        verbose = FALSE
      )
      
      # Get diagnostics and metrics
      diagnostics <- assess_convergence(fit, current_thin)
      
      test_metrics <- calculate_metrics(
        model_data$y_test,
        fit$yHat[model_data$test_idx],
        n_predictors = ncol(model_data$X),
        model_data$y_train,
        fit$yHat[model_data$train_idx]
      ) 
      
      # Determine convergence status
      convergence_status <- if (diagnostics$effective_size > 10 &&
                                abs(diagnostics$geweke) < 3 &&
                                diagnostics$autocorr < 0.9) {
        "converged"
      } else {
        "not_converged"
      }
      
      # Create current results
      current_results <- tibble(
        model = "BayesianRegression",
        method = method,
        convergence_status = convergence_status,
        test_metrics,
        convergence_attempts = attempt,
        final_nIter = current_nIter,
        final_burnIn = current_burnIn,
        final_thin = current_thin,
        effective_sample_size = diagnostics$effective_size,
        geweke_stat = diagnostics$geweke,
        autocorrelation = diagnostics$autocorr,
        seed = seed,
      )
      
      # Update all_results
      all_results <<- bind_rows(all_results, current_results)
      
      # Update best result if current iteration is better
      if (test_metrics$train_rmse < best_rmse) {
        best_rmse <<- test_metrics$train_rmse
        best_result <<- current_results
      }
      
      # Check if we should continue
      if (convergence_status == "converged") {
        cli::cli_inform("Convergence achieved at attempt {attempt}")
        return(TRUE)
      }
      
      # Adjust parameters for next attempt
      params <- adjust_parameters(diagnostics,
                                  current_nIter,
                                  current_burnIn,
                                  current_thin,
                                  attempt)
      if (!params$needs_adjustment) {
        cli::cli_inform("No further parameter adjustments needed at attempt {attempt}")
        return(TRUE)
      }
      
      # Update parameters for next iteration
      current_nIter <<- params$nIter
      current_burnIn <<- params$burnIn
      current_thin <<- params$thin
      attempt <<- attempt + 1
    }
    
    cli::cli_warn("Maximum attempts reached without convergence")
    return(FALSE)
  }
  
  tryCatch({
    run_model_attempts()
    
    # Add convergence history to best result
    best_result$convergence_history <- list(
      all_results %>%
        select(
          convergence_attempts,
          train_rmse,
          effective_sample_size,
          geweke_stat,
          autocorrelation
        )
    )
    
    return(best_result)
    
  }, error = function(e) {
    cli::cli_warn("Error in model fitting: {e$message}")
    return(
      tibble(
        model = "BayesianRegression",
        method = method,
        convergence_status = "failed",
        mae = NA_real_,
        mse = NA_real_,
        r2 = NA_real_,
        rmse = NA_real_,
        nrmse = NA_real_,
        n_samples = 0,
        convergence_attempts = attempt,
        final_nIter = current_nIter,
        final_burnIn = current_burnIn,
        final_thin = current_thin,
        effective_sample_size = 0,
        geweke_stat = Inf,
        autocorrelation = 1,
        seed = seed,
        convergence_history = list(NULL)
      )
    )
  })
}

# Random Forest function -------
fit_rf_model <- function(df,
                         seed,
                         ntree_values = c(500, 750, 1000, 1500, 2000)) {
  set.seed(seed)
  # Prepare model data
  model_data <- prepare_model_data(df, seed = seed)
  
  # Calculate mtry values dynamically
  n_features <- ncol(model_data$X_train)
  fractions <- seq(0.1, 0.4, length.out = 5)
  mtry_values <- floor(n_features * fractions)
  mtry_values <- unique(mtry_values)
  
  cli::cli_inform("Using mtry values: {paste(mtry_values, collapse = ', ')}")
  
  # Initialize variables outside tryCatch
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  best_params <- NULL
  
  result <- tryCatch({
    total_combinations <- length(mtry_values) * length(ntree_values)
    current_combination <- 0
    
    # Grid search
    for (mtry in mtry_values) {
      for (ntree in ntree_values) {
        current_combination <- current_combination + 1
        
        cli::cli_inform(
          "Testing combination {current_combination}/{total_combinations}: mtry={mtry},
ntree={ntree}"
        )
        
        rf <- randomForest(
          x = model_data$X_train,
          y = model_data$y_train,
          mtry = mtry,
          ntree = ntree,
          importance = TRUE
        )
        
        predictions <- predict(rf, as.matrix(model_data$X_test))
        
        test_metrics <- calculate_metrics(
          model_data$y_test,
          predictions,
          n_predictors = ncol(model_data$X_train),
          model_data$y_train,
          predict(rf, as.matrix(model_data$X_train))
        )
        
        oob_error <- mean(rf$mse, na.rm = TRUE)
        
        current_metrics <- tibble(
          model = "RandomForest",
          test_metrics,
          mtry = mtry,
          ntree = ntree,
          oob_error = oob_error,
          seed = seed,
          variable_importance = list(as.data.frame(randomForest::importance(rf))),
          performance_curve = list(tibble(
            trees = seq_len(ntree), oob_error = rf$mse
          ))
        )
        
        if (test_metrics$train_rmse < best_score) {
          best_metrics <- current_metrics
          best_score <- test_metrics$train_rmse
          best_model <- rf
          best_params <- list(mtry = mtry, ntree = ntree)
          
          cli::cli_inform("New best model found: RMSE={round(best_score, 3)}")
        }
      }
    }
    
    # Return the best metrics found
    if (!is.null(best_metrics)) {
      best_metrics
    } else {
      cli::cli_warn("No valid model found")
      tibble(
        model = "RandomForest",
        mae = NA_real_,
        mse = NA_real_,
        r2 = NA_real_,
        rmse = NA_real_,
        nrmse = NA_real_,
        n_samples = 0,
        mtry = NA_integer_,
        ntree = NA_integer_,
        oob_error = NA_real_,
        seed = seed,
        variable_importance = list(NULL),
        performance_curve = list(NULL)
      )
    }
  }, error = function(e) {
    cli::cli_warn("Error in Random Forest fitting: {e$message}")
    tibble(
      model = "RandomForest",
      mae = NA_real_,
      mse = NA_real_,
      r2 = NA_real_,
      rmse = NA_real_,
      nrmse = NA_real_,
      n_samples = 0,
      mtry = NA_integer_,
      ntree = NA_integer_,
      oob_error = NA_real_,
      seed = seed,
      variable_importance = list(NULL),
      performance_curve = list(NULL)
    )
  })
  
  return(result)
}

# XGBoost Model Function -------
fit_xgboost_model <- function(df,
                              seed,
                              nrounds = c(100, 200, 300),
                              max_depth = c(3, 6, 9),
                              eta = c(0.01, 0.1, 0.3),
                              subsample = c(0.5, 0.75, 1.0)) {
  # Input validation
  if (!is.data.frame(df))
    stop("df must be a data frame")
  if (!is.numeric(seed))
    stop("seed must be numeric")
  
  set.seed(seed)
  
  result <- tibble()
  # Prepare model data
  model_data <- prepare_model_data(df, seed = seed)
  
  # Create parameter grid
  param_grid <- expand.grid(
    nrounds = nrounds,
    max_depth = max_depth,
    eta = eta,
    subsample = subsample
  )
  
  # Initialize tracking variables
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  best_params <- NULL
  total_combinations <- nrow(param_grid)
  
  cli::cli_inform("Starting grid search with {total_combinations} parameter combinations")
  
  # Grid search
  for (i in seq_len(nrow(param_grid))) {
    params <- list(
      max_depth = param_grid$max_depth[i],
      eta = param_grid$eta[i],
      subsample = param_grid$subsample[i],
      objective = "reg:squarederror"
    )
    
    cli::cli_inform(
      c("i" = "Testing combination {i}/{total_combinations}:", "*" = "max_depth={params$max_depth}, eta={params$eta}, subsample={params$subsample}, nrounds={param_grid$nrounds[i]}")
    )
    
    # Create DMatrix objects
    dtrain <- xgb.DMatrix(data = as.matrix(model_data$X_train),
                          label = model_data$y_train)
    dtest <- xgb.DMatrix(data = as.matrix(model_data$X_test),
                         label = model_data$y_test)
    
    # Train model
    watchlist <- list(train = dtrain, test = dtest)
    tryCatch({
    model <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = param_grid$nrounds[i],
      verbose = 0,
      watchlist = watchlist,
      early_stopping_rounds = 10
    )}, error = function(e) {
      cli::cli_warn("Couldn't train the model with the specified parameters, skipping iteration. Error: {e$message}")
      next
    })
    

    # Make predictions
    test_pred <- predict(model, dtest)
    train_pred <- predict(model, dtrain)
    
    # Calculate metrics
    current_metrics <- calculate_metrics(
      test_y_true = model_data$y_test,
      test_y_pred = test_pred,
      n_predictors = ncol(model_data$X_train),
      train_y_true = model_data$y_train,
      train_y_pred = train_pred
    )
    
    if (current_metrics$train_rmse < best_score) {
      best_metrics <- current_metrics
      best_score <- current_metrics$train_rmse
      best_model <- model
      best_params <- param_grid[i, ]
      
      cli::cli_inform("New best model found: RMSE={round(best_score, 4)}")
    }
  }
  
  # Prepare final results
  result <- tibble(
    model = "XGBoost",
    best_metrics,
    best_nrounds = best_params$nrounds,
    best_max_depth = best_params$max_depth,
    best_eta = best_params$eta,
    best_subsample = best_params$subsample,
    seed = seed,
    variable_importance = list(xgb.importance(model = best_model))
  )
  return(result)
}


# SVM Model Function -------
fit_svm_model <- function(df,
                          seed,
                          kernel = "rbfdot",
                          C = c(0.1, 1, 10),
                          sigma = c(0.1, 0.5, 1)) {
  # Input validation
  if (!is.data.frame(df))
    stop("df must be a data frame")
  if (!is.numeric(seed))
    stop("seed must be numeric")
  if (!all(is.numeric(C)))
    stop("C values must be numeric")
  if (!all(is.numeric(sigma)))
    stop("sigma values must be numeric")
  if (!kernel %in% c(
    "rbfdot",
    "polydot",
    "vanilladot",
    "tanhdot",
    "laplacedot",
    "besseldot",
    "anovadot"
  )) {
    stop("Invalid kernel type")
  }
  
  set.seed(seed)
  
  result <-tibble()
    # Prepare model data
    model_data <- prepare_model_data(df, seed = seed)
  
  # Initialize tracking variables
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  best_params <- NULL
  
  # Calculate total combinations for progress tracking
  total_combinations <- length(C) * length(sigma)
  current_combination <- 0
  
  cli::cli_inform("Starting grid search with {total_combinations} parameter combinations")
  
  # Grid search
  for (c_val in C) {
    for (sigma_val in sigma) {
      current_combination <- current_combination + 1
      
      cli::cli_inform(
        c("i" = "Testing combination {current_combination}/{total_combinations}:", "*" = "C={c_val}, sigma={sigma_val}, kernel={kernel}")
      )
      
      # Train model
      svm_model <- ksvm(
        x = as.matrix(model_data$X_train),
        y = model_data$y_train,
        kernel = kernel,
        C = c_val,
        sigma = sigma_val,
        scaled = FALSE,
        type = "eps-svr"  # Explicitly set regression type
      )
      
      # Make predictions
      test_pred <- predict(svm_model, as.matrix(model_data$X_test))
      train_pred <- predict(svm_model, as.matrix(model_data$X_train))
      
      # Calculate metrics with updated function
      current_metrics <- calculate_metrics(
        test_y_true = model_data$y_test,
        test_y_pred = test_pred,
        n_predictors = ncol(model_data$X_train),
        train_y_true = model_data$y_train,
        train_y_pred = train_pred
      )
      
      if (current_metrics$train_rmse < best_score) {
        best_metrics <- current_metrics
        best_score <- current_metrics$train_rmse
        best_model <- svm_model
        best_params <- list(C = c_val, sigma = sigma_val)
        
        cli::cli_inform("New best model found: RMSE={round(best_score, 4)}")
      }
    }
  }
  
  # Prepare final results
  if (!is.null(best_metrics)) {
    result <- tibble(
      model = "SVM",
      best_metrics,
      kernel = kernel,
      best_C = best_params$C,
      best_sigma = best_params$sigma,
      seed = seed#,
      #model_object = list(best_model),  # Store the model object if needed
      #performance_curve = list(NULL)  # Could add learning curve if needed
    )
  }
  
  return(result)
}


# Neural Network Model Function -------
fit_nnet_model <- function(df,
                           hidden_units = c(5, 10, 15),
                           weight_decay = c(0.01, 0.1, 2.5),
                           maxit = 2000,
                           var_threshold = 0.95,
                           pca_threshold = 50,
                           k_folds = 5,
                           seed) {
  # Input validation
  if (!is.data.frame(df))
    stop("df must be a data frame")
  if (!is.numeric(seed))
    stop("seed must be numeric")
  
  
  set.seed(seed)
  result <- tibble()
  # Prepare model data - data is already scaled here
  model_data <- prepare_model_data(df, seed = seed)
  
  # Initialize tracking variables
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  best_params <- NULL
  cv_results <- tibble()
  
  # Feature processing
  feature_processing_result <- process_features(model_data, pca_threshold, var_threshold)
  
  X_train_final <- feature_processing_result$X_train
  X_test_final <- feature_processing_result$X_test
  feature_info <- feature_processing_result$info
  
  # Create cross-validation folds
  folds <- createFolds(model_data$y_train, k = k_folds, list = TRUE)
  
  total_combinations <- length(hidden_units) * length(weight_decay)
  current_combination <- 0
  
  cli::cli_inform("Starting grid search with {total_combinations} parameter combinations")
  
  # Grid search with cross-validation
  for (units in hidden_units) {
    for (decay in weight_decay) {
      current_combination <- current_combination + 1
      
      cli::cli_inform(
        c("i" = "Testing combination {current_combination}/{total_combinations}:", "*" = "hidden_units={units}, decay={decay}")
      )
      
      cv_result <- perform_cv(X_train_final,
                              model_data$y_train,
                              folds,
                              units,
                              decay,
                              maxit)
      
      if (!cv_result$valid) {
        cli::cli_warn(cv_result$message)
        next
      }
      
      # Update best model if necessary
      if (cv_result$rmse < best_score) {
        best_model_result <- fit_best_model(
          X_train_final,
          model_data$y_train,
          X_test_final,
          model_data$y_test,
          units,
          decay,
          maxit
        )
        
        if (best_model_result$valid) {
          best_metrics <- best_model_result$metrics
          best_score <- cv_result$rmse
          best_model <- best_model_result$model
          best_params <- list(
            hidden_units = units,
            decay = decay,
            n_predictors = ncol(X_train_final),
            using_pca = feature_info$using_pca,
            cv_score = cv_result$rmse
          )
          
          if (feature_info$using_pca) {
            best_params$n_components <- feature_info$n_components
            best_params$selected_variables <- feature_info$selected_vars
          }
          
          cli::cli_inform("New best model found: CV RMSE={round(best_score, 4)}")
        }
      }
      
      # Store CV results
      cv_results <- bind_rows(
        cv_results,
        tibble(
          hidden_units = units,
          decay = decay,
          cv_r2 = cv_result$mean_cv_score,
          cv_rmse = cv_result$rmse,
          cv_mae = cv_result$mae
        )
      )
    }
  }
  
  # Prepare final results
  if (!is.null(best_metrics)) {
    result <- tibble(
      model = "NeuralNetwork",
      best_metrics,
      best_hidden_units = best_params$hidden_units,
      best_decay = best_params$decay,
      n_predictors = best_params$n_predictors,
      using_pca = best_params$using_pca,
      n_components_used = if (best_params$using_pca)
        best_params$n_components
      else
        NA,
      selected_variables = if (best_params$using_pca)
        list(best_params$selected_variables)
      else
        NA,
      cv_results = list(cv_results),
      seed = seed
    )
  }
  
  return(result)
}

process_features <- function(model_data, pca_threshold, var_threshold) {
  n_predictors <- ncol(model_data$X_train)
  cli::cli_inform("Number of predictors: {n_predictors}")
  
  if (n_predictors > pca_threshold) {
    cli::cli_inform("Applying PCA dimension reduction...")
    pca_result <- prcomp(model_data$X_train, center = FALSE, scale. = FALSE)  # Data already scaled
    var_explained <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
    n_keep <- min(which(var_explained >= var_threshold)[1], 10)
    
    important_loadings <- abs(pca_result$rotation[,1:n_keep]) #%>% as.data.frame()
    var_importance <- important_loadings
    if (!is.null(dim(important_loadings))) {
      var_importance <- rowMeans(important_loadings)
    } 
    important_vars <- which(var_importance > mean(var_importance))
    
    X_train_final <- model_data$X_train[, important_vars]
    X_test_final <- model_data$X_test[, important_vars]
    
    info <- list(
      using_pca = TRUE,
      n_components = n_keep,
      selected_vars = colnames(model_data$X_train)[important_vars],
      var_explained = var_explained[n_keep]
    )
    return(list(
      X_train = X_train_final,
      X_test = X_test_final,
      info = info
    ))
    }
    
    cli::cli_inform("Using original variables...")
    X_train_final <- model_data$X_train
    X_test_final <- model_data$X_test
    
    info <- list(
      using_pca = FALSE,
      n_components = NA,
      selected_vars = NA,
      var_explained = NA
    )

  
  return(list(
    X_train = X_train_final,
    X_test = X_test_final,
    info = info
  ))
}

perform_cv <- function(X_train,
                       y_train,
                       folds,
                       units,
                       decay,
                       maxit) {
  cv_scores <- numeric(length(folds))
  cv_metrics <- list()
  
  for (fold in seq_along(folds)) {
    train_idx <- unlist(folds[-fold])
    valid_idx <- folds[[fold]]
    
    tryCatch({
      model <- nnet(
        x = X_train[train_idx, ],
        y = y_train[train_idx],
        size = units,
        decay = decay,
        maxit = maxit,
        linout = TRUE,
        trace = FALSE,
        MaxNWts = 10000
      )
      
      pred_cv <- predict(model, X_train[valid_idx, ])
      
      cv_metrics[[fold]] <- calculate_metrics(
        test_y_true = y_train[valid_idx],
        test_y_pred = pred_cv,
        n_predictors = ncol(X_train),
        train_y_true = y_train[train_idx],
        train_y_pred = predict(model, X_train[train_idx, ])
      )
      
      cv_scores[fold] <- cv_metrics[[fold]]$train_rmse
      
    }, error = function(e) {
      cv_scores[fold] <- NA
      cli::cli_warn("Error in fold {fold}: {e$message}")
    })
  }
  
  mean_cv_score <- mean(cv_scores, na.rm = TRUE)
  
  if (is.na(mean_cv_score) || mean_cv_score < 0) {
    return(
      list(
        valid = FALSE,
        message = "Invalid CV scores",
        mean_cv_score = NA,
        rmse = NA,
        mae = NA
      )
    )
  }
  
  return(list(
    valid = TRUE,
    mean_cv_score = mean_cv_score,
    rmse = mean(sapply(cv_metrics, function(x) x$test_rmse), na.rm = TRUE),
    mae = mean(sapply(cv_metrics, function(x) x$test_mae), na.rm = TRUE)
  ))
}

fit_best_model <- function(X_train,
                           y_train,
                           X_test,
                           y_test,
                           units,
                           decay,
                           maxit) {
  tryCatch({
    model <- nnet(
      x = X_train,
      y = y_train,
      size = units,
      decay = decay,
      maxit = maxit,
      linout = TRUE,
      trace = FALSE,
      MaxNWts = 10000
    )
    
    test_pred <- predict(model, X_test)
    train_pred <- predict(model, X_train)
    
    metrics <- calculate_metrics(
      test_y_true = y_test,
      test_y_pred = test_pred,
      n_predictors = ncol(X_train),
      train_y_true = y_train,
      train_y_pred = train_pred
    )
    
    return(list(
      valid = TRUE,
      metrics = metrics,
      model = model
    ))
    
  }, error = function(e) {
    cli::cli_warn("Error fitting best model: {e$message}")
    return(list(
      valid = FALSE,
      metrics = NULL,
      model = NULL
    ))
  })
}

# Ranger Model Function -------
fit_ranger_model <- function(df,
                             seed,
                             num_trees = c(500, 1000, 1500, 2000),
                             min_node_size = c(5, 10, 15, 20)) {
  # Input validation
  if (!is.data.frame(df))
    stop("df must be a data frame")
  if (!is.numeric(seed))
    stop("seed must be numeric")
  
  
  set.seed(seed)
  
  result <- tibble()
    # Prepare model data
    model_data <- prepare_model_data(df, seed = seed)
  
  # Calculate mtry values
  num_variables <- ncol(model_data$X_train)
  fractions <- seq(0.1, 0.4, length.out = 5)
  mtry_values <- unique(floor(num_variables * fractions))
  
  cli::cli_inform("Using mtry values: {paste(mtry_values, collapse = ', ')}")
  
  # Initialize tracking variables
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  best_params <- NULL
  
  total_combinations <- length(num_trees) * length(mtry_values) * length(min_node_size)
  current_combination <- 0
  
  cli::cli_inform("Starting grid search with {total_combinations} parameter combinations")
  
  # Grid search
  for (trees in num_trees) {
    for (mtry in mtry_values) {
      for (node_size in min_node_size) {
        current_combination <- current_combination + 1
        
        cli::cli_inform(
          c("i" = "Testing combination {current_combination}/{total_combinations}:", "*" = "trees={trees}, mtry={mtry}, min_node_size={node_size}")
        )
        
        tryCatch({
          ranger_model <- ranger(
            x = model_data$X_train,
            y = model_data$y_train,
            num.trees = trees,
            mtry = mtry,
            min.node.size = node_size,
            importance = 'impurity'
          )
          
          # Make predictions
          test_pred <- predict(ranger_model, data = model_data$X_test)$predictions
          train_pred <- predict(ranger_model, data = model_data$X_train)$predictions
          
          # Calculate metrics with updated function
          current_metrics <- calculate_metrics(
            test_y_true = model_data$y_test,
            test_y_pred = test_pred,
            n_predictors = ncol(model_data$X_train),
            train_y_true = model_data$y_train,
            train_y_pred = train_pred
          )
          
          if (current_metrics$train_rmse < best_score) {
            best_metrics <- current_metrics
            best_score <- current_metrics$train_rmse
            best_model <- ranger_model
            best_params <- list(
              num_trees = trees,
              mtry = mtry,
              min_node_size = node_size
            )
            
            cli::cli_inform("New best model found: RMSE={round(best_score, 4)}")
          }
          
        }, error = function(e) {
          cli::cli_warn("Error in model fitting: {e$message}")
        })
      }
    }
  }
  
  # Prepare final results
  if (!is.null(best_metrics)) {
    # Get variable importance from best model
    var_imp <- importance(best_model)
    var_imp_df <- data.frame(variable = names(var_imp),
                             importance = as.numeric(var_imp)) %>%
      arrange(desc(importance))
    
    result <- tibble(
      model = "Ranger",
      best_metrics,
      best_num_trees = best_params$num_trees,
      best_mtry = best_params$mtry,
      best_min_node_size = best_params$min_node_size,
      oob_error = best_model$prediction.error,
      variable_importance = list(var_imp_df),
      seed = seed
    )
  }
  return(result)
}


# Partial Least Squares Regression Function -------
fit_plsr_model <- function(df,
                           seed,
                           ncomp_values = c(5, 10, 15, 20, 25),
                           validation = "CV",
                           segments = 5) {
  # Input validation
  if (!is.data.frame(df))
    stop("df must be a data frame")
  if (!is.numeric(seed))
    stop("seed must be numeric")
  if (!validation %in% c("CV", "LOO"))
    stop("validation must be 'CV' or 'LOO'")
  
  set.seed(seed)
  
  result <- tibble()
    # Prepare model data
  model_data <- prepare_model_data(df, seed = seed)
  
  # Convert data to format required by pls
  df_train <- data.frame(y = model_data$y_train, model_data$X_train)
  df_test <- data.frame(model_data$X_test)
  
  # Initialize tracking variables
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  best_ncomp <- NULL
  validation_metrics <- list()
  
  # Calculate maximum possible components
  max_possible_ncomp <- min(ncol(model_data$X_train) - 1, nrow(model_data$X_train) - 1)
  
  # Filter valid ncomp values
  valid_ncomp_values <- ncomp_values[ncomp_values <= max_possible_ncomp]
  
  cli::cli_inform("Testing {length(valid_ncomp_values)} component values")
  
  # Grid search
  for (ncomp in valid_ncomp_values) {
    cli::cli_inform(c("i" = "Testing with {ncomp} components", "*" = "Validation method: {validation}"))
    
    tryCatch({
      # Fit model
      plsr_model <- plsr(
        y ~ .,
        data = df_train,
        ncomp = ncomp,
        scale = FALSE,
        # Already scaled in prepare_model_data
        validation = validation,
        segments = segments
      )
      
      # Make predictions
      test_pred <- predict(plsr_model, newdata = df_test, ncomp = ncomp)
      train_pred <- predict(plsr_model, ncomp = ncomp)
      
      # Calculate metrics
      current_metrics <- calculate_metrics(
        test_y_true = model_data$y_test,
        test_y_pred = test_pred,
        n_predictors = ncomp,
        # Use number of components as n_predictors
        train_y_true = model_data$y_train,
        train_y_pred = train_pred
      )
      
      # Store validation metrics
      validation_metrics[[length(validation_metrics) + 1]] <- list(
        ncomp = ncomp,
        metrics = current_metrics,
        cv_press = RMSEP(plsr_model)$val[1, 1, ncomp]
      )
      
      if (current_metrics$train_rmse < best_score) {
        best_metrics <- current_metrics
        best_score <- current_metrics$train_rmse
        best_model <- plsr_model
        best_ncomp <- ncomp
        
        cli::cli_inform("New best model found: RMSE={round(best_score, 4)}")
      }
      
    }, error = function(e) {
      cli::cli_warn("Error fitting model with {ncomp} components: {e$message}")
    })
  }
  
  # Prepare final results
  if (!is.null(best_metrics)) {
    # Get variable importance (loadings from best model)
    loadings <- loadings(best_model)[, 1:best_ncomp, drop = FALSE]
    var_imp <- apply(abs(loadings), 1, mean)
    var_imp_df <- data.frame(variable = names(var_imp),
                             importance = as.numeric(var_imp)) %>%
      arrange(desc(importance))
    
    # Calculate explained variance
    explained_var <- explvar(best_model)[1:best_ncomp]
    
    result <- tibble(
      model = "PLSR",
      best_metrics,
      optimal_ncomp = best_ncomp,
      validation_method = validation,
      explained_variance = list(explained_var),
      variable_importance = list(var_imp_df),
      validation_summary = list(validation_metrics),
      seed = seed
    )
  }
  
  return(result)
}


# Linear Support Vector Regression Function -------
fit_linear_svr_model <- function(df,
                                 seed,
                                 cost_values = c(0.1, 1, 10, 100),
                                 epsilon_values = c(0.01, 0.1, 0.5)) {
  # Input validation
  if (!is.data.frame(df))
    stop("df must be a data frame")
  if (!is.numeric(seed))
    stop("seed must be numeric")
  
  set.seed(seed)
  model_data <- prepare_model_data(df)
  result <- tibble()
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  
  # Calculate number of predictors
  n_predictors <- ncol(model_data$X_train)
  
  # Grid search
  for (cost in cost_values) {
    for (epsilon in epsilon_values) {
      tryCatch({
        # Type=11 for L2-regularized L2-loss support vector regression (primal)
        svr_model <- LiblineaR(
          data = model_data$X_train,
          target = model_data$y_train,
          type = 11,
          cost = cost,
          epsilon = epsilon,
          svr_eps = 0.1,
          verbose = FALSE
        )
        
        # Generate predictions for both test and training sets
        test_predictions <- predict(svr_model, model_data$X_test)$predictions
        train_predictions <- predict(svr_model, model_data$X_train)$predictions
        
        # Calculate metrics using updated function
        current_metrics <- calculate_metrics(
          test_y_true = model_data$y_test,
          test_y_pred = test_predictions,
          n_predictors = n_predictors,
          train_y_true = model_data$y_train,
          train_y_pred = train_predictions
        )
        
        if (current_metrics$train_rmse < best_score) {
          best_metrics <- current_metrics
          best_score <- current_metrics$r2
          best_params <- list(cost = cost, epsilon = epsilon)
          best_model <- svr_model
          
          cli::cli_inform("New best model found: R2={round(best_score, 4)}")
        }
      }, error = function(e) {
        cli::cli_warn("Error fitting model with cost={cost}, epsilon={epsilon}: {e$message}")
      })
    }
  }
  
  # Prepare final results
  if (!is.null(best_metrics)) {
    result <- tibble(
      model = "LinearSVR",
      best_metrics,
      best_cost = best_params$cost,
      best_epsilon = best_params$epsilon,
      seed = seed
    )
  }
  
  return(result)
}

# Elastic Net Regression Function -------
fit_elastic_net_model <- function(df,
                                  seed,
                                  alpha_values = c(1:9 / 10),
                                  n_folds = 10) {
  # Input validation
  if (!is.data.frame(df))
    stop("df must be a data frame")
  if (!is.numeric(seed))
    stop("seed must be numeric")
  
  set.seed(seed)
  model_data <- prepare_model_data(df, seed = seed)
  result <- tibble()
  best_metrics <- NULL
  best_score <- Inf
  best_model <- NULL
  
  # Get number of predictors
  n_predictors <- ncol(model_data$X_train)
  
  # Custom lambda sequence - make it wider and include smaller values
  lambda_seq <- exp(seq(log(0.0001), log(20), length.out = 100))
  
  # Grid search over alpha values
  cli::cli_inform("Testing {length(alpha_values)} alpha values for Elastic Net")
  
  for (alpha in alpha_values) {
    cli::cli_inform(c("i" = "Testing alpha = {alpha}"))
    
    tryCatch({
      # Fit model with cross-validation for lambda selection
      cv_fit <- cv.glmnet(
        x = model_data$X_train,
        y = model_data$y_train,
        alpha = alpha,
        nfolds = n_folds,
        lambda = lambda_seq
      )
      
      # Only proceed if we have non-zero coefficients
      n_nonzero <- as.numeric(cv_fit$nzero[which(cv_fit$lambda == cv_fit$lambda.min)])
      
      if (n_nonzero > 0) {
        # Get predictions for both test and training sets
        test_predictions <- predict(cv_fit,
                                    newx = model_data$X_test,
                                    s = cv_fit$lambda.min)
        train_predictions <- predict(cv_fit,
                                     newx = model_data$X_train,
                                     s = cv_fit$lambda.min)
        
        # Calculate metrics using updated function
        current_metrics <- calculate_metrics(
          test_y_true = model_data$y_test,
          test_y_pred = test_predictions,
          n_predictors = n_nonzero,
          # Use number of non-zero coefficients as effective predictors
          train_y_true = model_data$y_train,
          train_y_pred = train_predictions
        )
        
        if (!is.na(current_metrics$train_rmse)) {
          if (current_metrics$rmse < best_score) {
            best_metrics <- current_metrics
            best_score <- current_metrics$train_rmse
            best_model <- cv_fit
            best_params <- list(
              alpha = alpha,
              lambda = cv_fit$lambda.min,
              n_nonzero = n_nonzero
            )
            
            cli::cli_inform("New best model found: RMSE={round(best_score, 4)}")
          }
        }
      }
    }, error = function(e) {
      cli::cli_warn("Error fitting model with alpha={alpha}: {e$message}")
    })
  }
  
  # Ridge regression part (alpha = 0)
  cli::cli_inform("Fitting Ridge Regression model (alpha = 0)")
  
  
  cv_fit <- cv.glmnet(
    x = model_data$X_train,
    y = model_data$y_train,
    alpha = 0,
    nfolds = n_folds,
    lambda = lambda_seq
  )
  
  # Get predictions for both test and training sets
  test_predictions <- predict(cv_fit, newx = model_data$X_test, s = cv_fit$lambda.min)
  train_predictions <- predict(cv_fit, newx = model_data$X_train, s = cv_fit$lambda.min)
  
  # Calculate metrics for Ridge regression
  ridge_metrics <- calculate_metrics(
    test_y_true = model_data$y_test,
    test_y_pred = test_predictions,
    n_predictors = n_predictors,
    # Use all predictors for Ridge regression
    train_y_true = model_data$y_train,
    train_y_pred = train_predictions
  )
  
  # Combine results
  result <- bind_rows(
    if (!is.null(best_metrics)) {
      best_metrics %>%
        mutate(
          model = "ElasticNet",
          best_alpha = best_params$alpha,
          best_lambda = best_params$lambda,
          n_nonzero_coefficients = best_params$n_nonzero,
          seed = seed
        )
    },
    ridge_metrics %>%
      mutate(
        model = "RidgeRegression",
        best_alpha = 0,
        best_lambda = cv_fit$lambda.min,
        n_nonzero_coefficients = n_predictors,
        seed = seed
      )
  )
  
  return(result)
}


# MODEL REGISTRY -------
model_registry <- list(
  bayes = fit_bayes_model,
  rf = fit_rf_model,
  xgboost = fit_xgboost_model,
  svm = fit_svm_model,
  nnet = fit_nnet_model,
  ranger = fit_ranger_model,
  plsr = fit_plsr_model,
  linear_svr = fit_linear_svr_model,
  elastic_net = fit_elastic_net_model
)
