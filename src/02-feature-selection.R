# ===================================================================
# 02-feature-selection.R
# Feature selection algorithms (abbreviated for space)
# ===================================================================

#`` Include all feature selection functions from Script 2
# B-coefficient filtering from PLSR ----
#' @param df Input dataframe with GTR as first column, metadata (ID, Type), then spectral data
#' @param threshold Coefficient threshold for selection
#' @return List containing selected wavelengths and their indices
plsr_bcoef_filter <- function(df, threshold_quantile = 0.95, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract response and spectral data
  y <- colnames(df[, 1])
  X <- df[, (4:ncol(df))]
  plsr_data <- df[, c(1,4:ncol(df))] %>% rename("y"=all_of(y))
  
  # Fit PLSR model
  plsr_model <- plsr(y ~ .,
                     data = plsr_data,
                     scale = TRUE,
                     validation = "CV")
  
  # Coefficients Threshold
  coefs <- coef(plsr_model)
  abs_coefs <- abs(coef(plsr_model))
  coef_threshold <- quantile(abs_coefs, threshold_quantile)
  
  # Select vars that exceed threshold
  selected_idx <- which(abs_coefs > coef_threshold)
  selected_wavelengths <- colnames(X)[selected_idx]
  
  # Create summary information
  summary_info <- data.frame(
    wavelength = selected_wavelengths,
    coefficient = coefs[selected_idx],
    abs_coefficient = abs_coefs[selected_idx]
  ) %>%
    arrange(desc(abs_coefficient))
  
  return(
    list(
      indices = selected_idx,
      wavelengths = selected_wavelengths,
      coefficients = coefs[selected_idx],
      threshold_used = coef_threshold,
      summary = summary_info
    )
  )
}

# Successive Projections Algorithm (SPA) ----
#' @param df Data frame with:
#'          Column 1: Response variable
#'          Columns 2-3: Metadata
#'          Columns 4-465: Spectral data (462 wavelengths)
#' @param n_vars Number of wavelengths to select (default = 10)
#' @param initial_wave Starting wavelength index (optional)
#' @return Vector of selected wavelength indices
#' https://www.sciencedirect.com/science/article/pii/S0169743901001198

# SPA function that takes initial wavelength as parameter
spa_select <- function(df, n_vars = 10, initial_wave = NULL, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Extract spectral matrix and response matrix
  X <- as.matrix(df[, 4:ncol(df)])  # Spectral cols
  Y <- as.matrix(df[, 1])         # Response col
  
  n <- ncol(X)
  selected <- integer(n_vars)
  
  # If no initial wavelength provided, use variance method as fallback
  if(is.null(initial_wave)) {
    var_power <- apply(X, 2, var)
    selected[1] <- which.max(var_power)
  } else {
    selected[1] <- initial_wave
  }
  
  # Main SPA loop
  for(k in 2:n_vars) {
    proj <- X
    
    # Project data onto orthogonal space of selected variables
    for(j in 1:(k-1)) {
      v <- X[, selected[j]]
      v_norm <- sqrt(sum(v^2))
      if(v_norm > 1e-10) {
        v_normalized <- v / v_norm
        proj <- proj - tcrossprod(v_normalized) %*% proj
      }
    }
    
    # Select wavelength with maximum projection
    proj_norms <- apply(proj, 2, function(x) sqrt(sum(x^2)))
    proj_norms[selected[1:(k-1)]] <- -Inf
    selected[k] <- which.max(proj_norms)
  }
  
  return(selected)
}

# Function to calculate RMSEP
calculate_rmsep <- function(X_cal, Y_cal, X_test, Y_test, wavelengths) {
  # Build MLR model using selected wavelengths
  model_data <- as.data.frame(bind_cols(Y_cal, X_cal[, wavelengths]))
  formula_str <- paste("model_data[,1] ~",
                       paste("`", colnames(model_data[,-1]), "`",
                             sep="", collapse = " + "))      
  model <- lm(as.formula(formula_str), data = model_data)
  
  # Predict test set
  X_test_selected <- X_test[, wavelengths]
  Y_pred <- predict(model, newdata = as.data.frame(X_test_selected))
  
  # Calculate RMSEP
  rmsep <- sqrt(mean((Y_test - Y_pred)^2))
  return(rmsep)
}

# Function to find optimal initial wavelength
find_optimal_start <- function(df_cal, df_test, n_vars) {
  # Extract matrices
  X_cal <- as.matrix(df_cal[, 4:ncol(df_cal)])
  Y_cal <- as.matrix(df_cal[, 1])
  X_test <- as.matrix(df_test[, 4:ncol(df_test)])
  Y_test <- as.matrix(df_test[, 1])
  
  n_wavelengths <- ncol(X_cal)
  rmsep_values <- numeric(n_wavelengths)
  
  # Try each wavelength as starting point
  for(i in 1:n_wavelengths) {
    # Run SPA with this starting point
    selected_waves <- spa_select(df_cal, n_vars, initial_wave = i)
    
    cor_matrix <- cor(X_cal[,selected_waves])
    # Get high correlations (absolute value > threshold)
    high_cors <- abs(cor_matrix) > 0.80 & upper.tri(cor_matrix)
    
    # If any high correlations found
    if(any(high_cors)) {
      # Add penalty factor based on number/magnitude of high correlations
      n_high_cors <- sum(high_cors)
      # should exclude these from multiple linear regression
      rmsep_values[i] <- calculate_rmsep(X_cal, Y_cal, X_test, Y_test, selected_waves) * (1 + 0.1 * n_high_cors)
    } else {
      # Calculate RMSEP
      rmsep_values[i] <- calculate_rmsep(X_cal, Y_cal, X_test, Y_test, selected_waves)
    }
  }
  
  # Return wavelength with minimum RMSEP
  optimal_start <- which.min(rmsep_values)
  return(list(
    optimal_wavelength = optimal_start,
    rmsep = rmsep_values[optimal_start]
  ))
}

# Main function to run optimized SPA
run_optimized_spa <- function(df, n_vars = 10, test_size = 0.2, seed = NULL) {
  
  # Get indices for train-test split
  n_samples <- nrow(df)
  
  train_idx <- createDataPartition(y = pull(df[,1]), 
                                   p = 1 - test_size,
                                   list = FALSE)
  # Split data
  df_cal <- df[train_idx, ]
  df_test <- df[-train_idx, ]
  # Find optimal starting wavelength
  optimal_start <- find_optimal_start(df_cal, df_test, n_vars)
  
  # Run SPA with optimal starting point
  selected_waves <- spa_select(df_cal, n_vars, optimal_start$optimal_wavelength)
  
  return(list(
    indices = selected_waves,
    wavelengths = colnames(df[, selected_waves + 3]),
    initial_wavelength = optimal_start$optimal_wavelength,
    rmsep = optimal_start$rmsep
  ))
}

# Genetic Algorithm Feature Selection ----
#' @param df Input dataframe with GTR as first column, metadata (ID, Type), then spectral data
#' @return List containing selected wavelengths and their indices
ga_select <- function(df, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Extract response and spectral data
  y <- df[, 1] # response var
  X <- df[, 4:ncol(df)]  # Spectral data starts at column 4
  
  # Define fitness function
  fitness <- function(x) {
    n_selected <- sum(x)
    
    # Base cases
    if (n_selected == 0)
      return(-1e6)
    if (n_selected > ncol(X) / 2)
      return(-1e6)  # Prevent selecting too many features
    
    selected <- which(x == 1)
    X_selected <- as.matrix(X[, selected])
    
    tryCatch({
      # Cross-validated model fit
      cv_fit <- cv.glmnet(
        X_selected,
        y,
        alpha = 0.5,
        #standardize = TRUE,
        nfolds = 5
      )
      
      # Calculate base performance
      mse <- min(cv_fit$cvm)
      rmse <- sqrt(mse)
      
      # Complexity penalty based on number of features
      complexity_penalty <- log(n_selected) / log(ncol(X))
      
      # Final score balancing performance and complexity
      score <- -rmse * (1 + complexity_penalty)
      
      return(score)
    }, error = function(e) {
      return(-1e6)
    })
  }
  
  # Initialize population with varying feature counts
  initial_pop <- matrix(0, nrow = 50, ncol = ncol(X))
  
  # Initialize with varying feature counts
  for (i in 1:50) {
    n_sel <- sample(2:ceiling(ncol(X) / 4), 1)  # Random selection between 2 and 25% of features
    selected <- sample(ncol(X), size = n_sel)
    initial_pop[i, selected] <- 1
  }
  
  # Run GA
  ga_result <- ga(
    type = "binary",
    # For binary feature selection
    fitness = fitness,
    # Our defined fitness function
    nBits = ncol(X),
    # One bit per wavelength
    names = colnames(X),
    # Names for wavelengths
    popSize = 50,
    # Default population size
    pcrossover = 0.8,
    # Default crossover probability
    pmutation = 0.1,
    # Default mutation probability
    elitism = 3,
    # Keep top 3 solutions
    maxiter = 100,
    # Maximum iterations
    run = 50,
    # Convergence criterion
    suggestions = initial_pop,
    parallel = FALSE,
    # Sequential processing
    monitor = gaMonitor      # Default monitoring
  )
  
  # Extract results
  selected_idx <- which(ga_result@solution[1, ] == 1)
  selected_wavelengths <- colnames(X)[selected_idx]
  
  # Calculate selection frequencies
  selection_frequency <- colMeans(ga_result@population)
  
  # Create summary information
  summary_info <- data.frame(
    wavelength = selected_wavelengths,
    selection_frequency = selection_frequency[selected_idx],
    index = selected_idx
  ) %>%
    arrange(desc(selection_frequency))
  
  return(
    list(
      indices = selected_idx,
      wavelengths = selected_wavelengths,
      fitness_value = ga_result@fitnessValue,
      summary = summary_info,
      ga_object = ga_result,
      n_selected = length(selected_idx)
    )
  )
}

# Uninformative Variable Elimination ----
#' @param df Input dataframe with response as first column
#' @param threshold Reliability threshold
#' @return List containing selected wavelengths and their indices
uve_select <- function(df, threshold = 0.95, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  X <- as.matrix(df[, 4:ncol(df)])  # Spectral data
  y <- as.matrix(df[, 1])
  
  # Calculate regression coefficients for each variable
  coefs <- apply(X, 2, function(x) {
    model <- lm(y ~ x)
    return(coef(model)[2])  # Get slope coefficient
  })
  
  # Calculate reliability score
  reliability <- abs(coefs) / sd(coefs)
  
  # Select variables above threshold
  selected_idx <- which(reliability > quantile(reliability, threshold))
  selected_wavelengths <- colnames(X)[selected_idx]
  
  return(
    list(
      indices = selected_idx,
      wavelengths = selected_wavelengths,
      reliability_scores = reliability[selected_idx]
    )
  )
}

# Competitive Adaptive Re-weighted Sampling ---------
#' @param df Input dataframe with GTR as first column
#' @param n_samples Number of sampling runs
#' @param nLV Number of latent variables for PLS
#' @param fold Number of CV folds
#' @param partition_type Cross validation partition type ("random", "consecutive", "interleaved")
#' @return List containing selected wavelengths and their indices
#' @references https://github.com/lhdcsu/libPLS

cars_select <- function(df,
                        # nLV = 2,
                        fold = 10,
                        partition_type = "interleaved",
                        seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  all_seeds <- tibble()
  seed_list <- sample(1:1e+06, 10)
 
  # Seed loop
  for (current_seed in seed_list) {
    set.seed(current_seed)
    
    # Data preparation
    X <- as.matrix(df[, 4:ncol(df)])
    y <- as.matrix(df[, 1])
    n_samples <- nrow(X)
    
    # Sort data by response variable
    new_order <- order(y)
    X <- X[new_order, ]
    y <- y[new_order]
    
    # Initialize storage
    rmsecv <- numeric(n_samples)
    num_lv <- numeric(n_samples)
    coef_matrix <- matrix(0, ncol(X), n_samples)
    subset_vars <- seq_len(ncol(X))
    
    # Parameters for exponential decay
    ratio0 <- 1
    ratio1 <- 2 / ncol(X)
    b <- log(ratio0 / ratio1) / (n_samples - 1)
    a <- ratio0 * exp(b)
    
    # Initialize weights
    weights <- rep(1, ncol(X))
    
    # Main loop
    for (sample_idx in seq_len(n_samples)) {
      x_temp <- X[, subset_vars, drop = FALSE]
      data_temp <- tibble(y = y) %>% bind_cols(as_tibble(x_temp))
      
      # n_comp <- min(nLV, ncol(x_temp))
      
      # Fit PLS model
      model <- mvr(
        y ~ .,
        ncomp = 2,
        data = data_temp,
        method = "simpls",
        scale = T
      )
      
      # Cross validation
      cv_result <- pls::crossval(model, segments = fold, segment.type = partition_type)
      
      # Calculate RMSECV
      press <- cv_result$validation$PRESS
      rmsecv[sample_idx] <- min(sqrt(press / nrow(X)))
      num_lv[sample_idx] <- which.min(sqrt(press / nrow(X)))
      
      # Update coefficients
      coef_temp <- numeric(ncol(X))
      coef_iter <- model$coefficients[, , num_lv[sample_idx]]
      coef_temp[subset_vars] <- coef_iter
      coef_matrix[, sample_idx] <- coef_temp
      
      # Update weights
      weights <- abs(coef_temp)
      
      # Calculate ratio for variable selection
      ratio_variable <- a * exp(-b * (sample_idx + 1))
      k <- ceiling(ncol(X) * ratio_variable)
      
      # Perform variable selection
      if (k < ncol(X)) {
        # Normalize weights to probabilities
        sample_probs <- weights / sum(weights)
        
        # Sample without sorting to maintain randomness
        subset_vars <- sample(
          seq_len(ncol(X)),
          size = k,
          prob = sample_probs,
          replace = FALSE
        )
      }
    }
    
    # Process results
    min_error <- min(rmsecv)
    opt_iter <- which.max(rmsecv == min_error)
    
    selected_idx <- which(coef_matrix[, opt_iter] != 0)
    selected_wavelengths <- colnames(X)[selected_idx]
    final_weights <- abs(coef_matrix[selected_idx, opt_iter])
    
    # Store results using tidyverse approach
    seed_results <- tibble(wavelengths = selected_wavelengths, seed = current_seed)
    
    all_seeds <- bind_rows(all_seeds, seed_results)
  }
  
  wavelength_counts <- all_seeds %>%
    count(wavelengths) %>%
    arrange(desc(n))
  
  # Wavelengths that appear more than once
  selected_wavelengths <- wavelength_counts %>%
    filter(n > 1) %>%
    pull(wavelengths)
  
  # column indices
  spectral_colnames <- colnames(df)[4:ncol(df)]
  selected_indices <- which(spectral_colnames %in% selected_wavelengths)
  
  
  return(
    list(
      wavelengths = selected_wavelengths,
      indices = selected_indices,
      wavelength_counts = wavelength_counts
    )
  )
}

# PCA ------
#'
#' @references
#' @reference Roger, J. M., et al. (2011). Pre-processing Methods-Wavelength Selection
#'           Methods in Spectroscopy: Sequential and Global Methods.
#' @reference Xiaobo, Z., et al. (2010). Variables selection methods in near-infrared
#'           spectroscopy. Analytica Chimica Acta, 667(1-2), 14-32.
#'
#' @param df Input dataframe with response variable as first column
#' @param n_components Number of principal components to use (default: NULL for automatic selection)
#' @param method Selection method: 'cumulative', 'eigenvalue', or 'loading_weight'
#' @param threshold Selection threshold (interpretation depends on method)
#' @param scale_method Method for scaling: "auto", "pareto", or "none"
#' @return List containing selected wavelengths and their indices
pca_select <- function(df,
                       n_components = NULL,
                       method = "loading_weight",
                       threshold = 0.75,
                       scale_method = "none",
                       seed = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Input validation
  if (!is.data.frame(df)) {
    stop('Input must be a dataframe')
  }
  
  if (ncol(df) < 5) {
    stop('Input dataframe must have at least 5 columns (including response variables)')
  }
  
  # Extract spectral data
  X <- try(as.matrix(df[, 4:ncol(df)]), silent = TRUE)
  if (inherits(X, 'try-error')) {
    stop('Failed to extract spectral data from columns 4 onwards')
  }
  
  # Data preprocessing based on scale_method
  X_processed <- switch(
    scale_method,
    "auto" = scale(X, center = TRUE, scale = TRUE),
    # Standardization
    "pareto" = scale(X, center = TRUE, scale = sqrt(apply(X, 2, sd))),
    # Pareto scaling
    "none" = scale(X, center = TRUE, scale = FALSE),
    # Only centering
    stop("Invalid scaling method")
  )
  
  # Perform PCA
  pca_result <- prcomp(X_processed, center = FALSE, scale = FALSE)
  
  # Calculate all eigenvalues
  all_eigenvalues <- pca_result$sdev^2
  max_possible_components <- length(all_eigenvalues)
  
  # Automatic component selection if n_components is NULL
  if (is.null(n_components)) {
    # Apply Kaiser criterion (eigenvalues > 1)
    n_components <- sum(all_eigenvalues > threshold)
    # Ensure at least 2 components
    n_components <- max(2, min(n_components, max_possible_components))
  } else {
    # Ensure n_components doesn't exceed maximum possible
    n_components <- min(n_components, max_possible_components)
  }
  
  # Get loadings and calculate their importance
  loadings <- pca_result$rotation[, 1:n_components, drop = FALSE]
  var_explained <- (all_eigenvalues / sum(all_eigenvalues))[1:n_components]
  
  # Variable selection based on method
  selected_idx <- switch(
    method,
    "eigenvalue" = {
      # Select based on eigenvalue criterion
      loading_weights <- sweep(loadings^2, 2, var_explained, "*")
      total_weights <- rowSums(loading_weights)
      which(total_weights > mean(total_weights))
    },
    "loading_weight" = {
      # Select based on loading weights
      loading_weights <- abs(loadings)
      max_weights <- apply(loading_weights, 1, max)
      which(max_weights > quantile(max_weights, threshold))
    },
    "cumulative" = {
      # Cumulative variance contribution
      loading_importance <- rowSums(sweep(loadings^2, 2, var_explained, "*"))
      cumsum_importance <- cumsum(sort(loading_importance, decreasing = TRUE))
      cutoff <- which(cumsum_importance >= threshold * sum(loading_importance))[1]
      order(loading_importance, decreasing = TRUE)[1:cutoff]
    },
    stop("Invalid selection method")
  )
  
  if (length(selected_idx) == 0) {
    stop("No variables were selected with the current criteria")
  }
  
  # Calculate metrics for selected variables
  selected_loadings <- loadings[selected_idx, , drop = FALSE]
  selected_wavelengths <- colnames(X)[selected_idx]
  
  # Calculate importance scores using weighted loadings
  importance_scores <- rowSums(sweep(selected_loadings^2, 2, var_explained, "*"))
  
  # Sort by importance
  order_idx <- order(importance_scores, decreasing = TRUE)
  selected_idx <- selected_idx[order_idx]
  selected_wavelengths <- selected_wavelengths[order_idx]
  importance_scores <- importance_scores[order_idx]
  
  # Calculate additional metrics
  results <- list(
    indices = selected_idx,
    wavelengths = selected_wavelengths,
    loadings = selected_loadings,
    importance = importance_scores,
    explained_variance = var_explained,
    metrics = list(
      n_components = n_components,
      n_selected = length(selected_idx),
      total_variance_explained = sum(var_explained),
      eigenvalues = all_eigenvalues[1:n_components],
      condition_number = all_eigenvalues[1] / all_eigenvalues[n_components],
      selection_ratio = length(selected_idx) / ncol(X)
    ),
    parameters = list(
      method = method,
      threshold = threshold,
      scale_method = scale_method
    )
  )
  
  class(results) <- c('pca_select', 'list')
  return(results)
}

# Create list of feature selection functions ----
feature_selection_methods <- list(
  plsr_bcoef = plsr_bcoef_filter,
  spa = run_optimized_spa,
  ga = ga_select,
  uve = uve_select,
  cars = cars_select,
  pca = pca_select,
  none = function(df)
    list(indices = 4:ncol(df), wavelengths = colnames(df)[4:ncol(df)])
)
