# ===================================================================
# 01-data-preparation.R
# Data preprocessing and train-test splitting functions
# ===================================================================

#`` STREAMLINED DATASET PREPARATION
prepare_datasets <- function(df) {
  stopifnot(is.data.frame(df), "Type" %in% names(df), nrow(df) > 0)
  
  type_summary <- df |>
    count(Type, name = "n_samples") |>
    filter(Type %in% c('F', 'S'))
  
  cli::cli_inform(c(
    "i" = "Dataset composition:",
    " " = glue::glue_data(type_summary, "{Type}: {n_samples} samples")
  ))
  
  list(
    combined = df,
    flesh = df |> filter(Type == 'F'),
    skin = df |> filter(Type == 'S')
  )
}

#`` MODEL DATA PREPARATION
prepare_model_data <- function(df,
                               response_col = 1,
                               test_prop = 0.2,
                               seed = 123) {
  predictor_cols <- 4:ncol(df)
  
  scaled_data <- df |>
    mutate(across(all_of(c(response_col,predictor_cols)), ~ as.numeric(scale(.x))))
  
  X <- as.matrix(scaled_data[, predictor_cols])
  y <- pull(scaled_data, response_col)
  
  set.seed(seed)
  train_indices <- caret::createDataPartition(y, p = 1 - test_prop, list = FALSE)[, 1]
  
  list(
    X = X,
    y = y,
    X_train = X[train_indices, ],
    X_test = X[-train_indices, ],
    y_train = y[train_indices],
    y_test = y[-train_indices],
    train_idx = train_indices,
    test_idx = setdiff(seq_len(nrow(df)), train_indices)
  )
}

#`` METRICS CALCULATION
calculate_metrics <- function(test_y_true,
                              test_y_pred,
                              n_predictors = NULL,
                              train_y_true,
                              train_y_pred) {
  # Input validation
  stopifnot(
    length(test_y_true) == length(test_y_pred),
    length(train_y_true) == length(train_y_pred)
  )
  
  # Convert predictions to numeric
  test_y_pred <- as.numeric(test_y_pred)
  train_y_pred <- as.numeric(train_y_pred)
  
  # Calculate residuals
  test_residuals <- test_y_true - test_y_pred
  train_residuals <- train_y_true - train_y_pred
  
  # Sample sizes
  n_test <- length(test_y_true)
  n_train <- length(train_y_true)
  
  # Helper function for R-squared
  calc_r2 <- function(actual, predicted) {
    1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
  }
  
  # Helper function for all metrics
  calc_metrics <- function(actual,
                           predicted,
                           residuals,
                           dataset_name) {
    r2 <- calc_r2(actual, predicted)
    n <- length(actual)
    
    # Adjusted R-squared
    adj_r2 <- if (!is.null(n_predictors) &&
                  n_predictors < (n - 1)) {
      1 - ((1 - r2) * (n - 1) / (n - n_predictors - 1))
    } else {
      if (!is.null(n_predictors)) {
        warning(
          paste(
            "Number of predictors too large relative to",
            dataset_name,
            "sample size"
          )
        )
      }
      NA
    }
    
    rmse <- sqrt(mean(residuals^2))
    sd_actual <- sd(actual)
    rpd <- sd_actual/rmse
    
    list(
      n = n,
      mae = mean(abs(residuals)),
      mse = mean(residuals^2),
      rmse = rmse,
      nrmse = sqrt(mean(residuals^2)) / diff(range(actual)),
      r2 = r2,
      adj_r2 = as.numeric(adj_r2),
      corr = cor(predicted, actual),
      rpd = rpd,
      residual_std = sd(residuals),
      mean_residual = mean(residuals),
      median_residual = median(residuals)
    )
  }
  
  # Calculate metrics for both datasets
  test_metrics <- calc_metrics(test_y_true, test_y_pred, test_residuals, "test")
  train_metrics <- calc_metrics(train_y_true, train_y_pred, train_residuals, "train")
  
  # Return comprehensive metrics with clear naming
  tibble(
    # Sample sizes
    n_test = test_metrics$n,
    n_train = train_metrics$n,
    
    # Test metrics
    test_mae = test_metrics$mae,
    test_mse = test_metrics$mse,
    test_rmse = test_metrics$rmse,
    test_nrmse = test_metrics$nrmse,
    test_r2 = test_metrics$r2,
    test_adj_r2 = test_metrics$adj_r2,
    test_corr = test_metrics$corr,
    test_rpd = test_metrics$rpd,
    test_residual_std = test_metrics$residual_std,
    test_mean_residual = test_metrics$mean_residual,
    test_median_residual = test_metrics$median_residual,
    
    # Train metrics
    train_mae = train_metrics$mae,
    train_mse = train_metrics$mse,
    train_rmse = train_metrics$rmse,
    train_nrmse = train_metrics$nrmse,
    train_r2 = train_metrics$r2,
    train_adj_r2 = train_metrics$adj_r2,
    train_corr = train_metrics$corr,
    train_rpd = train_metrics$rpd,
    train_residual_std = train_metrics$residual_std,
    train_mean_residual = train_metrics$mean_residual,
    train_median_residual = train_metrics$median_residual
  )
}
