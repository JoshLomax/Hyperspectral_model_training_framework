#!/usr/bin/env Rscript
# ===================================================================
# 04-main-pipeline.R
# Main HPC pipeline orchestration script
# Version: 0.1.2
# ===================================================================

# LOAD ----
# suppressPackageStartupMessages({
#   if (!require("pak", quietly = TRUE)) install.packages("pak")
  
  required_packages <- c(
    "tidyverse",
    "glue",
    "tidymodels",
    "pacman",
    "furrr",
    "BGLR",
    "randomForest",
    "caret",
    "glmnet",
    "coda",
    "xgboost",
    "kernlab",
    "ranger",
    "nnet",
    "pls",
    "GetoptLong",
    "LiblineaR",
    "cli",
    "tictoc",
    "GA",
    "prospectr",
    "keras",
    "mdatools",
    "car"
  )
  # pak::pak(required_packages, ask=FALSE)
  pacman::p_load(char = basename(required_packages), install = FALSE)
# })

# SOURCE MODULES -------
source("Scripts/01-data-preparation.R")
source("Scripts/02-feature-selection.R")
source("Scripts/03-model-functions.R")

# COMMAND LINE INTERFACE ----
VERSION <- "0.1.2"
# PARAMETER DEFAULTS ----
#GetoptLong handles defaults automatically
opt = new.env()
opt$cpus <-  1
opt$models <- "all" # c("bayes","rf","xgboost","svm","nnet","ranger",  "plsr", "linear_svr", "elastic_net")
opt$features <- "all"  
opt$seeds <- 5
# opt$verbose <- FALSE


GetoptLong(
  "response=s",       "Response dataset (CSV file path)",
  "var|v=s",          "Response variable column name",
  "input=s",          "Path to preprocessed hyperspectral data (RDS)",
  "output=s",         "Output file path for results (RDS/TXT/CSV)",
  "cpus|c=i",         "Number of CPU cores to use (default: 1)",
  "models|m=s@",      "Models to run (comma-separated or 'all')",
  "features|f=s@",    "Feature selection methods (comma-separated, 'all', or 'none')", 
  "seeds|s=i",        "Number of random seeds (default: 5)",
  "verbose",          "Print detailed progress messages",
  envir = opt
)

# LOGGING TIME ----
tic("Pipeline runtime is ")
# LOGGING SETUP ----
log_message <- function(msg, verbose=TRUE, ...) {
  if (verbose) {
    cli::cli_alert_info(sprintf("%s: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), ...)
  }
}

# FEATURE SELECTION PROCESSING LOGIC ----
# note that the options will be returned as vectors from GetoptLong, so you don't need the last 2 elses 
process_feature_methods <- function(feature_input, fs_methods) {
  # Handle single string input (backward compatibility)
  if (length(feature_input) == 1) {
    if (feature_input == "all") {
      return(names(fs_methods)[names(fs_methods) != "none"])
    } else if (feature_input == "none") {
      return("none")
    } else {
      # Split comma-separated string
      methods <- as.character(str_split(feature_input, ",", simplify = TRUE)) |> trimws()
        # trimws(strsplit(feature_input, ",")[[1]])
      return(methods)
    }
  } else {
    # Handle multiple arguments passed directly
    return(feature_input)
  }
}


#`` Apply feature processing
# opt$features <- c("none")
features_to_run <- process_feature_methods(opt$features, feature_selection_methods)

#`` Validate feature methods
available_features <- names(feature_selection_methods)
invalid_features <- setdiff(features_to_run, available_features)

if (length(invalid_features) > 0) {
  cli::cli_abort("Invalid feature selection methods: {paste(invalid_features, collapse = ', ')}")
}




# MAIN PIPELINE FUNCTION ----
# run_hpc_pipeline <- function() {
  

  
  #`` Step 1: Load and validate data
  log_message("Loading data from {opt$input}", verbose = opt$verbose)
  
  hs_data <- tryCatch({
    readRDS(opt$input)
  }, error = function(e) {
    cli::cli_abort("Failed to load hyperspectral data: {e$message}")
  })
  
  response_data <- tryCatch({
    read_csv(opt$response, show_col_types = FALSE)
  }, error = function(e) {
    cli::cli_abort("Failed to load response data: {e$message}")
  })
  
  #`` Validate response variable exists
  if (!opt$var %in% names(response_data)) {
    cli::cli_abort("Response variable '{opt$var}' not found in response data")
  }
  
  #`` Step 2: Merge response with hyperspectral data
  log_message("Merging response variable '{opt$var}' with hyperspectral data", verbose = opt$verbose)
  
  combined_data <- response_data |>
    select(all_of(opt$var), ID) |>
    inner_join(hs_data, by = "ID")
  
  if (nrow(combined_data) == 0) {
    cli::cli_abort("No matching samples found between response and hyperspectral data")
  }
  
  log_message("Combined dataset: {nrow(combined_data)} samples × {ncol(combined_data)} columns",
              verbose = opt$verbose)
  
  # #`` Step 3: Apply feature selection (if specified)
  # if (opt$features != "none") {
  #   log_message("Applying {opt$features} feature selection")
  #   
  #   if (!opt$features %in% names(feature_selection_methods)) {
  #     cli::cli_abort("Unknown feature selection method: {opt$features}")
  #   }
  #   
  #   selected_features <- feature_selection_methods[[opt$features]](combined_data) #e.g., opt$features = "uve"
  #   
  #   #`` Apply feature selection
  #   feature_cols <- selected_features$indices + 3  # Adjust for metadata columns
  #   combined_data <- combined_data[, c(1:3, feature_cols)]
  #   
  #   log_message("Selected {length(selected_features$indices)} features")
  #   # toc()
  # }
  
  #`` Step 4: Prepare datasets (flesh, skin, combined)
  datasets <- prepare_datasets(combined_data)
  
  #`` Step 5: Determine models to run
  if (opt$models == "all") {
    models_to_run <- names(model_registry)
  } else {
    # opt$models=names(model_registry[4:5])
    # opt$models=c("svm, nnet")
    models_to_run <- as.character(str_split(opt$models, ",", simplify = TRUE)) |> trimws()
    invalid_models <- setdiff(models_to_run, names(model_registry))
    
    if (length(invalid_models) > 0) {
      cli::cli_abort("Invalid models: {paste(invalid_models, collapse = ', ')}")
    }
  }
  
  log_message("Models to run: {paste(models_to_run, collapse = ', ')}", verbose = opt$verbose)
  
  #`` Step 6: Generate random seeds
  set.seed(1262753)  # Master seed for reproducibility
  seeds <- sample(1:1e6, size = opt$seeds) #e.g., opt$seeds = 5
  
  log_message("Using {length(seeds)} random seeds", verbose = opt$verbose)
  
  #`` Step 7: Setup parallel processing
  if (opt$cpus > 1) {
    log_message("Setting up parallel processing with {opt$cpus} cores", verbose = opt$verbose)
    plan(multisession, workers = opt$cpus)
  } else {
    plan(sequential)
  }
  
  # #`` Step 8: Create parameter grid with multiple datasets and feature selection methods
  # data_prep_grid <- expand_grid(dataset_type = names(datasets),
  #                               feature_method = features_to_run)
  # data_list <- data_prep_grid |>
  #   pmap(.f = function(dataset_type, feature_method) {
  #     cli::cli_h1(
  #       "Preparing data: {.field {dataset_type}} | {.field {feature_method}}"
  #     )
  #     # cli_progress_bar("Cleaning data", total = 100)
  #     source("Scripts/01-data-preparation.R")
  #     source("Scripts/02-feature-selection.R")
  #     
  #     #`` Get dataset and apply feature selection
  #     current_data <- datasets[[dataset_type]]
  #     # tic(quiet = T)
  #     set.seed(seed)
  #     # Apply feature selection
  #     
  #     if (feature_method != "none") {
  #       cli::cli_progress_step("Running {.field {feature_method}} feature selection")
  #       set.seed(seed)
  #       selected_features <- feature_selection_methods[[feature_method]](current_data, seed=seed)
  #       feature_cols <- selected_features$indices + 3
  #       current_data <- current_data[, c(1:3, feature_cols)]
  #       n_features_selected <- length(selected_features$indices)
  #     } else {
  #       n_features_selected <- ncol(current_data) - 3
  #     }
  #     
  #     
  #   })
  # 
  
  
  #`` Step 8: Create parameter grid with multiple feature methods
  #' 810 combinations
  param_grid <- expand_grid(
    # model = models_to_run,
    seed = seeds,
    dataset_type = names(datasets),
    feature_method = features_to_run
  )
  
  # log_message("Total model runs: {nrow(param_grid)} ({length(models_to_run)} model{?s} x {length(features_to_run)} feature selection method{?s}) x {length(seeds)} seed{?s} x {length(datasets)} dataset type{?s}",
  #             verbose = opt$verbose)

  log_message("Total dataset runs: {nrow(param_grid)} ({length(features_to_run)} feature selection method{?s}) x {length(seeds)} seed{?s} x {length(datasets)} dataset type{?s}",
              verbose = opt$verbose)
  
  #`` Step 9: Run models with feature selection
  results <- param_grid |>
    
    # mutate(row_id = row_number()) |>
    
    pmap(
      .f = function(seed, dataset_type, feature_method) {
        cli::cli_h1(
          "Processing data: {.field {dataset_type}} | {.field {feature_method}} | {.field seed
{seed}}"
        )
        # cli_progress_bar("Cleaning data", total = 100)
        source("Scripts/01-data-preparation.R")
        source("Scripts/02-feature-selection.R")
        source("Scripts/03-model-functions.R")
        # tryCatch({
          #`` Get dataset and apply feature selection
          current_data <- datasets[[dataset_type]]
          
          set.seed(seed)
          # Apply feature selection
          tic(quiet = T)
          if (feature_method != "none") {
            cli::cli_progress_step("Running {.field {feature_method}} feature selection")
            # set.seed(seed)
            selected_features <- feature_selection_methods[[feature_method]](current_data, seed = seed)
            feature_cols <- selected_features$indices + 3
            current_data <- current_data[, c(1:3, feature_cols)]
            n_features_selected <- length(selected_features$indices)
          } else {
            n_features_selected <- ncol(current_data) - 3
          }
          fs_run_time <- toc(quiet = T)
          
          # Run model
          model_results <- models_to_run |> 
            map(.f = function(model){
              cli::cli_progress_step("Running {.field {model}} model with {.val {n_features_selected}} features")
              tic(quiet = T)
              # cli::cli_alert_info("Running {.field {model}} model with {.val {n_features_selected}} features")
              model_func <- model_registry[[model]]
              result <- model_func(current_data, seed = seed)
              model_run_time <- toc(quiet = T)
              cli::cli_progress_step("Processing model results")
              # Process results
              result |>
                # filter(data_type == dataset_type) |>
                # select(model, data_type, mae, mse, r2, rmse, nrmse, accuracy, n_samples, seed) %>%
                # mutate(across(where(is.numeric), ~ round(.x, digits = 2))) %>%
                # arrange(desc(accuracy))
                mutate(
                  # model = model,
                  # seed = seed,
                  # feature_method = feature_method,
                  # data_type = dataset_type,
                  # n_features = n_features_selected,
                  model_run_time = model_run_time$toc - model_run_time$tic,
                  # r2 = as.numeric(r2),
                  # corr = as.numeric(corr),
                  across(where(is.numeric), ~ round(.x, digits = 2))
                )
              
            })
              model_results |>
                list_rbind() |> 
                mutate(
                  feature_method = feature_method,
                  data_type = dataset_type,
                  n_features = n_features_selected,
                  fs_run_time = fs_run_time$toc - fs_run_time$tic,
                )
          
        # }, error = function(e) {
        #   cli::cli_alert_danger("Error in {.field {dataset_type}}/{.field {feature_method}}")
        #   cli::cli_alert_warning("Message: {.val {e$message}}")
        #   
        #   tibble(
        #     model = NA,
        #     data_type = dataset_type,
        #     seed = seed,
        #     feature_method = feature_method,
        #     n_features = NA,
        #     error = as.character(e$message)
        #   )
        })
    #   }# , .options = furrr_options(seed = TRUE, packages = required_packages)
    # ) 
  
  
  #`` Step 10: Process and save results
  log_message("Combining all model results", verbose = opt$verbose)
  
  preprocessing_method <- basename(opt$input) |> 
    tools::file_path_sans_ext()
  
  final_results <- results |>
    list_rbind() |> 
    mutate(
      preprocessing = preprocessing_method,
      response_var = opt$var,
      timestamp = Sys.time()
    ) 
  
  # #`` Summary statistics
  # if (nrow(final_results) > 0 && "r2" %in% names(final_results)) {
  #   summary_stats <- final_results |>
  #     group_by(model, data_type) |>
  #     summarise(
  #       mean_r2 = mean(r2, na.rm = TRUE),
  #       sd_r2 = sd(r2, na.rm = TRUE),
  #       .groups = "drop"
  #     )
  #   
  #   cli::cli_h2("Top performing models:")
  #   print(head(summary_stats |> arrange(desc(mean_r2)), 10))
  # }
  
  #`` Save results
  log_message("Saving results to {opt$output}", verbose = opt$verbose)
  # opt$output <- "tests/test_results.out"
  # final_results <- response_data
  out_filetype <- tolower(fs::path_ext(opt$output))
  # out_filetype <- fs::path_ext(outputfn)
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  switch (out_filetype,
    rds = write_rds(final_results, file = opt$output),
    csv = write_csv(final_results, file = opt$output),
    tsv = write_tsv(final_results, file = opt$output),
    txt = write_tsv(final_results, file = opt$output),
    cli::cli_abort("Unable to parse filetype extension for {.field {opt$output}}, see documentation for valid file types")
  )
  
  # if (!grepl("rds|csv|tsv|txt", out_filetype , ignore.case = TRUE)) {
  #   cli::cli_abort("Unable to parse filetype extension for '{opt$output}', see documentation for valid file types")
  # }
  # final_results_table <- case_when(grepl("rds", out_filetype, ignore.case = TRUE) ~ write_rds(final_results, 
  #                                                                                             file = opt$output),
  #           grepl("csv", out_filetype, ignore.case = TRUE) ~ write_csv(final_results, file = opt$output),
  #           grepl("tsv", out_filetype, ignore.case = TRUE) ~ write_tsv(final_results, file = opt$output),
  #           grepl("txt", out_filetype, ignore.case = TRUE) ~ write_tsv(final_results, file = opt$output)) 
  # 
  toc()
  
  # #`` Save session info for reproducibility
  # session_info_file <- sub("\\.rds$", "_session_info.txt", opt$output)
  # writeLines(capture.output(sessionInfo()), session_info_file)
  
  # toc()
  
  # cli::cli_alert_success("Pipeline completed successfully!")
  # cli::cli_alert_info("Results saved to: {opt$output}")
  # 
  # return(invisible(final_results))
# }
# 
# #`` ERROR HANDLING WRAPPER
# tryCatch({
#   run_hpc_pipeline()
# }, error = function(e) {
#   cli::cli_alert_danger("Pipeline failed: {e$message}")
#   quit(status = 1)
# })