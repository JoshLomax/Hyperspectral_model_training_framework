# Functions ----
required_packages <- c(
  "tidyverse",  # For dplyr, purrr, tibble operations
  "glue",       # For string interpolation
  "hyperSpec",  # For spectral data handling and baseline corrections
  "signal",     # For sgolayfilt function used in Derivative
  "prospectr"   # For snv (Standard Normal Variate) transformation
)
pak::pak(required_packages)
pacman::p_load(char = required_packages, install = F)
#' Process spectral data with multiple transformations
#' @param df Input dataframe containing spectral data
#' @param transformations Vector of transformation names to apply
#' @param spectra spectral data without meta data columns
#' @param smooth_window smoothing function parameter
#' @return Processed dataframe and saves results

process_spectra <- function(df,
                            spectra = df[, -c(1:2)], # account for meta data columns
                            transformations = c("smooth",
                                                "snv",
                                                "baseline1",
                                                "baseline2",
                                                "derivative1",
                                                "derivative2"),
                            smooth_window = 5) {
  # Input validation
  validate_inputs <- function(spectra) {
    if(any(is.na(spectra))) stop("Spectra contain NA values")
    if(any(spectra < 0)) stop("Negative values in spectra")
    if(!is.numeric(as.matrix(spectra))) stop("Spectra must be numeric")
    invisible(TRUE)
  }
  
  # Define transformation functions
  transformation_functions <- list(
    smooth = function(data) {
      suppressWarnings({
        smoothed_data <- SmoothFast(data, smooth_window)
        colnames(smoothed_data) <- colnames(data)
        return(as_tibble(smoothed_data))
      })
    },
    
    baseline1 = function(data) {
      if(any(is.na(data))) stop("Cannot perform baseline correction on NA values")
      Baseline(data, Order = 1)
    },
    
    baseline2 = function(data) {
      if(any(is.na(data))) stop("Cannot perform baseline correction on NA values")
      Baseline(data, Order = 2)
    },
    
    derivative1 = function(data) {
      derived_data <- t(Derivative(data, Order = 1))
      colnames(derived_data) <- colnames(data)
      as_tibble(derived_data)
    },
    
    derivative2 = function(data) {
      derived_data <- t(Derivative(data, Order = 2))
      colnames(derived_data) <- colnames(data)
      as_tibble(derived_data)
    },
    
    snv = function(data) {
      if(any(is.na(data))) stop("Cannot perform SNV on NA values")
      as_tibble(standardNormalVariate(data))
    }
  )
  # Helper function to generate all possible combinations
  get_transformation_combinations <- function(transformations, max_depth = 4) {
    # Helper function to check valid combination
    is_valid_combination <- function(combo) {
      # Count derivatives and baselines
      derivative_count <- sum(str_detect(combo, "derivative"))
      baseline_count <- sum(str_detect(combo, "baseline"))
      
      # Valid if at most one of each type
      return(derivative_count <= 1 && baseline_count <= 1)
    }
    
    # Generate all combinations and filter
    map(1:max_depth, function(n) {
      combn(transformations, n, simplify = FALSE) %>%
        keep(is_valid_combination)
    }) %>%
      unlist(recursive = FALSE)
  }
  
  # Main processing pipeline
  process_and_save <- function(df, transformation_sequence) {
    
    # Create transformation name
    transform_name <- paste(transformation_sequence, collapse = "_")
    
    message(glue::glue("\nProcessing transformation: {transform_name}"))
    
    tryCatch({
      # Apply transformations
      result <- spectra
      for (transform in transformation_sequence){
        message(glue::glue("applying: {transform}..."))
        result <- transformation_functions[[transform]](result)
        
        if(any(is.na(result))) {
          stop(glue::glue("{transform} produced NA values"))
        }
        if(any(!is.finite(results))) {
          stop(glue:glue("{transform} produced infinite values"))
        }
      }
      
     # Combine with metadata
      final_df <- bind_cols(
        df[, 1:2],
        as_tibble(result)
      )
      save_path <- glue::glue("input_data/{transform_name}_data")
      dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
      saveRDS(final_df, glue::glue("{save_path}/{transform_name}.RDS"))

      message(glue::glue("  ✓ Successfully saved to {save_path}"))
    },
    error = function(e) {
      message(glue::glue("  ✗ Failed: {e$message}"))
    })
  }

  # Main execution
  validate_inputs(spectra)
  
  combinations <- get_transformation_combinations(transformations)

  # Process all combinations with progress reporting
  message(glue::glue("\nProcessing {length(combinations)} transformation combinations..."))
  walk(combinations, ~process_and_save(df, .))

  message("\nProcessing complete. Check the 'input_data' directory for results.")
}

# smoothing function:
SmoothFast <- function(spectra,
                       smooth_window){                         
  #Create smoothing matrix: 
  Mat <- matrix(0,
                length((smooth_window+1):
                         (ncol(spectra)-smooth_window)),
                2*smooth_window+1)
  
  for(j in 1:nrow(Mat)){
    Mat[j, ] <- seq(j, j + 2 * smooth_window, 1)
  }
  
  #Smoothing spectra using matrix operations:
  newspectra<-matrix(0,
                     nrow(spectra),
                     length((smooth_window+1):
                              (ncol(spectra)-
                                 smooth_window)))
  
  for(i in 1:nrow(Mat)) {
    newspectra[, i] <- apply(spectra[, Mat[i, ]], 
                             1, 
                             mean)
  }
  
  #Add front and end tails (not smoothed):
  fronttail <- newspectra[, 1]
  endtail <- newspectra[, ncol(newspectra)]
  
  for(k in 1:(smooth_window - 1)){
    fronttail <- data.frame(fronttail, newspectra[, 1])
    endtail <- data.frame(endtail, 
                          newspectra[, ncol(newspectra)])
  }
  
  newspectra<-data.frame(fronttail,newspectra,endtail)
  
  return(newspectra)
}

# baseline function:
Baseline <- function(df, Order = 1){
  
  # Convert mydata to an hyperSpec S4 object:
  mydataHS <- new("hyperSpec",
                  spc = as.matrix(df[, -c(1:2)]),
                  wavelength = as.numeric(colnames(df[, -c(1:2)])))
  
  baseline <- spc.fit.poly.below(fit.to = mydataHS, 
                                 poly.order = Order)
  
  # Extract baseline spectra preserving column names
  wl <- mydataHS@wavelength
  mybaseline <- data.frame(setNames(as.data.frame(baseline[[]], baseline@data$spc), wl))
  
  # Baseline removal with preserved column names
  oldspectra <- as.data.frame(mydataHS[[]])
  baselinespectra <- as.data.frame(baseline[[]])
  newspectra <- oldspectra - baselinespectra
  
  return(newspectra)
}

# derivative function:
Derivative <- function(df, Order = 1){
  newspectra <- apply(
    df,
    1,
    FUN = sgolayfilt,
    p = 2,
    n = 5,
    m = Order,
    ts = 1
  )
  return(newspectra)
}

df2 <- readRDS('data/absorbance_data.RDS')

results <- process_spectra(
  df = df2,
  spectra = df2[,-c(1:2)]
  )
