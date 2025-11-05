#' Helper Functions for Fire Recovery Analysis
#' 
#' Collection of utility functions for data preprocessing and manipulation
#' 
#' @author Your Name
#' @date 2025

library(terra)
library(sf)
library(dplyr)
library(stringr)

#' Convert Fire Year Format to Date Indices
#'
#' Converts fire year format (e.g., 200102) to start and end indices in time series
#'
#' @param fire_year Integer fire year in YYYYMM format (e.g., 200102 = June 2001 to July 2002)
#' @param n_per_year Number of observations per year (default = 46 for 8-day MODIS)
#' @return List containing start_idx, end_idx, start_year, end_year
#' @export
#' @examples
#' fire_indices <- convert_fire_year_to_indices(200607, n_per_year = 46)
convert_fire_year_to_indices <- function(fire_year, n_per_year = 46) {
  # Parse fire year format (e.g., 200102 = June 2001 to July 2002)
  fire_year_str <- sprintf("%06d", fire_year)
  start_year <- as.numeric(substr(fire_year_str, 1, 4))
  end_year_suffix <- substr(fire_year_str, 5, 6)
  
  # Handle special case where end year is "00" (represents 2000)
  if (end_year_suffix == "00") {
    end_year <- 2000
  } else {
    century <- floor(start_year / 100) * 100
    end_year <- century + as.numeric(end_year_suffix)
  }
  
  # Fire season: June 1st of start_year to July 31st of end_year
  # Assuming time series starts from 2000-01-01 with 8-day intervals
  
  # Calculate start index (June 1st of start year)
  years_from_2000_start <- start_year - 2000
  june_start_step <- round(((6-1) * 30 + 1) / 8)  # June 1st ≈ day 152, step ≈ 19
  fire_start_idx <- years_from_2000_start * n_per_year + june_start_step
  
  # Calculate end index (July 31st of end year)  
  years_from_2000_end <- end_year - 2000
  july_end_step <- round(((7-1) * 30 + 31) / 8)  # July 31st ≈ day 212, step ≈ 27
  fire_end_idx <- years_from_2000_end * n_per_year + july_end_step
  
  return(list(
    start_idx = max(1, fire_start_idx),
    end_idx = fire_end_idx,
    start_year = start_year,
    end_year = end_year
  ))
}

#' Extract Dates from Column Names
#'
#' Extracts and converts dates from MODIS raster column names
#'
#' @param column_names Character vector of column names
#' @return Vector of Date objects
#' @export
extract_dates_from_names <- function(column_names) {
  cat("Extracting dates from", length(column_names), "column names\n")
  
  # Try multiple patterns to match different formats
  date_parts <- NULL
  
  # Pattern 1: BU_LAI_FPAR_2000_2000057_b1 format
  if (is.null(date_parts) || all(is.na(date_parts))) {
    date_parts <- stringr::str_extract(column_names, "\\d{4}_\\d{7}")
    if (!all(is.na(date_parts))) cat("Matched pattern YYYY_YYYYJJJ\n")
  }
  
  # Pattern 2: Direct YYYYJJJ format (e.g., 2000057)
  if (is.null(date_parts) || all(is.na(date_parts))) {
    date_parts <- stringr::str_extract(column_names, "\\d{7}")
    if (!all(is.na(date_parts))) cat("Matched pattern YYYYJJJ\n")
  }
  
  # Pattern 3: X followed by date (e.g., X2000.01.01)
  if (is.null(date_parts) || all(is.na(date_parts))) {
    date_parts <- stringr::str_extract(column_names, "X(\\d{4}[._]\\d{2,3}[._]?\\d{0,2})")
    if (!all(is.na(date_parts))) cat("Matched pattern X date format\n")
  }
  
  # Convert to dates based on the pattern matched
  dates <- sapply(seq_along(date_parts), function(i) {
    x <- date_parts[i]
    if (is.na(x)) return(NA)
    
    tryCatch({
      if (grepl("_", x) && nchar(x) > 8) {
        # Format: 2000_2000057
        parts <- strsplit(x, "_")[[1]]
        if (length(parts) == 2) {
          year_part <- parts[2]
          if (nchar(year_part) == 7) {
            year <- as.numeric(substr(year_part, 1, 4))
            julian_day <- as.numeric(substr(year_part, 5, 7))
            date <- as.Date(julian_day - 1, origin = paste0(year, "-01-01"))
            return(as.character(date))
          }
        }
      } else if (nchar(x) == 7 && !grepl("[A-Za-z]", x)) {
        # Format: 2000057
        year <- as.numeric(substr(x, 1, 4))
        julian_day <- as.numeric(substr(x, 5, 7))
        date <- as.Date(julian_day - 1, origin = paste0(year, "-01-01"))
        return(as.character(date))
      }
      return(NA)
    }, error = function(e) {
      return(NA)
    })
  })
  
  result_dates <- as.Date(dates)
  valid_dates <- sum(!is.na(result_dates))
  
  cat("Successfully converted", valid_dates, "out of", length(dates), "dates\n")
  
  # If no dates could be extracted, create a simple sequence
  if (valid_dates == 0) {
    cat("No dates extracted, creating simple sequence from 2000-01-01\n")
    start_date <- as.Date("2000-01-01")
    result_dates <- seq(start_date, by = 8, length.out = length(column_names))
  }
  
  return(result_dates)
}

#' Validate Input Data
#'
#' Checks that required data objects exist and have correct structure
#'
#' @param cluster_data Data frame with cluster assignments
#' @param fire_polygons SF object with fire polygons
#' @param fpar_stack SpatRaster with FPAR time series
#' @param unburnt_reference Data frame with unburnt mean time series
#' @return List with validation results
#' @export
validate_input_data <- function(cluster_data, fire_polygons, fpar_stack, unburnt_reference) {
  
  validation_results <- list(
    valid = TRUE,
    messages = c()
  )
  
  # Check cluster data
  if (!("class" %in% names(cluster_data)) || !("x" %in% names(cluster_data)) || !("y" %in% names(cluster_data))) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages, 
                                     "Cluster data must contain 'class', 'x', 'y' columns")
  }
  
  # Check fire polygons
  if (!inherits(fire_polygons, "sf")) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages, 
                                     "Fire polygons must be an sf object")
  }
  
  if (!("FireYear" %in% names(fire_polygons))) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages, 
                                     "Fire polygons must contain 'FireYear' column")
  }
  
  # Check FPAR stack
  if (!inherits(fpar_stack, "SpatRaster")) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages, 
                                     "FPAR stack must be a SpatRaster object")
  }
  
  # Check unburnt reference
  required_cols <- c("cluster", "fireyear")
  if (!all(required_cols %in% names(unburnt_reference))) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages, 
                                     "Unburnt reference must contain 'cluster' and 'fireyear' columns")
  }
  
  # Check for time series columns
  ts_cols <- grep("^BU_LAI_FPAR_|^X[0-9]", names(unburnt_reference), value = TRUE)
  if (length(ts_cols) == 0) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(validation_results$messages, 
                                     "Unburnt reference must contain time series columns")
  }
  
  # Print validation summary
  if (validation_results$valid) {
    cat("✓ All input data validation checks passed\n")
    cat("  - Clusters:", length(unique(cluster_data$class)), "\n")
    cat("  - Fire years:", length(unique(fire_polygons$FireYear)), "\n")
    cat("  - FPAR layers:", nlyr(fpar_stack), "\n")
    cat("  - Time series columns:", length(ts_cols), "\n")
  } else {
    cat("✗ Input data validation FAILED:\n")
    for (msg in validation_results$messages) {
      cat("  -", msg, "\n")
    }
  }
  
  return(validation_results)
}

#' Create Output Directories
#'
#' Creates necessary output directory structure
#'
#' @param base_path Base output path
#' @return List of created directory paths
#' @export
create_output_directories <- function(base_path = "outputs") {
  
  dirs <- list(
    base = base_path,
    rasters = file.path(base_path, "rasters"),
    impact = file.path(base_path, "rasters", "impact"),
    recovery = file.path(base_path, "rasters", "recovery"),
    plots = file.path(base_path, "plots"),
    timeseries = file.path(base_path, "plots", "timeseries"),
    summary = file.path(base_path, "plots", "summary"),
    diagnostics = file.path(base_path, "diagnostics"),
    tables = file.path(base_path, "tables")
  )
  
  # Create all directories
  for (dir_name in names(dirs)) {
    dir.create(dirs[[dir_name]], recursive = TRUE, showWarnings = FALSE)
  }
  
  cat("Created output directory structure:\n")
  for (dir_name in names(dirs)) {
    cat("  -", dir_name, ":", dirs[[dir_name]], "\n")
  }
  
  return(dirs)
}

#' Log Message to File
#'
#' Appends timestamped message to log file
#'
#' @param message Message to log
#' @param log_file Path to log file
#' @param level Log level (INFO, WARNING, ERROR)
#' @export
log_message <- function(message, log_file, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s\n", timestamp, level, message)
  cat(log_entry)
  write(log_entry, file = log_file, append = TRUE)
}