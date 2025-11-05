#' Reference Value Calculation Functions
#' 
#' Functions for identifying modal peak periods and calculating reference values
#' 
#' @author Your Name
#' @date 2025

#' Find Reference Values Using Modal Peak Method
#'
#' Identifies modal peak position across years and extracts reference values
#'
#' @param unburnt_ts Numeric vector of unburnt FPAR time series
#' @param n_per_year Number of observations per year (default = 46)
#' @param window_size Size of peak detection window in observations (default = 5)
#' @return List containing reference values and peak information
#' @export
#' @examples
#' ref_info <- find_reference_values(unburnt_timeseries, n_per_year = 46, window_size = 5)
find_reference_values <- function(unburnt_ts, n_per_year = 46, window_size = 5) {
  
  # Handle short time series
  if (length(unburnt_ts) < n_per_year) {
    return(list(
      reference_values = rep(mean(unburnt_ts, na.rm = TRUE), window_size),
      reference_indices = 1:min(window_size, length(unburnt_ts)),
      modal_peak_position = 1,
      yearly_peak_indices = c(),
      peak_position_in_year = 1,
      n_years = 0,
      yearly_peaks = list()
    ))
  }
  
  # Calculate number of complete years
  n_years <- floor(length(unburnt_ts) / n_per_year)
  
  # Find peak windows for each year
  yearly_peaks <- list()
  yearly_peak_positions <- c()
  
  for (year in 1:n_years) {
    year_start <- (year - 1) * n_per_year + 1
    year_end <- year * n_per_year
    
    if (year_end > length(unburnt_ts)) {
      year_end <- length(unburnt_ts)
    }
    
    year_data <- unburnt_ts[year_start:year_end]
    
    # Find the best window with highest mean
    best_mean <- -Inf
    best_start <- 1
    
    # Try all possible windows of size window_size within the year
    max_start <- length(year_data) - window_size + 1
    
    for (start_pos in 1:max_start) {
      window_data <- year_data[start_pos:(start_pos + window_size - 1)]
      
      # Skip if any NA values
      if (any(is.na(window_data))) next
      
      window_mean <- mean(window_data, na.rm = TRUE)
      
      if (window_mean > best_mean) {
        best_mean <- window_mean
        best_start <- start_pos
      }
    }
    
    yearly_peaks[[year]] <- list(
      position = best_start,
      absolute_start = year_start + best_start - 1,
      values = year_data[best_start:(best_start + window_size - 1)],
      mean_value = best_mean
    )
    
    yearly_peak_positions <- c(yearly_peak_positions, best_start)
  }
  
  # Find the modal (most common) peak position across years
  if (length(yearly_peak_positions) > 0) {
    position_table <- table(yearly_peak_positions)
    modal_position <- as.numeric(names(position_table)[which.max(position_table)])
  } else {
    modal_position <- 1
  }
  
  # Collect all reference values from the modal position across all years
  all_reference_values <- c()
  all_reference_indices <- c()
  
  for (year in 1:n_years) {
    year_start <- (year - 1) * n_per_year + 1
    
    # Calculate absolute indices for the modal position in this year
    modal_start_abs <- year_start + modal_position - 1
    modal_end_abs <- modal_start_abs + window_size - 1
    
    # Check if indices are valid
    if (modal_end_abs <= length(unburnt_ts)) {
      year_reference_values <- unburnt_ts[modal_start_abs:modal_end_abs]
      
      # Only include if no NA values
      if (!any(is.na(year_reference_values))) {
        all_reference_values <- c(all_reference_values, year_reference_values)
        all_reference_indices <- c(all_reference_indices, modal_start_abs:modal_end_abs)
      }
    }
  }
  
  # If we have reference values, use them; otherwise fallback
  if (length(all_reference_values) > 0) {
    reference_values <- all_reference_values
    reference_indices <- all_reference_indices
  } else {
    # Fallback: use the first year's modal position
    modal_start_abs <- modal_position
    modal_end_abs <- modal_position + window_size - 1
    if (modal_end_abs <= length(unburnt_ts)) {
      reference_values <- unburnt_ts[modal_start_abs:modal_end_abs]
      reference_indices <- modal_start_abs:modal_end_abs
    } else {
      reference_values <- rep(mean(unburnt_ts, na.rm = TRUE), window_size)
      reference_indices <- 1:window_size
    }
  }
  
  # Create yearly peak indices for all years (for recovery comparison)
  yearly_peak_indices <- c()
  for (year in 1:n_years) {
    year_start <- (year - 1) * n_per_year + 1
    modal_indices_year <- (year_start + modal_position - 1):(year_start + modal_position + window_size - 2)
    
    # Only include valid indices
    valid_indices <- modal_indices_year[modal_indices_year <= length(unburnt_ts)]
    yearly_peak_indices <- c(yearly_peak_indices, valid_indices)
  }
  
  return(list(
    reference_values = reference_values,
    reference_indices = reference_indices,
    modal_peak_position = modal_position,
    yearly_peak_indices = yearly_peak_indices,
    peak_position_in_year = modal_position,
    n_years = n_years,
    yearly_peaks = yearly_peaks
  ))
}

#' Calculate Reference Statistics
#'
#' Calculates summary statistics for reference values
#'
#' @param reference_info Output from find_reference_values()
#' @return Data frame with reference statistics
#' @export
calculate_reference_statistics <- function(reference_info) {
  
  stats <- data.frame(
    modal_peak_position = reference_info$modal_peak_position,
    n_years_used = reference_info$n_years,
    n_reference_values = length(reference_info$reference_values),
    reference_mean = mean(reference_info$reference_values, na.rm = TRUE),
    reference_sd = sd(reference_info$reference_values, na.rm = TRUE),
    reference_min = min(reference_info$reference_values, na.rm = TRUE),
    reference_max = max(reference_info$reference_values, na.rm = TRUE),
    reference_median = median(reference_info$reference_values, na.rm = TRUE),
    cv = sd(reference_info$reference_values, na.rm = TRUE) / 
      mean(reference_info$reference_values, na.rm = TRUE) * 100
  )
  
  return(stats)
}