#' Fire Impact Detection
#' 
#' Detects and quantifies fire impact on vegetation FPAR
#' 
#' @author Your Name
#' @date 2025

source("scripts/utils/helper_functions.R")

#' Detect Fire Impact
#'
#' Quantifies fire impact as maximum FPAR decline during fire period
#'
#' @param burnt_ts Numeric vector of burnt pixel FPAR time series
#' @param unburnt_ts Numeric vector of unburnt reference FPAR time series
#' @param fire_indices List with start_idx and end_idx from convert_fire_year_to_indices()
#' @param pre_fire_baseline Number of observations for baseline calculation (default = 5)
#' @return List containing impact metrics
#' @export
#' @examples
#' fire_indices <- convert_fire_year_to_indices(200607)
#' impact <- detect_fire_impact(burnt_ts, unburnt_ts, fire_indices)
detect_fire_impact <- function(burnt_ts, unburnt_ts, fire_indices, pre_fire_baseline = 5) {
  
  # Validate inputs
  if (length(burnt_ts) == 0 || all(is.na(burnt_ts)) || 
      fire_indices$start_idx > length(burnt_ts)) {
    return(list(
      impact_magnitude = NA, 
      impact_date_idx = NA,
      pre_fire_mean = NA, 
      impact_value = NA,
      unburnt_reference = NA
    ))
  }
  
  fire_start_idx <- fire_indices$start_idx
  fire_end_idx <- min(fire_indices$end_idx, length(burnt_ts))
  
  # Calculate pre-fire baseline
  if (fire_start_idx <= pre_fire_baseline) {
    pre_fire_indices <- 1:(fire_start_idx-1)
  } else {
    pre_fire_indices <- (fire_start_idx-pre_fire_baseline):(fire_start_idx-1)
  }
  
  if (length(pre_fire_indices) > 0) {
    pre_fire_mean <- mean(burnt_ts[pre_fire_indices], na.rm = TRUE)
    unburnt_reference <- mean(unburnt_ts[pre_fire_indices], na.rm = TRUE)
  } else {
    pre_fire_mean <- burnt_ts[1]
    unburnt_reference <- unburnt_ts[1]
  }
  
  # Find maximum drop during fire period
  fire_period_ts <- burnt_ts[fire_start_idx:fire_end_idx]
  min_idx_relative <- which.min(fire_period_ts)
  
  if (length(min_idx_relative) == 0) {
    return(list(
      impact_magnitude = NA,
      impact_date_idx = NA,
      pre_fire_mean = pre_fire_mean,
      impact_value = NA,
      unburnt_reference = unburnt_reference
    ))
  }
  
  impact_idx_absolute <- fire_start_idx + min_idx_relative - 1
  impact_value <- fire_period_ts[min_idx_relative]
  impact_magnitude <- pre_fire_mean - impact_value
  
  return(list(
    impact_magnitude = max(0, impact_magnitude),  # Ensure non-negative
    impact_date_idx = impact_idx_absolute,
    pre_fire_mean = pre_fire_mean,
    impact_value = impact_value,
    unburnt_reference = unburnt_reference
  ))
}

#' Batch Process Fire Impact Detection
#'
#' Processes multiple burnt pixels for a cluster-fire year combination
#'
#' @param burnt_matrix Matrix of burnt pixel time series (rows = pixels, cols = time)
#' @param unburnt_ts Numeric vector of unburnt reference time series
#' @param fire_year Fire year in YYYYMM format
#' @param n_per_year Number of observations per year
#' @return Data frame with impact metrics for all pixels
#' @export
batch_detect_fire_impact <- function(burnt_matrix, unburnt_ts, fire_year, n_per_year = 46) {
  
  # Get fire indices
  fire_indices <- convert_fire_year_to_indices(fire_year, n_per_year)
  
  # Process each pixel
  results_list <- list()
  
  for (i in 1:nrow(burnt_matrix)) {
    burnt_ts <- as.numeric(burnt_matrix[i, ])
    
    # Skip if all NA
    if (all(is.na(burnt_ts))) {
      results_list[[i]] <- data.frame(
        pixel_id = i,
        impact_magnitude = NA,
        impact_date_idx = NA,
        pre_fire_mean = NA,
        impact_value = NA,
        unburnt_reference = NA
      )
      next
    }
    
    # Detect impact
    impact_result <- detect_fire_impact(burnt_ts, unburnt_ts, fire_indices)
    
    results_list[[i]] <- data.frame(
      pixel_id = i,
      impact_magnitude = impact_result$impact_magnitude,
      impact_date_idx = impact_result$impact_date_idx,
      pre_fire_mean = impact_result$pre_fire_mean,
      impact_value = impact_result$impact_value,
      unburnt_reference = impact_result$unburnt_reference
    )
  }
  
  # Combine results
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}