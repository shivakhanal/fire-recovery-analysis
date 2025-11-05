#' Post-Fire Recovery Estimation
#' 
#' Two-stage recovery detection with peak-period restriction
#' 
#' @author Your Name
#' @date 2025

source("scripts/utils/helper_functions.R")
source("scripts/utils/reference_value_calculation.R")

#' Calculate Recovery Time with Two-Stage Detection
#'
#' Implements hierarchical recovery detection with peak-period restriction
#'
#' @param burnt_ts Numeric vector of burnt pixel FPAR time series
#' @param unburnt_ts Numeric vector of unburnt reference FPAR time series
#' @param fire_impact_result Output from detect_fire_impact()
#' @param recovery_threshold_stage1 Stage 1 threshold (proportion of unburnt, default = 0.95)
#' @param recovery_threshold_stage2 Stage 2 threshold (proportion of reference mean, default = 0.90)
#' @param n_per_year Number of observations per year (default = 46)
#' @param strict_peak_only If TRUE, only detect recovery during peak periods (default = TRUE)
#' @return List containing recovery metrics
#' @export
#' @examples
#' recovery <- calculate_recovery_time(burnt_ts, unburnt_ts, fire_impact, 
#'                                     recovery_threshold_stage1 = 0.95,
#'                                     recovery_threshold_stage2 = 0.90,
#'                                     strict_peak_only = TRUE)
calculate_recovery_time <- function(burnt_ts, unburnt_ts, fire_impact_result, 
                                    recovery_threshold_stage1 = 0.95,
                                    recovery_threshold_stage2 = 0.90,
                                    n_per_year = 46, 
                                    strict_peak_only = TRUE) {
  
  # Check if impact was detected
  if (is.na(fire_impact_result$impact_date_idx)) {
    return(list(
      recovery_days = -1, 
      recovery_idx = -1, 
      recovered = FALSE,
      reference_values = NA, 
      reference_mean = NA, 
      recovery_method = "no_impact"
    ))
  }
  
  impact_idx <- fire_impact_result$impact_date_idx
  
  # Get reference values using modal peak method
  reference_info <- find_reference_values(unburnt_ts, n_per_year = n_per_year, window_size = 5)
  reference_mean <- mean(reference_info$reference_values, na.rm = TRUE)
  
  # Get all yearly peak indices that occur after the impact
  future_peak_indices <- reference_info$yearly_peak_indices[
    reference_info$yearly_peak_indices > impact_idx
  ]
  
  if (length(future_peak_indices) == 0) {
    return(list(
      recovery_days = -1, 
      recovery_idx = -1, 
      recovered = FALSE,
      reference_values = reference_info$reference_values,
      reference_mean = reference_mean,
      recovery_method = "no_future_peaks"
    ))
  }
  
  # STAGE 1: Dynamic Catchup - Check if burnt reaches 95% of contemporary unburnt
  recovery_idx_stage1 <- NA
  stage1_checks <- data.frame(
    peak_idx = integer(),
    burnt_value = numeric(),
    unburnt_value = numeric(),
    catchup_threshold = numeric(),
    meets_criterion = logical()
  )
  
  for (peak_idx in future_peak_indices) {
    if (peak_idx <= length(burnt_ts) && peak_idx <= length(unburnt_ts) && 
        !is.na(burnt_ts[peak_idx]) && !is.na(unburnt_ts[peak_idx])) {
      
      burnt_value <- burnt_ts[peak_idx]
      unburnt_value <- unburnt_ts[peak_idx]
      catchup_threshold <- recovery_threshold_stage1 * unburnt_value
      meets_criterion <- burnt_value >= catchup_threshold
      
      # Store check results
      stage1_checks <- rbind(stage1_checks, data.frame(
        peak_idx = peak_idx,
        burnt_value = burnt_value,
        unburnt_value = unburnt_value,
        catchup_threshold = catchup_threshold,
        meets_criterion = meets_criterion
      ))
      
      if (meets_criterion) {
        recovery_idx_stage1 <- peak_idx
        break
      }
    }
  }
  
  # If Stage 1 successful, return
  if (!is.na(recovery_idx_stage1)) {
    recovery_days <- (recovery_idx_stage1 - impact_idx) * 8
    
    return(list(
      recovery_days = recovery_days,
      recovery_idx = recovery_idx_stage1,
      recovered = TRUE,
      reference_values = reference_info$reference_values,
      reference_mean = reference_mean,
      modal_peak_position = reference_info$modal_peak_position,
      recovery_method = "stage1_dynamic_catchup",
      future_peak_indices = future_peak_indices,
      stage1_checks = stage1_checks,
      stage2_checks = NULL
    ))
  }
  
  # STAGE 2: Static Threshold - Check if burnt reaches 90% of reference mean
  if (!strict_peak_only) {
    # This section kept for backwards compatibility but typically not used
    recovery_idx_stage2 <- NA
    static_threshold_value <- recovery_threshold_stage2 * reference_mean
    
    for (peak_idx in future_peak_indices) {
      if (peak_idx <= length(burnt_ts) && !is.na(burnt_ts[peak_idx])) {
        burnt_value <- burnt_ts[peak_idx]
        
        if (burnt_value >= static_threshold_value) {
          recovery_idx_stage2 <- peak_idx
          break
        }
      }
    }
    
    if (!is.na(recovery_idx_stage2)) {
      recovery_days <- (recovery_idx_stage2 - impact_idx) * 8
      
      return(list(
        recovery_days = recovery_days,
        recovery_idx = recovery_idx_stage2,
        recovered = TRUE,
        reference_values = reference_info$reference_values,
        reference_mean = reference_mean,
        modal_peak_position = reference_info$modal_peak_position,
        recovery_method = "stage2_static_threshold",
        future_peak_indices = future_peak_indices,
        stage1_checks = stage1_checks,
        stage2_checks = data.frame(
          static_threshold = static_threshold_value,
          recovery_idx = recovery_idx_stage2
        )
      ))
    }
  }
  
  # No recovery detected
  max_peak_value <- ifelse(nrow(stage1_checks) > 0, 
                           max(stage1_checks$burnt_value, na.rm = TRUE), 
                           NA)
  
  return(list(
    recovery_days = -1,
    recovery_idx = -1,
    recovered = FALSE,
    reference_values = reference_info$reference_values,
    reference_mean = reference_mean,
    modal_peak_position = reference_info$modal_peak_position,
    recovery_method = "no_recovery",
    future_peak_indices = future_peak_indices,
    stage1_checks = stage1_checks,
    stage2_checks = NULL,
    max_peak_value = max_peak_value
  ))
}

#' Batch Process Recovery Estimation
#'
#' Processes multiple burnt pixels for recovery estimation
#'
#' @param burnt_matrix Matrix of burnt pixel time series (rows = pixels, cols = time)
#' @param unburnt_ts Numeric vector of unburnt reference time series
#' @param impact_results Data frame from batch_detect_fire_impact()
#' @param recovery_threshold_stage1 Stage 1 threshold (default = 0.95)
#' @param recovery_threshold_stage2 Stage 2 threshold (default = 0.90)
#' @param n_per_year Number of observations per year
#' @return Data frame with recovery metrics for all pixels
#' @export
batch_calculate_recovery <- function(burnt_matrix, unburnt_ts, impact_results, 
                                     recovery_threshold_stage1 = 0.95,
                                     recovery_threshold_stage2 = 0.90,
                                     n_per_year = 46) {
  
  results_list <- list()
  
  for (i in 1:nrow(burnt_matrix)) {
    burnt_ts <- as.numeric(burnt_matrix[i, ])
    
    # Get impact result for this pixel
    impact_result <- list(
      impact_magnitude = impact_results$impact_magnitude[i],
      impact_date_idx = impact_results$impact_date_idx[i],
      pre_fire_mean = impact_results$pre_fire_mean[i],
      impact_value = impact_results$impact_value[i],
      unburnt_reference = impact_results$unburnt_reference[i]
    )
    
    # Calculate recovery
    recovery_result <- calculate_recovery_time(
      burnt_ts, unburnt_ts, impact_result,
      recovery_threshold_stage1 = recovery_threshold_stage1,
      recovery_threshold_stage2 = recovery_threshold_stage2,
      n_per_year = n_per_year,
      strict_peak_only = TRUE
    )
    
    results_list[[i]] <- data.frame(
      pixel_id = i,
      recovery_days = recovery_result$recovery_days,
      recovery_idx = recovery_result$recovery_idx,
      recovered = recovery_result$recovered,
      reference_mean = recovery_result$reference_mean,
      modal_peak_position = ifelse(is.null(recovery_result$modal_peak_position), 
                                   NA, recovery_result$modal_peak_position),
      recovery_method = recovery_result$recovery_method
    )
  }
  
  # Combine results
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}