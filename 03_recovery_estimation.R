# Estimate post-fire canopy recovery time for each burnt pixel.
#
# Recovery is checked only during the annual peak-growth period of each
# phenology type, not throughout the full year. Two criteria are applied:
#
#   Stage 1 (primary): the burnt pixel's FPAR reaches at least 95% of the
#     current-year unburnt average at the same peak-growth observation.
#     This accounts for year-to-year variation in climate conditions.
#
#   Stage 2 (fallback): if Stage 1 is never met, the burnt pixel's FPAR
#     reaches at least 90% of the long-term average peak FPAR for that
#     phenology type. This captures substantial but incomplete recovery for
#     pixels severely affected by high-severity fire.
#
# Recovery time is calculated as:
#   (recovery observation number - impact observation number) x 8 days
#
# Pixels meeting neither criterion by the end of 2022 are returned with

source("utils/helper_functions.R")
source("utils/reference_value_calculation.R")

# ------------------------------------------------------------------------------
# Estimate recovery time for a single pixel
# ------------------------------------------------------------------------------
# unburnt_ts must be the mean FPAR across all unburnt pixels in the same
# phenology type and fire year — not a single raw pixel time series.

calculate_recovery_time <- function(burnt_ts,
                                    unburnt_ts,
                                    fire_impact_result,
                                    recovery_threshold_stage1 = 0.95,
                                    recovery_threshold_stage2 = 0.90,
                                    n_per_year = 46) {

  # No impact found — cannot assess recovery
  if (is.na(fire_impact_result$impact_date_idx)) {
    return(list(
      recovery_days       = -1L,
      recovery_idx        = -1L,
      recovered           = FALSE,
      reference_mean      = NA_real_,
      recovery_method     = "no_impact",
      modal_peak_position = NA_integer_,
      future_peak_indices = integer(0),
      stage1_checks       = NULL,
      stage2_checks       = NULL
    ))
  }

  impact_idx <- fire_impact_result$impact_date_idx

  # Identify the typical peak-growth period and long-term reference level
  ref_info       <- find_reference_values(unburnt_ts, n_per_year = n_per_year,
                                          window_size = 5)
  reference_mean <- mean(ref_info$reference_values, na.rm = TRUE)

  # Only check recovery at peak-growth observations after the impact date
  future_peaks <- ref_info$yearly_peak_indices[
    ref_info$yearly_peak_indices > impact_idx
  ]

  if (length(future_peaks) == 0) {
    return(list(
      recovery_days       = -1L,
      recovery_idx        = -1L,
      recovered           = FALSE,
      reference_mean      = reference_mean,
      recovery_method     = "no_future_peaks",
      modal_peak_position = ref_info$modal_peak_position,
      future_peak_indices = integer(0),
      stage1_checks       = NULL,
      stage2_checks       = NULL
    ))
  }

  # --- Stage 1: does the burnt pixel reach 95% of the current-year unburnt? -
  stage1_checks   <- data.frame(peak_idx = integer(0), burnt_value = numeric(0),
                                unburnt_value = numeric(0),
                                threshold = numeric(0), met = logical(0))
  recovery_idx_s1 <- NA_integer_

  for (pk in future_peaks) {
    if (pk > length(burnt_ts) || pk > length(unburnt_ts)) next
    b <- burnt_ts[pk];  u <- unburnt_ts[pk]
    if (is.na(b) || is.na(u)) next
    thr <- recovery_threshold_stage1 * u
    met <- b >= thr
    stage1_checks <- rbind(stage1_checks,
                           data.frame(peak_idx = pk, burnt_value = b,
                                      unburnt_value = u, threshold = thr, met = met))
    if (met) { recovery_idx_s1 <- pk; break }
  }

  if (!is.na(recovery_idx_s1)) {
    return(list(
      recovery_days       = (recovery_idx_s1 - impact_idx) * 8L,
      recovery_idx        = recovery_idx_s1,
      recovered           = TRUE,
      reference_mean      = reference_mean,
      recovery_method     = "stage1_dynamic_catchup",
      modal_peak_position = ref_info$modal_peak_position,
      future_peak_indices = future_peaks,
      stage1_checks       = stage1_checks,
      stage2_checks       = NULL
    ))
  }

  # --- Stage 2: does the burnt pixel reach 90% of the long-term average? ----
  static_thr    <- recovery_threshold_stage2 * reference_mean
  stage2_checks <- data.frame(peak_idx = integer(0), burnt_value = numeric(0),
                              static_threshold = numeric(0), met = logical(0))
  recovery_idx_s2 <- NA_integer_

  for (pk in future_peaks) {
    if (pk > length(burnt_ts)) next
    b <- burnt_ts[pk]
    if (is.na(b)) next
    met <- b >= static_thr
    stage2_checks <- rbind(stage2_checks,
                           data.frame(peak_idx = pk, burnt_value = b,
                                      static_threshold = static_thr, met = met))
    if (met) { recovery_idx_s2 <- pk; break }
  }

  if (!is.na(recovery_idx_s2)) {
    return(list(
      recovery_days       = (recovery_idx_s2 - impact_idx) * 8L,
      recovery_idx        = recovery_idx_s2,
      recovered           = TRUE,
      reference_mean      = reference_mean,
      recovery_method     = "stage2_static_threshold",
      modal_peak_position = ref_info$modal_peak_position,
      future_peak_indices = future_peaks,
      stage1_checks       = stage1_checks,
      stage2_checks       = stage2_checks
    ))
  }

  # --- Neither criterion met -------------------------------------------------
  return(list(
    recovery_days       = -1L,
    recovery_idx        = -1L,
    recovered           = FALSE,
    reference_mean      = reference_mean,
    recovery_method     = "no_recovery",
    modal_peak_position = ref_info$modal_peak_position,
    future_peak_indices = future_peaks,
    stage1_checks       = stage1_checks,
    stage2_checks       = stage2_checks,
    max_peak_value      = if (nrow(stage1_checks) > 0)
                            max(stage1_checks$burnt_value, na.rm = TRUE) else NA_real_
  ))
}


# ------------------------------------------------------------------------------
# Process all burnt pixels for one cluster–fire-year combination
# ------------------------------------------------------------------------------

batch_calculate_recovery <- function(burnt_matrix,
                                     unburnt_ts,
                                     impact_results,
                                     recovery_threshold_stage1 = 0.95,
                                     recovery_threshold_stage2 = 0.90,
                                     n_per_year = 46) {

  n_pixels     <- nrow(burnt_matrix)
  results_list <- vector("list", n_pixels)

  for (i in seq_len(n_pixels)) {
    burnt_ts <- as.numeric(burnt_matrix[i, ])

    impact <- list(
      impact_magnitude = impact_results$impact_magnitude[i],
      impact_date_idx  = impact_results$impact_date_idx[i],
      onset_date_idx   = impact_results$onset_date_idx[i],
      pre_fire_mean    = impact_results$pre_fire_mean[i],
      impact_value     = impact_results$impact_value[i]
    )

    rec <- calculate_recovery_time(
      burnt_ts                  = burnt_ts,
      unburnt_ts                = unburnt_ts,
      fire_impact_result        = impact,
      recovery_threshold_stage1 = recovery_threshold_stage1,
      recovery_threshold_stage2 = recovery_threshold_stage2,
      n_per_year                = n_per_year
    )

    results_list[[i]] <- data.frame(
      pixel_id            = i,
      recovery_days       = rec$recovery_days,
      recovery_idx        = rec$recovery_idx,
      recovered           = rec$recovered,
      reference_mean      = rec$reference_mean,
      modal_peak_position = ifelse(is.null(rec$modal_peak_position),
                                   NA_integer_, rec$modal_peak_position),
      recovery_method     = rec$recovery_method
    )
  }

  do.call(rbind, results_list)
}
