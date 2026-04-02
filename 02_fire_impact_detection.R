# Detect and measure fire impact on FPAR for each burnt pixel.
#
# For each burnt pixel within a fire year (1 July – 30 June) we:
#   1. Find when FPAR stopped declining and started to drop sharply
#      (fire-onset point, determined from day-to-day differences)
#   2. Find the date of the lowest FPAR value within the fire year (impact point)
#   3. Calculate the pre-fire average from the 5 observations before onset
#   4. Impact magnitude = pre-fire average minus the FPAR at impact
#      (set to zero if the minimum did not fall below the pre-fire average)


source("utils/helper_functions.R")


# ------------------------------------------------------------------------------
# Detect fire impact for a single pixel
# ------------------------------------------------------------------------------

detect_fire_impact <- function(burnt_ts, fire_indices, pre_fire_baseline = 5) {

  # Return NAs if the time series is empty or entirely missing
  na_result <- list(
    impact_magnitude = NA_real_,
    impact_date_idx  = NA_integer_,
    onset_date_idx   = NA_integer_,
    pre_fire_mean    = NA_real_,
    impact_value     = NA_real_
  )

  if (length(burnt_ts) == 0 || all(is.na(burnt_ts))) return(na_result)

  fire_start <- fire_indices$start_idx
  fire_end   <- min(fire_indices$end_idx, length(burnt_ts))

  if (fire_start > length(burnt_ts) || fire_start > fire_end) return(na_result)

  fire_period <- burnt_ts[fire_start:fire_end]

  # --- Impact point: observation with the lowest FPAR in the fire year -------
  min_rel <- which.min(fire_period)
  if (length(min_rel) == 0) return(na_result)

  impact_idx   <- fire_start + min_rel - 1L
  impact_value <- fire_period[min_rel]

  # --- Fire-onset point: last observation still declining before impact -------
  # We look at day-to-day differences in the series up to the impact point.
  # The onset is the last step where FPAR was still going down.
  onset_idx <- NA_integer_

  if (impact_idx > fire_start) {
    seg      <- burnt_ts[fire_start:impact_idx]
    diffs    <- diff(seg)             # negative value = FPAR declining
    dec_pos  <- which(diffs < 0)
    if (length(dec_pos) > 0) {
      onset_idx <- fire_start + max(dec_pos) - 1L
    }
  }

  # If no declining step was found, set onset one step before impact
  if (is.na(onset_idx)) onset_idx <- max(fire_start, impact_idx - 1L)

  # --- Pre-fire baseline: mean of 5 observations immediately before onset ----
  if (onset_idx <= 1L) {
    pre_indices <- integer(0)
  } else {
    pre_end     <- onset_idx - 1L
    pre_start   <- max(1L, pre_end - pre_fire_baseline + 1L)
    pre_indices <- pre_start:pre_end
  }

  pre_fire_mean <- if (length(pre_indices) > 0)
    mean(burnt_ts[pre_indices], na.rm = TRUE) else burnt_ts[1]

  # --- Impact magnitude -------------------------------------------------------
  impact_magnitude <- max(0, pre_fire_mean - impact_value)

  return(list(
    impact_magnitude = impact_magnitude,
    impact_date_idx  = impact_idx,
    onset_date_idx   = onset_idx,
    pre_fire_mean    = pre_fire_mean,
    impact_value     = impact_value
  ))
}


# ------------------------------------------------------------------------------
# Process all burnt pixels for one cluster–fire-year combination
# ------------------------------------------------------------------------------

batch_detect_fire_impact <- function(burnt_matrix, fire_year,
                                     n_per_year = 46, pre_fire_baseline = 5) {

  fire_indices <- convert_fire_year_to_indices(fire_year, n_per_year)
  n_pixels     <- nrow(burnt_matrix)
  results_list <- vector("list", n_pixels)

  for (i in seq_len(n_pixels)) {
    burnt_ts <- as.numeric(burnt_matrix[i, ])

    if (all(is.na(burnt_ts))) {
      results_list[[i]] <- data.frame(
        pixel_id         = i,
        impact_magnitude = NA_real_,
        impact_date_idx  = NA_integer_,
        onset_date_idx   = NA_integer_,
        pre_fire_mean    = NA_real_,
        impact_value     = NA_real_
      )
      next
    }

    res <- detect_fire_impact(burnt_ts, fire_indices, pre_fire_baseline)

    results_list[[i]] <- data.frame(
      pixel_id         = i,
      impact_magnitude = res$impact_magnitude,
      impact_date_idx  = res$impact_date_idx,
      onset_date_idx   = res$onset_date_idx,
      pre_fire_mean    = res$pre_fire_mean,
      impact_value     = res$impact_value
    )
  }

  do.call(rbind, results_list)
}
