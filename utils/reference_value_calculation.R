# Functions to identify the typical peak-growth period for each phenology type
# and calculate long-term reference FPAR values from unburnt pixels.
#
# The input time series (unburnt_ts) should be the mean FPAR across all unburnt
# pixels in the same phenology type and fire year — not a single raw pixel.


# ------------------------------------------------------------------------------
# Find peak-growth periods and long-term reference values
# ------------------------------------------------------------------------------
# For each year in the unburnt reference time series, we find the 5-observation
# window (~40 days) with the highest mean FPAR. We then identify the most common
# timing of this peak across all years (the "modal peak position"). Reference
# FPAR values are collected from this same window position in every year and
# averaged to give the long-term reference level.
#
# The modal peak indices are also returned so that recovery can be checked only
# at these peak-growth observations rather than throughout the whole year.

find_reference_values <- function(unburnt_ts, n_per_year = 46, window_size = 5) {

  # If the series is shorter than one year, return a simple fallback
  if (length(unburnt_ts) < n_per_year) {
    fallback <- mean(unburnt_ts, na.rm = TRUE)
    return(list(
      reference_values      = rep(fallback, window_size),
      reference_indices     = seq_len(min(window_size, length(unburnt_ts))),
      modal_peak_position   = 1L,
      yearly_peak_indices   = integer(0),
      peak_position_in_year = 1L,
      n_years               = 0L,
      yearly_peaks          = list()
    ))
  }

  n_years <- floor(length(unburnt_ts) / n_per_year)

  # Step 1: find the best 5-observation window within each year
  yearly_peaks          <- vector("list", n_years)
  yearly_peak_positions <- integer(n_years)

  for (yr in seq_len(n_years)) {
    yr_start  <- (yr - 1L) * n_per_year + 1L
    yr_end    <- min(yr * n_per_year, length(unburnt_ts))
    yr_data   <- unburnt_ts[yr_start:yr_end]

    best_mean  <- -Inf
    best_start <- 1L
    max_start  <- length(yr_data) - window_size + 1L

    if (max_start >= 1L) {
      for (pos in seq_len(max_start)) {
        w <- yr_data[pos:(pos + window_size - 1L)]
        if (any(is.na(w))) next
        if (mean(w) > best_mean) {
          best_mean  <- mean(w)
          best_start <- pos
        }
      }
    }

    abs_start <- yr_start + best_start - 1L
    w_end     <- min(abs_start + window_size - 1L, length(unburnt_ts))

    yearly_peaks[[yr]] <- list(
      position       = best_start,
      absolute_start = abs_start,
      values         = unburnt_ts[abs_start:w_end],
      mean_value     = if (best_mean == -Inf) NA_real_ else best_mean
    )
    yearly_peak_positions[yr] <- best_start
  }

  # Step 2: find the most common peak position across all years
  pos_table      <- table(yearly_peak_positions)
  modal_position <- as.integer(names(pos_table)[which.max(pos_table)])

  # Step 3: collect FPAR values from that modal window in every year
  all_ref_values  <- numeric(0)
  all_ref_indices <- integer(0)

  for (yr in seq_len(n_years)) {
    yr_start  <- (yr - 1L) * n_per_year + 1L
    modal_abs <- yr_start + modal_position - 1L
    modal_end <- modal_abs + window_size - 1L
    if (modal_end > length(unburnt_ts)) next
    vals <- unburnt_ts[modal_abs:modal_end]
    if (any(is.na(vals))) next
    all_ref_values  <- c(all_ref_values,  vals)
    all_ref_indices <- c(all_ref_indices, modal_abs:modal_end)
  }

  # Fallback if every window had missing data
  if (length(all_ref_values) == 0) {
    all_ref_values  <- rep(mean(unburnt_ts, na.rm = TRUE), window_size)
    all_ref_indices <- seq_len(window_size)
  }

  # Step 4: collect all peak-period band indices (used to restrict recovery checks)
  yearly_peak_indices <- integer(0)
  for (yr in seq_len(n_years)) {
    yr_start  <- (yr - 1L) * n_per_year + 1L
    modal_abs <- yr_start + modal_position - 1L
    modal_end <- modal_abs + window_size - 1L
    valid_idx <- (modal_abs:modal_end)[(modal_abs:modal_end) <= length(unburnt_ts)]
    yearly_peak_indices <- c(yearly_peak_indices, valid_idx)
  }

  return(list(
    reference_values      = all_ref_values,
    reference_indices     = all_ref_indices,
    modal_peak_position   = modal_position,
    yearly_peak_indices   = yearly_peak_indices,
    peak_position_in_year = modal_position,
    n_years               = n_years,
    yearly_peaks          = yearly_peaks
  ))
}


# ------------------------------------------------------------------------------
# Simple summary statistics for reference values
# ------------------------------------------------------------------------------

calculate_reference_statistics <- function(reference_info) {
  v <- reference_info$reference_values
  data.frame(
    modal_peak_position = reference_info$modal_peak_position,
    n_years_used        = reference_info$n_years,
    n_reference_values  = length(v),
    reference_mean      = mean(v,   na.rm = TRUE),
    reference_sd        = sd(v,     na.rm = TRUE),
    reference_min       = min(v,    na.rm = TRUE),
    reference_max       = max(v,    na.rm = TRUE),
    reference_median    = median(v, na.rm = TRUE),
    cv_pct              = sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE) * 100
  )
}
