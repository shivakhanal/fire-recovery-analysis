# Helper functions used across the fire recovery pipeline.
# Author: Shiva Khanal, 2025

library(terra)
library(sf)
library(dplyr)
library(stringr)


# ------------------------------------------------------------------------------
# Convert a fire year code to position indices in the FPAR time series
# ------------------------------------------------------------------------------
# Australian fire years run from 1 July to 30 June.
# Fire years are stored as a 6-digit code: YYYYMM, where YYYY is the start
# calendar year and MM is the last two digits of the end year.
# For example, 200607 means 1 July 2006 to 30 June 2007.
#
# The MODIS FPAR stack has 1,012 bands (8-day composites, 2001-2022).
# This function returns the band numbers for 1 July and 30 June of a given
# fire year so the rest of the code knows where to look.

convert_fire_year_to_indices <- function(fire_year, n_per_year = 46,
                                         stack_start_year = 2001) {

  fire_year_str   <- sprintf("%06d", fire_year)
  start_year      <- as.numeric(substr(fire_year_str, 1, 4))
  end_year_suffix <- substr(fire_year_str, 5, 6)

  century  <- floor(start_year / 100) * 100
  end_year <- century + as.numeric(end_year_suffix)
  if (end_year < start_year) end_year <- end_year + 100

  # 1 July is approximately the 23rd 8-day step of the year (day of year ~182)
  july1_step  <- ceiling(182 / 8)
  june30_step <- ceiling(181 / 8)

  yrs_start <- start_year - stack_start_year
  yrs_end   <- end_year   - stack_start_year

  fire_start_idx <- yrs_start * n_per_year + july1_step
  fire_end_idx   <- yrs_end   * n_per_year + june30_step

  return(list(
    start_idx  = max(1, fire_start_idx),
    end_idx    = fire_end_idx,
    start_year = start_year,
    end_year   = end_year
  ))
}


# ------------------------------------------------------------------------------
# Extract dates from MODIS FPAR raster band names
# ------------------------------------------------------------------------------
# MODIS bands are named in the format BU_LAI_FPAR_YYYY_YYYYJJJ_b1 where JJJ is
# the day of year. This function converts those names into R Date objects so
# results can be plotted against real dates.

extract_dates_from_names <- function(column_names) {

  cat("Extracting dates from", length(column_names), "band names\n")

  date_parts <- NULL

  # Try format: YYYY_YYYYJJJ (e.g. 2001_2001009)
  date_parts <- stringr::str_extract(column_names, "\\d{4}_\\d{7}")
  if (all(is.na(date_parts))) {
    # Try bare YYYYJJJ (e.g. 2001009)
    date_parts <- stringr::str_extract(column_names, "\\d{7}")
  }

  dates <- sapply(seq_along(date_parts), function(i) {
    x <- date_parts[i]
    if (is.na(x)) return(NA_character_)
    tryCatch({
      if (grepl("_", x) && nchar(x) > 8) {
        parts <- strsplit(x, "_")[[1]]
        year  <- as.numeric(substr(parts[2], 1, 4))
        doy   <- as.numeric(substr(parts[2], 5, 7))
      } else {
        year  <- as.numeric(substr(x, 1, 4))
        doy   <- as.numeric(substr(x, 5, 7))
      }
      as.character(as.Date(doy - 1, origin = paste0(year, "-01-01")))
    }, error = function(e) NA_character_)
  })

  result_dates <- as.Date(dates)
  valid_n <- sum(!is.na(result_dates))
  cat("Converted", valid_n, "of", length(dates), "dates successfully\n")

  if (valid_n == 0) {
    cat("Could not parse dates — using 8-day sequence from 2001-01-01\n")
    result_dates <- seq(as.Date("2001-01-01"), by = 8,
                        length.out = length(column_names))
  }

  return(result_dates)
}


# ------------------------------------------------------------------------------
# Check that all input data objects have the expected structure
# ------------------------------------------------------------------------------

validate_input_data <- function(cluster_data, fire_polygons,
                                fpar_stack, unburnt_reference) {

  ok  <- TRUE
  msg <- character(0)

  if (!all(c("class", "x", "y") %in% names(cluster_data))) {
    ok  <- FALSE
    msg <- c(msg, "cluster_data needs columns: class, x, y")
  }
  if (!inherits(fire_polygons, "sf")) {
    ok  <- FALSE
    msg <- c(msg, "fire_polygons must be an sf object")
  }
  if (!("FireYear" %in% names(fire_polygons))) {
    ok  <- FALSE
    msg <- c(msg, "fire_polygons must have a FireYear column (format YYYYMM)")
  }
  if (!inherits(fpar_stack, "SpatRaster")) {
    ok  <- FALSE
    msg <- c(msg, "fpar_stack must be a SpatRaster")
  }
  if (!all(c("cluster", "fireyear") %in% names(unburnt_reference))) {
    ok  <- FALSE
    msg <- c(msg, "unburnt_reference needs columns: cluster, fireyear")
  }

  ts_cols <- grep("^BU_LAI_FPAR_|^X[0-9]", names(unburnt_reference), value = TRUE)
  if (length(ts_cols) == 0) {
    ok  <- FALSE
    msg <- c(msg, "unburnt_reference must contain FPAR time series columns")
  }

  if (ok) {
    cat("All inputs look good\n")
    cat("  Phenology types:", length(unique(cluster_data$class)), "\n")
    cat("  Fire years:     ", length(unique(fire_polygons$FireYear)), "\n")
    cat("  FPAR bands:     ", nlyr(fpar_stack), "\n")
  } else {
    cat("Input problems found:\n")
    for (m in msg) cat(" -", m, "\n")
  }

  return(list(valid = ok, messages = msg))
}


# ------------------------------------------------------------------------------
# Create output folder structure
# ------------------------------------------------------------------------------

create_output_directories <- function(base_path = "outputs") {

  dirs <- list(
    base        = base_path,
    rasters     = file.path(base_path, "rasters"),
    impact      = file.path(base_path, "rasters", "impact"),
    recovery    = file.path(base_path, "rasters", "recovery"),
    plots       = file.path(base_path, "plots"),
    timeseries  = file.path(base_path, "plots", "timeseries"),
    summary     = file.path(base_path, "plots", "summary"),
    diagnostics = file.path(base_path, "diagnostics"),
    tables      = file.path(base_path, "tables")
  )

  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  cat("Output folders created under:", base_path, "\n")
  return(dirs)
}


# ------------------------------------------------------------------------------
# Write a timestamped line to a log file
# ------------------------------------------------------------------------------

log_message <- function(message, log_file, level = "INFO") {
  entry <- sprintf("[%s] %s: %s\n",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, message)
  cat(entry)
  write(entry, file = log_file, append = TRUE)
}
