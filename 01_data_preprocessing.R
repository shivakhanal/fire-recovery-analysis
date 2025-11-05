#' Data Preprocessing Script
#' 
#' Prepares input data for fire recovery analysis
#' 
#' @author Your Name
#' @date 2025

library(terra)
library(sf)
library(dplyr)

source("scripts/utils/helper_functions.R")

#' Load and Prepare FPAR Stack
#'
#' Loads MODIS FPAR raster stack and performs quality checks
#'
#' @param fpar_file Path to FPAR raster stack
#' @param start_year Start year for subsetting (optional)
#' @param end_year End year for subsetting (optional)
#' @return SpatRaster object
#' @export
load_fpar_stack <- function(fpar_file, start_year = NULL, end_year = NULL) {
  
  cat("Loading FPAR raster stack from:", fpar_file, "\n")
  
  # Load raster
  fpar_stack <- rast(fpar_file)
  
  cat("Loaded stack with", nlyr(fpar_stack), "layers\n")
  cat("CRS:", crs(fpar_stack, describe = TRUE)$name, "\n")
  cat("Extent:", paste(as.vector(ext(fpar_stack)), collapse = ", "), "\n")
  cat("Resolution:", paste(res(fpar_stack), collapse = " x "), "\n")
  
  # Subset by year if requested
  if (!is.null(start_year) && !is.null(end_year)) {
    layer_names <- names(fpar_stack)
    dates <- extract_dates_from_names(layer_names)
    
    keep_layers <- which(
      year(dates) >= start_year & year(dates) <= end_year
    )
    
    if (length(keep_layers) > 0) {
      fpar_stack <- fpar_stack[[keep_layers]]
      cat("Subset to", nlyr(fpar_stack), "layers from", start_year, "to", end_year, "\n")
    }
  }
  
  return(fpar_stack)
}

#' Load and Prepare Cluster Data
#'
#' Loads cluster assignments and coordinates
#'
#' @param cluster_file Path to cluster data CSV
#' @param fpar_stack SpatRaster to match CRS
#' @return Data frame with cluster assignments
#' @export
load_cluster_data <- function(cluster_file, fpar_stack = NULL) {
  
  cat("Loading cluster data from:", cluster_file, "\n")
  
  cluster_data <- read.csv(cluster_file)
  
  # Check required columns
  required_cols <- c("class", "x", "y")
  if (!all(required_cols %in% names(cluster_data))) {
    stop("Cluster data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  cat("Loaded", nrow(cluster_data), "pixels with", 
      length(unique(cluster_data$class)), "clusters\n")
  
  # Print cluster size distribution
  cluster_sizes <- table(cluster_data$class)
  cat("Cluster size range:", min(cluster_sizes), "-", max(cluster_sizes), "pixels\n")
  
  return(cluster_data)
}

#' Load and Prepare Fire Polygons
#'
#' Loads fire history polygons and performs validation
#'
#' @param fire_file Path to fire polygon file (GeoPackage or Shapefile)
#' @param fpar_stack Optional SpatRaster to check CRS match
#' @param start_year Filter fires from this year onwards (optional)
#' @param end_year Filter fires up to this year (optional)
#' @return SF object with fire polygons
#' @export
load_fire_polygons <- function(fire_file, fpar_stack = NULL, 
                               start_year = NULL, end_year = NULL) {
  
  cat("Loading fire polygons from:", fire_file, "\n")
  
  fire_polygons <- st_read(fire_file, quiet = TRUE)
  
  # Check required columns
  if (!("FireYear" %in% names(fire_polygons))) {
    stop("Fire polygons must contain 'FireYear' column")
  }
  
  # Filter by year if requested
  if (!is.null(start_year)) {
    fire_polygons <- fire_polygons %>%
      filter(FireYear >= start_year * 100)  # Assumes YYYYMM format
  }
  
  if (!is.null(end_year)) {
    fire_polygons <- fire_polygons %>%
      filter(FireYear <= (end_year + 1) * 100)
  }
  
  cat("Loaded", nrow(fire_polygons), "fire polygons\n")
  cat("Fire years:", paste(sort(unique(fire_polygons$FireYear)), collapse = ", "), "\n")
  cat("CRS:", st_crs(fire_polygons)$input, "\n")
  
  # Check CRS match with raster if provided
  if (!is.null(fpar_stack)) {
    if (st_crs(fire_polygons)$input != crs(fpar_stack, describe = TRUE)$code) {
      cat("WARNING: Fire polygon CRS does not match raster CRS\n")
      cat("  Fire CRS:", st_crs(fire_polygons)$input, "\n")
      cat("  Raster CRS:", crs(fpar_stack, describe = TRUE)$name, "\n")
      cat("  Will transform fire polygons during analysis\n")
    }
  }
  
  return(fire_polygons)
}

#' Load Unburnt Reference Data
#'
#' Loads pre-calculated mean FPAR for unburnt pixels by cluster and fire year
#'
#' @param reference_file Path to unburnt reference CSV
#' @return Data frame with unburnt reference time series
#' @export
load_unburnt_reference <- function(reference_file) {
  
  cat("Loading unburnt reference data from:", reference_file, "\n")
  
  unburnt_ref <- read.csv(reference_file)
  
  # Check required columns
  required_cols <- c("cluster", "fireyear")
  if (!all(required_cols %in% names(unburnt_ref))) {
    stop("Unburnt reference must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Check for time series columns
  ts_cols <- grep("^BU_LAI_FPAR_|^X[0-9]", names(unburnt_ref), value = TRUE)
  
  if (length(ts_cols) == 0) {
    stop("Unburnt reference must contain time series columns")
  }
  
  cat("Loaded reference data for", nrow(unburnt_ref), "cluster-fire year combinations\n")
  cat("Time series length:", length(ts_cols), "observations\n")
  
  return(unburnt_ref)
}

#' Prepare All Input Data
#'
#' Convenience function to load and validate all required data
#'
#' @param fpar_file Path to FPAR raster stack
#' @param cluster_file Path to cluster data
#' @param fire_file Path to fire polygons
#' @param reference_file Path to unburnt reference data
#' @param start_year Optional start year
#' @param end_year Optional end year
#' @return List containing all input data objects
#' @export
prepare_all_data <- function(fpar_file, cluster_file, fire_file, reference_file,
                             start_year = NULL, end_year = NULL) {
  
  cat("\n=== LOADING INPUT DATA ===\n\n")
  
  # Load data
  fpar_stack <- load_fpar_stack(fpar_file, start_year, end_year)
  cluster_data <- load_cluster_data(cluster_file, fpar_stack)
  fire_polygons <- load_fire_polygons(fire_file, fpar_stack, start_year, end_year)
  unburnt_ref <- load_unburnt_reference(reference_file)
  
  # Validate
  cat("\n=== VALIDATING INPUT DATA ===\n\n")
  validation <- validate_input_data(cluster_data, fire_polygons, fpar_stack, unburnt_ref)
  
  if (!validation$valid) {
    stop("Input data validation failed")
  }
  
  cat("\n=== DATA LOADING COMPLETE ===\n\n")
  
  return(list(
    fpar_stack = fpar_stack,
    cluster_data = cluster_data,
    fire_polygons = fire_polygons,
    unburnt_reference = unburnt_ref,
    validation = validation
  ))
}