#' Example Workflow for Fire Recovery Analysis
#' 
#' Demonstrates complete analysis pipeline
#' 
#' @author Your Name
#' @date 2025

# Load required libraries
library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# Source all functions
source("scripts/01_data_preprocessing.R")
source("scripts/02_fire_impact_detection.R")
source("scripts/03_recovery_estimation.R")
source("scripts/04_visualization.R")
source("scripts/utils/helper_functions.R")
source("scripts/utils/reference_value_calculation.R")

# ============================================================================
# STEP 1: SETUP
# ============================================================================

# Create output directories
output_dirs <- create_output_directories("outputs")

# Initialize log file
log_file <- file.path(output_dirs$diagnostics, "analysis_log.txt")
log_message("Starting fire recovery analysis", log_file, "INFO")

# ============================================================================
# STEP 2: LOAD AND PREPARE DATA
# ============================================================================

cat("\n=== STEP 2: LOADING DATA ===\n")

# Load all input data
input_data <- prepare_all_data(
  fpar_file = "data/modis_fpar_stack_2000_2022.tif",
  cluster_file = "data/cluster_assignments.csv",
  fire_file = "data/fire_history.gpkg",
  reference_file = "data/unburnt_mean_timeseries.csv",
  start_year = 2000,
  end_year = 2022
)

# Extract objects
fpar_stack <- input_data$fpar_stack
cluster_data <- input_data$cluster_data
fire_polygons <- input_data$fire_polygons
unburnt_ref <- input_data$unburnt_reference

log_message("Data loading completed successfully", log_file, "INFO")

# ============================================================================
# STEP 3: EXAMPLE ANALYSIS FOR SINGLE CLUSTER-FIRE YEAR
# ============================================================================

cat("\n=== STEP 3: EXAMPLE ANALYSIS (Cluster 23, Fire Year 200607) ===\n")

# Select example cluster and fire year
example_cluster <- 23
example_fire_year <- 200607

# Get unburnt reference for this combination
unburnt_data <- unburnt_ref %>%
  filter(cluster == example_cluster, fireyear == example_fire_year)

if (nrow(unburnt_data) == 0) {
  stop("No unburnt reference data found for example cluster/fire year")
}

# Extract time series columns
ts_cols <- grep("^BU_LAI_FPAR_|^X[0-9]", names(unburnt_data), value = TRUE)
unburnt_ts <- as.numeric(unburnt_data[1, ts_cols])

cat("Unburnt reference time series length:", length(unburnt_ts), "\n")

# Get cluster points
cluster_points <- cluster_data %>%
  filter(class == example_cluster)

cat("Cluster", example_cluster, "has", nrow(cluster_points), "pixels\n")

# Convert to sf object
cluster_points_sf <- st_as_sf(cluster_points, 
                              coords = c("x", "y"), 
                              crs = crs(fpar_stack))

# Get fire polygons for this fire year
fire_polys <- fire_polygons %>%
  filter(FireYear == example_fire_year)

cat("Found", nrow(fire_polys), "fire polygons for fire year", example_fire_year, "\n")

# Ensure same CRS
if (st_crs(cluster_points_sf) != st_crs(fire_polys)) {
  cluster_points_sf <- st_transform(cluster_points_sf, st_crs(fire_polys))
}

# Find burnt points
burnt_intersection <- st_intersects(cluster_points_sf, fire_polys)
burnt_points_idx <- which(lengths(burnt_intersection) > 0)

cat("Found", length(burnt_points_idx), "burnt pixels\n")

if (length(burnt_points_idx) == 0) {
  stop("No burnt pixels found for example")
}

# Extract burnt pixel data
burnt_points_sf <- cluster_points_sf[burnt_points_idx, ]
original_coords <- st_coordinates(burnt_points_sf)

# Transform back to raster CRS for extraction
burnt_points_sf <- st_transform(burnt_points_sf, crs(fpar_stack))

# Extract time series from raster
burnt_vect <- vect(burnt_points_sf)
extracted_values <- terra::extract(fpar_stack, burnt_vect)

# Remove ID column
if ("ID" %in% names(extracted_values)) {
  extracted_values <- extracted_values[, -which(names(extracted_values) == "ID")]
}

cat("Extracted time series for", nrow(extracted_values), "burnt pixels\n")

# ============================================================================
# STEP 4: DETECT FIRE IMPACT
# ============================================================================

cat("\n=== STEP 4: DETECTING FIRE IMPACT ===\n")

# Batch detect fire impact
impact_results <- batch_detect_fire_impact(
  burnt_matrix = as.matrix(extracted_values),
  unburnt_ts = unburnt_ts,
  fire_year = example_fire_year,
  n_per_year = 46
)

# Summary statistics
cat("\nFire Impact Summary:\n")
cat("  Mean impact:", round(mean(impact_results$impact_magnitude, na.rm = TRUE), 3), "\n")
cat("  Median impact:", round(median(impact_results$impact_magnitude, na.rm = TRUE), 3), "\n")
cat("  SD impact:", round(sd(impact_results$impact_magnitude, na.rm = TRUE), 3), "\n")
cat("  Pixels with valid impact:", sum(!is.na(impact_results$impact_magnitude)), "\n")

log_message(sprintf("Fire impact detected for %d pixels", 
                    sum(!is.na(impact_results$impact_magnitude))), 
            log_file, "INFO")

# ============================================================================
# STEP 5: CALCULATE RECOVERY TIME
# ============================================================================

cat("\n=== STEP 5: CALCULATING RECOVERY TIME ===\n")

# Batch calculate recovery
recovery_results <- batch_calculate_recovery(
  burnt_matrix = as.matrix(extracted_values),
  unburnt_ts = unburnt_ts,
  impact_results = impact_results,
  recovery_threshold_stage1 = 0.95,  # Stage 1: 95% of unburnt
  recovery_threshold_stage2 = 0.90,  # Stage 2: 90% of reference mean
  n_per_year = 46
)

# Summary statistics
recovered_pixels <- recovery_results %>%
  filter(recovered == TRUE)

cat("\nRecovery Summary:\n")
cat("  Total pixels analyzed:", nrow(recovery_results), "\n")
cat("  Pixels recovered:", nrow(recovered_pixels), "\n")
cat("  Recovery rate:", round(nrow(recovered_pixels)/nrow(recovery_results)*100, 1), "%\n")

if (nrow(recovered_pixels) > 0) {
  cat("  Mean recovery time:", round(mean(recovered_pixels$recovery_days), 1), "days\n")
  cat("  Median recovery time:", round(median(recovered_pixels$recovery_days), 1), "days\n")
  cat("  SD recovery time:", round(sd(recovered_pixels$recovery_days), 1), "days\n")
  cat("  Min recovery time:", round(min(recovered_pixels$recovery_days), 1), "days\n")
  cat("  Max recovery time:", round(max(recovered_pixels$recovery_days), 1), "days\n")
  
  # Recovery method breakdown
  method_table <- table(recovered_pixels$recovery_method)
  cat("\nRecovery Methods Used:\n")
  print(method_table)
}

log_message(sprintf("Recovery calculated for %d pixels (%d recovered)", 
                    nrow(recovery_results), nrow(recovered_pixels)), 
            log_file, "INFO")

# ============================================================================
# STEP 6: COMBINE RESULTS
# ============================================================================

cat("\n=== STEP 6: COMBINING RESULTS ===\n")

# Combine all results
final_results <- cbind(
  x = original_coords[, 1],
  y = original_coords[, 2],
  cluster = example_cluster,
  fire_year = example_fire_year,
  impact_results,
  recovery_results[, -1]  # Remove duplicate pixel_id
)

# Save results
results_file <- file.path(output_dirs$tables, 
                          sprintf("results_cluster_%d_fire_%d.csv", 
                                  example_cluster, example_fire_year))
write.csv(final_results, results_file, row.names = FALSE)
cat("Results saved to:", results_file, "\n")

log_message(sprintf("Results saved for cluster %d fire year %d", 
                    example_cluster, example_fire_year), 
            log_file, "INFO")

# ============================================================================
# STEP 7: CREATE VISUALIZATIONS
# ============================================================================

cat("\n=== STEP 7: CREATING VISUALIZATIONS ===\n")

# Extract dates from column names
dates <- extract_dates_from_names(names(extracted_values))

# Create sample diagnostic plots
cat("\nCreating time series diagnostic plots...\n")
sample_plots <- save_sample_diagnostic_plots(
  burnt_matrix = as.matrix(extracted_values),
  unburnt_ts = unburnt_ts,
  dates = dates,
  impact_results = impact_results,
  recovery_results = recovery_results,
  output_dir = output_dirs$timeseries,
  cluster_id = example_cluster,
  fire_year = example_fire_year,
  sample_size = 3
)

# Create recovery histogram
cat("Creating recovery histogram...\n")
hist_plot <- create_recovery_histogram(
  recovery_data = final_results,
  title = sprintf("Recovery Time Distribution - Cluster %d, Fire Year %d", 
                  example_cluster, example_fire_year),
  binwidth = 100
)

if (!is.null(hist_plot)) {
  ggsave(file.path(output_dirs$summary, 
                   sprintf("recovery_histogram_cluster_%d_fire_%d.png", 
                           example_cluster, example_fire_year)),
         hist_plot, width = 10, height = 6, dpi = 300)
  cat("Saved recovery histogram\n")
}

# Create impact vs recovery scatter plot
cat("Creating impact vs recovery scatter plot...\n")
scatter_plot <- create_impact_recovery_scatter(
  combined_data = final_results,
  add_smoothing = TRUE
)

if (!is.null(scatter_plot)) {
  ggsave(file.path(output_dirs$summary, 
                   sprintf("impact_recovery_scatter_cluster_%d_fire_%d.png", 
                           example_cluster, example_fire_year)),
         scatter_plot, width = 10, height = 8, dpi = 300)
  cat("Saved impact vs recovery scatter plot\n")
}

# Create recovery method comparison
cat("Creating recovery method comparison plot...\n")
method_plot <- create_recovery_method_plot(
  recovery_data = final_results
)

if (!is.null(method_plot)) {
  ggsave(file.path(output_dirs$summary, 
                   sprintf("recovery_methods_cluster_%d_fire_%d.png", 
                           example_cluster, example_fire_year)),
         method_plot, width = 10, height = 6, dpi = 300)
  cat("Saved recovery method comparison plot\n")
}

log_message("Visualizations created successfully", log_file, "INFO")

# ============================================================================
# STEP 8: GENERATE SUMMARY REPORT
# ============================================================================

cat("\n=== STEP 8: GENERATING SUMMARY REPORT ===\n")

# Calculate comprehensive statistics
report_stats <- final_results %>%
  summarise(
    total_pixels = n(),
    pixels_with_impact = sum(!is.na(impact_magnitude)),
    pixels_recovered = sum(recovered == TRUE, na.rm = TRUE),
    recovery_rate_pct = pixels_recovered / total_pixels * 100,
    mean_impact = mean(impact_magnitude, na.rm = TRUE),
    sd_impact = sd(impact_magnitude, na.rm = TRUE),
    median_impact = median(impact_magnitude, na.rm = TRUE),
    mean_recovery_days = mean(recovery_days[recovery_days > 0], na.rm = TRUE),
    sd_recovery_days = sd(recovery_days[recovery_days > 0], na.rm = TRUE),
    median_recovery_days = median(recovery_days[recovery_days > 0], na.rm = TRUE),
    min_recovery_days = min(recovery_days[recovery_days > 0], na.rm = TRUE),
    max_recovery_days = max(recovery_days[recovery_days > 0], na.rm = TRUE),
    q25_recovery = quantile(recovery_days[recovery_days > 0], 0.25, na.rm = TRUE),
    q75_recovery = quantile(recovery_days[recovery_days > 0], 0.75, na.rm = TRUE)
  )

# Create report text
report_content <- sprintf("
================================================================================
FIRE RECOVERY ANALYSIS REPORT
================================================================================

Analysis Date: %s
Cluster ID: %d
Fire Year: %d

--------------------------------------------------------------------------------
DATA SUMMARY
--------------------------------------------------------------------------------
Total pixels analyzed: %d
Pixels with detected impact: %d
Pixels recovered: %d
Recovery rate: %.1f%%

--------------------------------------------------------------------------------
FIRE IMPACT STATISTICS
--------------------------------------------------------------------------------
Mean impact magnitude: %.4f (±%.4f)
Median impact magnitude: %.4f
Range: %.4f - %.4f

--------------------------------------------------------------------------------
RECOVERY TIME STATISTICS (for recovered pixels only)
--------------------------------------------------------------------------------
Mean recovery time: %.1f days (±%.1f)
Median recovery time: %.1f days
Range: %.1f - %.1f days
Interquartile range: %.1f - %.1f days

--------------------------------------------------------------------------------
RECOVERY METHOD BREAKDOWN
--------------------------------------------------------------------------------
%s

--------------------------------------------------------------------------------
FILES GENERATED
--------------------------------------------------------------------------------
Results CSV: %s
Diagnostic plots: %s
Summary plots: %s

--------------------------------------------------------------------------------
REFERENCE INFORMATION
--------------------------------------------------------------------------------
Modal peak position: %d (observation within year)
Reference mean FPAR: %.4f
Number of reference values: %d

================================================================================
END OF REPORT
================================================================================
",
                          Sys.time(),
                          example_cluster,
                          example_fire_year,
                          report_stats$total_pixels,
                          report_stats$pixels_with_impact,
                          report_stats$pixels_recovered,
                          report_stats$recovery_rate_pct,
                          report_stats$mean_impact,
                          report_stats$sd_impact,
                          report_stats$median_impact,
                          min(final_results$impact_magnitude, na.rm = TRUE),
                          max(final_results$impact_magnitude, na.rm = TRUE),
                          report_stats$mean_recovery_days,
                          report_stats$sd_recovery_days,
                          report_stats$median_recovery_days,
                          report_stats$min_recovery_days,
                          report_stats$max_recovery_days,
                          report_stats$q25_recovery,
                          report_stats$q75_recovery,
                          paste(capture.output(print(table(final_results$recovery_method))), collapse = "\n"),
                          results_file,
                          output_dirs$timeseries,
                          output_dirs$summary,
                          ifelse(nrow(recovered_pixels) > 0, 
                                 recovered_pixels$modal_peak_position[1], NA),
                          ifelse(nrow(recovered_pixels) > 0, 
                                 recovered_pixels$reference_mean[1], NA),
                          sum(!is.na(unburnt_ts))
)

# Save report
report_file <- file.path(output_dirs$diagnostics, 
                         sprintf("analysis_report_cluster_%d_fire_%d.txt", 
                                 example_cluster, example_fire_year))
writeLines(report_content, report_file)
cat("Analysis report saved to:", report_file, "\n")

# Print report to console
cat(report_content)

log_message("Analysis completed successfully", log_file, "INFO")

# ============================================================================
# STEP 9: CLEANUP AND FINAL MESSAGE
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All outputs saved to:", output_dirs$base, "\n")
cat("Check the following directories:\n")
cat("  - Tables:", output_dirs$tables, "\n")
cat("  - Time series plots:", output_dirs$timeseries, "\n")
cat("  - Summary plots:", output_dirs$summary, "\n")
cat("  - Diagnostics:", output_dirs$diagnostics, "\n")
cat("\nLog file:", log_file, "\n")
cat("\n================================================================================\n")