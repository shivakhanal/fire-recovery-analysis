#' Cluster Homogeneity Analysis
#' 
#' Assesses clustering quality using burnt pixel time series

library(data.table)
library(terra)
library(sf)
library(ggplot2)

#' Calculate Cluster Homogeneity Metrics
#'
#' Calculates homogeneity and consistency scores for a cluster
#'
#' @param homogeneity_results Data table with individual pixel correlations
#' @param cluster_id Cluster identifier
#' @return Data table with cluster summary statistics
#' @export
calculate_cluster_summary_stats <- function(homogeneity_results, cluster_id) {
  
  # Filter valid data
  valid_data <- homogeneity_results[consistency_category != "Insufficient_Data"]
  
  if (nrow(valid_data) == 0) {
    return(data.table(
      cluster = cluster_id,
      total_valid_pixels = 0,
      mean_correlation = 0,
      homogeneity_score = 0,
      consistency_score = 0,
      overall_quality = 0
    ))
  }
  
  # Calculate summary metrics
  summary_stats <- data.table(
    cluster = cluster_id,
    total_valid_pixels = nrow(valid_data),
    
    # Consistency counts
    highly_consistent = sum(valid_data$consistency_category == "Highly_Consistent"),
    moderately_consistent = sum(valid_data$consistency_category == "Moderately_Consistent"),
    weakly_consistent = sum(valid_data$consistency_category == "Weakly_Consistent"),
    consistently_above = sum(valid_data$consistency_category == "Consistently_Above"),
    consistently_below = sum(valid_data$consistency_category == "Consistently_Below"),
    inconsistent = sum(valid_data$consistency_category == "Inconsistent"),
    
    # Summary metrics
    mean_correlation = mean(valid_data$correlation_with_mean, na.rm = TRUE),
    median_correlation = median(valid_data$correlation_with_mean, na.rm = TRUE),
    sd_correlation = sd(valid_data$correlation_with_mean, na.rm = TRUE),
    mean_deviation = mean(valid_data$mean_absolute_deviation, na.rm = TRUE),
    median_deviation = median(valid_data$mean_absolute_deviation, na.rm = TRUE)
  )
  
  # Calculate quality scores (0-1, higher = better)
  total_pixels <- summary_stats$total_valid_pixels
  
  summary_stats[, `:=`(
    # Homogeneity Score: Mean correlation (higher = better clustering)
    homogeneity_score = mean_correlation,
    
    # Consistency Score: Proportion with at least moderate correlation (r > 0.5)
    consistency_score = (highly_consistent + moderately_consistent) / total_pixels,
    
    # Proportion breakdowns
    prop_highly_consistent = highly_consistent / total_pixels,
    prop_moderately_consistent = moderately_consistent / total_pixels,
    prop_weakly_consistent = weakly_consistent / total_pixels,
    prop_inconsistent = inconsistent / total_pixels
  )]
  
  # Overall Quality Score: Mean of Homogeneity and Consistency
  summary_stats[, overall_quality := (homogeneity_score + consistency_score) / 2]
  
  return(summary_stats)
}

#' Analyze Cluster Homogeneity
#'
#' Analyzes homogeneity for a single cluster using burnt pixels
#'
#' @param cluster_id Cluster identifier
#' @param time_stack SpatRaster stack with FPAR time series
#' @param cluster_points Data frame with cluster assignments and coordinates
#' @param fire_polygons SF object with fire polygons
#' @param dates Vector of dates corresponding to time series
#' @param min_burnt_pixels Minimum burnt pixels required (default = 10)
#' @return List with homogeneity results and cluster summary
#' @export
analyze_cluster_homogeneity_revised <- function(cluster_id,
                                                time_stack,
                                                cluster_points,
                                                fire_polygons,
                                                dates = NULL,
                                                min_burnt_pixels = 10) {
  
  cat(sprintf("\n=== Analyzing Cluster %d ===\n", cluster_id))
  
  # Get cluster points
  cluster_pts <- cluster_points[cluster_points$class == cluster_id, ]
  
  if (nrow(cluster_pts) == 0) {
    cat("No points in cluster\n")
    return(NULL)
  }
  
  cat(sprintf("Total cluster pixels: %d\n", nrow(cluster_pts)))
  
  # Convert to sf
  cluster_sf <- st_as_sf(cluster_pts, coords = c("x", "y"), crs = crs(time_stack))
  
  # Ensure CRS match with fire polygons
  if (st_crs(cluster_sf) != st_crs(fire_polygons)) {
    cluster_sf <- st_transform(cluster_sf, st_crs(fire_polygons))
  }
  
  # Find burnt points
  burnt_intersection <- st_intersects(cluster_sf, fire_polygons)
  burnt_idx <- which(lengths(burnt_intersection) > 0)
  
  if (length(burnt_idx) < min_burnt_pixels) {
    cat(sprintf("Insufficient burnt pixels: %d (minimum: %d)\n", 
                length(burnt_idx), min_burnt_pixels))
    return(NULL)
  }
  
  cat(sprintf("Burnt pixels: %d\n", length(burnt_idx)))
  
  # Get burnt points
  burnt_sf <- cluster_sf[burnt_idx, ]
  original_coords <- st_coordinates(burnt_sf)
  
  # Transform to raster CRS for extraction
  burnt_sf <- st_transform(burnt_sf, crs(time_stack))
  burnt_vect <- vect(burnt_sf)
  
  # Extract time series
  extracted <- terra::extract(burnt_vect, time_stack)
  
  if ("ID" %in% names(extracted)) {
    extracted <- extracted[, -which(names(extracted) == "ID")]
  }
  
  # Convert to matrix
  ts_matrix <- as.matrix(extracted)
  
  # Calculate cluster mean time series
  cluster_mean_ts <- colMeans(ts_matrix, na.rm = TRUE)
  
  # Calculate correlation for each pixel
  correlations <- apply(ts_matrix, 1, function(pixel_ts) {
    if (sum(!is.na(pixel_ts)) < 10) return(NA)  # Need at least 10 valid observations
    cor(pixel_ts, cluster_mean_ts, use = "pairwise.complete.obs")
  })
  
  # Calculate mean absolute deviation
  deviations <- apply(ts_matrix, 1, function(pixel_ts) {
    mean(abs(pixel_ts - cluster_mean_ts), na.rm = TRUE)
  })
  
  # Classify pixels by consistency
  consistency_category <- sapply(correlations, function(r) {
    if (is.na(r)) return("Insufficient_Data")
    if (r > 0.8) return("Highly_Consistent")
    if (r > 0.5) return("Moderately_Consistent")
    if (r > 0.3) return("Weakly_Consistent")
    if (r < 0) return("Inconsistent")
    return("Weakly_Consistent")
  })
  
  # Create results data table
  homogeneity_results <- data.table(
    cluster = cluster_id,
    pixel_id = 1:length(burnt_idx),
    x = original_coords[, 1],
    y = original_coords[, 2],
    correlation_with_mean = correlations,
    mean_absolute_deviation = deviations,
    consistency_category = consistency_category
  )
  
  # Calculate cluster summary
  cluster_summary <- calculate_cluster_summary_stats(homogeneity_results, cluster_id)
  
  # Sample data for plotting
  valid_pixels <- which(!is.na(correlations))
  if (length(valid_pixels) > 0) {
    sample_size <- min(5, length(valid_pixels))
    sample_idx <- sample(valid_pixels, sample_size)
    
    sample_data <- list(
      cluster_id = cluster_id,
      sample_indices = sample_idx,
      sample_ts = ts_matrix[sample_idx, , drop = FALSE],
      cluster_mean_ts = cluster_mean_ts,
      sample_correlations = correlations[sample_idx],
      dates = dates
    )
  } else {
    sample_data <- NULL
  }
  
  cat(sprintf("Homogeneity Score: %.3f\n", cluster_summary$homogeneity_score))
  cat(sprintf("Consistency Score: %.3f\n", cluster_summary$consistency_score))
  cat(sprintf("Overall Quality: %.3f\n", cluster_summary$overall_quality))
  
  return(list(
    homogeneity_results = homogeneity_results,
    cluster_summary = cluster_summary,
    sample_data = sample_data
  ))
}

#' Run Complete Homogeneity Analysis
#'
#' Analyzes homogeneity for multiple clusters
#'
#' @param clusters_to_analyze Vector of cluster IDs to analyze
#' @param time_stack SpatRaster with FPAR time series
#' @param cluster_points Data frame with cluster assignments
#' @param fire_polygons SF object with fire polygons
#' @param dates Vector of dates (optional)
#' @param output_path Output directory path
#' @param cache_path Cache directory path
#' @return List with combined results
#' @export
run_revised_homogeneity_analysis <- function(clusters_to_analyze = 1:50,
                                             time_stack = stk,
                                             cluster_points = clus,
                                             fire_polygons = fire,
                                             dates = NULL,
                                             output_path = "outputs/cluster_homogeneity",
                                             cache_path = "outputs/cluster_homogeneity/cache") {
  
  cat("=== REVISED CLUSTER HOMOGENEITY ANALYSIS ===\n")
  cat("Analyzing clusters:", paste(range(clusters_to_analyze), collapse = "-"), "\n")
  cat("Total clusters to process:", length(clusters_to_analyze), "\n\n")
  
  # Create directories
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
  
  # Storage
  all_homogeneity_results <- list()
  all_cluster_summaries <- list()
  successful_clusters <- c()
  
  # Process each cluster
  for (cluster_id in clusters_to_analyze) {
    tryCatch({
      result <- analyze_cluster_homogeneity_revised(
        cluster_id = cluster_id,
        time_stack = time_stack,
        cluster_points = cluster_points,
        fire_polygons = fire_polygons,
        dates = dates
      )
      
      if (!is.null(result)) {
        # Store results
        all_homogeneity_results[[as.character(cluster_id)]] <- result$homogeneity_results
        all_cluster_summaries[[as.character(cluster_id)]] <- result$cluster_summary
        
        successful_clusters <- c(successful_clusters, cluster_id)
        
        # Save individual results
        fwrite(result$homogeneity_results, 
               paste0(cache_path, "/cluster_", cluster_id, "_homogeneity.csv"))
        fwrite(result$cluster_summary,
               paste0(cache_path, "/cluster_", cluster_id, "_summary.csv"))
        
        # Create time series plot if sample data available
        if (!is.null(result$sample_data)) {
          tryCatch({
            create_homogeneity_plot(result$sample_data, output_path)
          }, error = function(e) {
            cat("Warning: Could not create plot for cluster", cluster_id, "\n")
          })
        }
      }
      
      # Memory cleanup
      rm(result)
      if (cluster_id %% 5 == 0) {
        gc(verbose = FALSE)
        cat("Memory cleanup after cluster", cluster_id, "\n")
      }
      
    }, error = function(e) {
      cat("Error in cluster", cluster_id, ":", e$message, "\n")
    })
  }
  
  if (length(successful_clusters) == 0) {
    cat("No clusters successfully analyzed\n")
    return(NULL)
  }
  
  # Combine results
  combined_homogeneity <- rbindlist(all_homogeneity_results, fill = TRUE)
  combined_summaries <- rbindlist(all_cluster_summaries, fill = TRUE)
  
  # Save combined results
  fwrite(combined_homogeneity, paste0(output_path, "/all_homogeneity_results.csv"))
  fwrite(combined_summaries, paste0(output_path, "/cluster_homogeneity_summary.csv"))
  
  # Generate summary report
  generate_homogeneity_report(combined_summaries, output_path)
  
  cat("\n=== HOMOGENEITY ANALYSIS COMPLETE ===\n")
  cat("Successfully analyzed clusters:", length(successful_clusters), "\n")
  cat("Total valid pixels:", nrow(combined_homogeneity), "\n")
  
  # Show data variation
  cat("\nMetric ranges across clusters:\n")
  cat(sprintf("  Homogeneity Score: %.3f - %.3f\n", 
              min(combined_summaries$homogeneity_score, na.rm = TRUE),
              max(combined_summaries$homogeneity_score, na.rm = TRUE)))
  cat(sprintf("  Consistency Score: %.3f - %.3f\n",
              min(combined_summaries$consistency_score, na.rm = TRUE),
              max(combined_summaries$consistency_score, na.rm = TRUE)))
  cat(sprintf("  Overall Quality: %.3f - %.3f\n",
              min(combined_summaries$overall_quality, na.rm = TRUE),
              max(combined_summaries$overall_quality, na.rm = TRUE)))
  
  # Top performers
  cat("\nTop 10 clusters by overall quality:\n")
  top_clusters <- combined_summaries[order(-overall_quality)][1:min(10, nrow(combined_summaries))]
  print(top_clusters[, .(cluster, overall_quality, homogeneity_score, 
                         consistency_score, total_valid_pixels)])
  
  # Bottom performers
  cat("\nBottom 10 clusters by overall quality:\n")
  bottom_clusters <- combined_summaries[order(overall_quality)][1:min(10, nrow(combined_summaries))]
  print(bottom_clusters[, .(cluster, overall_quality, homogeneity_score, 
                            consistency_score, total_valid_pixels)])
  
  cat("\nResults saved to:", output_path, "\n")
  
  return(list(
    individual_results = combined_homogeneity,
    cluster_summaries = combined_summaries,
    successful_clusters = successful_clusters
  ))
}

#' Create Homogeneity Visualization
#'
#' Creates time series plot showing cluster homogeneity
#'
#' @param sample_data List with sample data from analyze_cluster_homogeneity_revised
#' @param output_path Output directory
#' @export
create_homogeneity_plot <- function(sample_data, output_path) {
  
  if (is.null(sample_data) || is.null(sample_data$dates)) {
    return(NULL)
  }
  
  cluster_id <- sample_data$cluster_id
  
  # Prepare data
  n_samples <- nrow(sample_data$sample_ts)
  plot_data <- data.frame()
  
  for (i in 1:n_samples) {
    pixel_data <- data.frame(
      date = sample_data$dates,
      value = sample_data$sample_ts[i, ],
      pixel = paste0("Pixel ", i, " (r=", 
                     round(sample_data$sample_correlations[i], 2), ")"),
      type = "Individual"
    )
    plot_data <- rbind(plot_data, pixel_data)
  }
  
  # Add cluster mean
  mean_data <- data.frame(
    date = sample_data$dates,
    value = sample_data$cluster_mean_ts,
    pixel = "Cluster Mean",
    type = "Mean"
  )
  plot_data <- rbind(plot_data, mean_data)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = date, y = value, color = pixel, linetype = type)) +
    geom_line(linewidth = ifelse(plot_data$type == "Mean", 1.2, 0.6),
              alpha = ifelse(plot_data$type == "Mean", 1, 0.7)) +
    scale_linetype_manual(values = c("Individual" = "solid", "Mean" = "dashed")) +
    labs(
      title = sprintf("Cluster %d: Sample Pixel Time Series", cluster_id),
      subtitle = "Individual burnt pixels vs. cluster mean",
      x = "Date",
      y = "FPAR",
      color = "Series",
      linetype = "Type"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save plot
  ggsave(paste0(output_path, "/cluster_", cluster_id, "_homogeneity.png"),
         p, width = 12, height = 6, dpi = 300)
  
  return(p)
}

#' Generate Homogeneity Analysis Report
#'
#' Creates text report summarizing homogeneity analysis
#'
#' @param cluster_summaries Data table with cluster summary statistics
#' @param output_path Output directory
#' @export
generate_homogeneity_report <- function(cluster_summaries, output_path) {
  
  report_text <- sprintf("
=== CLUSTER HOMOGENEITY ANALYSIS REPORT ===
Generated: %s

BASIC STATISTICS:
- Total clusters analyzed: %d
- Total pixels analyzed: %d
- Average pixels per cluster: %.1f

QUALITY SCORE STATISTICS:
- homogeneity_score:
  Mean: %.3f
  Range: %.3f - %.3f
  Std Dev: %.3f

- consistency_score:
  Mean: %.3f
  Range: %.3f - %.3f
  Std Dev: %.3f

- overall_quality:
  Mean: %.3f
  Range: %.3f - %.3f
  Std Dev: %.3f

TOP 10 PERFORMERS (by overall quality):
%s

BOTTOM 10 PERFORMERS (by overall quality):
%s

CLUSTER SIZE DISTRIBUTION:
- Small clusters (<1000 pixels): %d
- Medium clusters (1000-5000 pixels): %d
- Large clusters (>5000 pixels): %d
",
                         Sys.time(),
                         nrow(cluster_summaries),
                         sum(cluster_summaries$total_valid_pixels),
                         mean(cluster_summaries$total_valid_pixels),
                         mean(cluster_summaries$homogeneity_score, na.rm = TRUE),
                         min(cluster_summaries$homogeneity_score, na.rm = TRUE),
                         max(cluster_summaries$homogeneity_score, na.rm = TRUE),
                         sd(cluster_summaries$homogeneity_score, na.rm = TRUE),
                         mean(cluster_summaries$consistency_score, na.rm = TRUE),
                         min(cluster_summaries$consistency_score, na.rm = TRUE),
                         max(cluster_summaries$consistency_score, na.rm = TRUE),
                         sd(cluster_summaries$consistency_score, na.rm = TRUE),
                         mean(cluster_summaries$overall_quality, na.rm = TRUE),
                         min(cluster_summaries$overall_quality, na.rm = TRUE),
                         max(cluster_summaries$overall_quality, na.rm = TRUE),
                         sd(cluster_summaries$overall_quality, na.rm = TRUE),
                         paste(capture.output(
                           print(cluster_summaries[order(-overall_quality)][1:10, 
                                                                            .(cluster, overall_quality, total_valid_pixels)])
                         ), collapse = "\n"),
                         paste(capture.output(
                           print(cluster_summaries[order(overall_quality)][1:10, 
                                                                           .(cluster, overall_quality, total_valid_pixels)])
                         ), collapse = "\n"),
                         sum(cluster_summaries$total_valid_pixels < 1000),
                         sum(cluster_summaries$total_valid_pixels >= 1000 & 
                               cluster_summaries$total_valid_pixels <= 5000),
                         sum(cluster_summaries$total_valid_pixels > 5000)
  )
  
  # Write report
  writeLines(report_text, paste0(output_path, "/homogeneity_analysis_report.txt"))
  
  cat(report_text)
}