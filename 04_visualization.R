#' Visualization Functions for Fire Recovery Analysis
#' 
#' Creates diagnostic plots and summary visualizations
#' 
#' @author Your Name
#' @date 2025

library(ggplot2)
library(gridExtra)
library(viridis)
library(dplyr)

source("scripts/utils/helper_functions.R")
source("scripts/utils/reference_value_calculation.R")

#' Create Time Series Diagnostic Plot
#'
#' Creates a diagnostic plot showing burnt vs unburnt time series with recovery markers
#'
#' @param burnt_ts Numeric vector of burnt pixel FPAR time series
#' @param unburnt_ts Numeric vector of unburnt reference time series
#' @param dates Vector of Date objects corresponding to time series
#' @param fire_impact_result Output from detect_fire_impact()
#' @param recovery_result Output from calculate_recovery_time()
#' @param pixel_coords Optional vector c(x, y) with pixel coordinates
#' @param cluster_id Cluster identifier
#' @param fire_year Fire year
#' @return ggplot object
#' @export
create_timeseries_plot <- function(burnt_ts, unburnt_ts, dates, 
                                   fire_impact_result, recovery_result,
                                   pixel_coords = NULL, cluster_id = NA, fire_year = NA) {
  
  # Create data frame
  ts_data <- data.frame(
    date = dates,
    burnt = burnt_ts,
    unburnt = unburnt_ts,
    time_idx = 1:length(dates)
  )
  
  # Remove missing data
  ts_data <- ts_data %>%
    filter(!is.na(date), !is.na(burnt), !is.na(unburnt))
  
  if (nrow(ts_data) < 10) {
    warning("Insufficient data points for plotting")
    return(NULL)
  }
  
  # Base plot
  p <- ggplot(ts_data, aes(x = date)) +
    geom_line(aes(y = unburnt, color = "Unburnt"), linewidth = 0.8, alpha = 0.8) +
    geom_line(aes(y = burnt, color = "Burnt"), linewidth = 0.8) +
    scale_color_manual(values = c("Unburnt" = "blue", "Burnt" = "orange"))
  
  # Add reference peak points if available
  if (!is.null(recovery_result$reference_values) && 
      length(recovery_result$future_peak_indices) > 0) {
    
    peak_indices <- recovery_result$future_peak_indices
    peak_indices <- peak_indices[peak_indices <= nrow(ts_data)]
    
    if (length(peak_indices) > 0) {
      ref_df <- data.frame(
        date = ts_data$date[peak_indices],
        value = unburnt_ts[peak_indices]
      )
      
      p <- p + geom_point(
        data = ref_df,
        aes(x = date, y = value),
        color = "purple", size = 2, alpha = 0.6, shape = 15
      )
    }
  }
  
  # Mark fire impact point
  if (!is.na(fire_impact_result$impact_date_idx) && 
      fire_impact_result$impact_date_idx <= nrow(ts_data) &&
      !is.na(fire_impact_result$impact_value)) {
    
    impact_date <- ts_data$date[fire_impact_result$impact_date_idx]
    impact_value <- fire_impact_result$impact_value
    
    p <- p +
      geom_vline(xintercept = impact_date, color = "red", 
                 linetype = "dashed", alpha = 0.7) +
      annotate("point", x = impact_date, y = impact_value,
               color = "red", size = 3, shape = 17)
  }
  
  # Mark recovery point
  if (!is.na(recovery_result$recovery_idx) && 
      recovery_result$recovery_idx > 0 && 
      recovery_result$recovery_idx <= nrow(ts_data)) {
    
    recovery_date <- ts_data$date[recovery_result$recovery_idx]
    recovery_value <- ts_data$burnt[recovery_result$recovery_idx]
    
    if (!is.na(recovery_value)) {
      p <- p +
        geom_vline(xintercept = recovery_date, color = "green", 
                   linetype = "dashed", alpha = 0.7) +
        annotate("point", x = recovery_date, y = recovery_value,
                 color = "green", size = 3, shape = 16)
    }
  }
  
  # Add recovery threshold line
  if (!is.null(recovery_result$reference_mean) && 
      !is.na(recovery_result$reference_mean)) {
    p <- p + geom_hline(
      yintercept = recovery_result$reference_mean * 0.95,
      color = "darkgreen", linetype = "dotted", alpha = 0.7
    )
  }
  
  # Create subtitle
  coord_text <- if (!is.null(pixel_coords)) {
    sprintf("Coords: (%.3f, %.3f)", pixel_coords[1], pixel_coords[2])
  } else {
    ""
  }
  
  subtitle_text <- sprintf(
    "Impact: %.3f, Recovery: %s days, Method: %s, Ref Mean: %.3f\n%s",
    ifelse(is.na(fire_impact_result$impact_magnitude), 0, 
           fire_impact_result$impact_magnitude),
    ifelse(recovery_result$recovery_days < 0, "No recovery", 
           as.character(round(recovery_result$recovery_days))),
    recovery_result$recovery_method,
    ifelse(is.null(recovery_result$reference_mean), NA, 
           recovery_result$reference_mean),
    coord_text
  )
  
  # Final formatting
  p <- p +
    labs(
      title = sprintf("Fire Recovery - Cluster %s, Fire Year %s", 
                      cluster_id, fire_year),
      subtitle = subtitle_text,
      x = "Date", 
      y = "FPAR", 
      color = "Series Type",
      caption = "Purple squares: Peak periods | Red: Fire impact | Green: Recovery | Dotted: 95% threshold"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 7, hjust = 0)
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  
  return(p)
}

#' Create Summary Histogram of Recovery Times
#'
#' Creates histogram showing distribution of recovery times
#'
#' @param recovery_data Data frame with recovery_days column
#' @param title Plot title
#' @param binwidth Histogram bin width in days (default = 100)
#' @return ggplot object
#' @export
create_recovery_histogram <- function(recovery_data, title = "Recovery Time Distribution", 
                                      binwidth = 100) {
  
  # Filter valid recovery times
  valid_recovery <- recovery_data %>%
    filter(recovery_days > 0, recovery_days < 9999)
  
  if (nrow(valid_recovery) == 0) {
    warning("No valid recovery data for histogram")
    return(NULL)
  }
  
  p <- ggplot(valid_recovery, aes(x = recovery_days)) +
    geom_histogram(binwidth = binwidth, fill = "skyblue", 
                   color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = median(recovery_days)), 
               color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", 
             x = median(valid_recovery$recovery_days) + 200,
             y = Inf,
             label = sprintf("Median: %.0f days", median(valid_recovery$recovery_days)),
             vjust = 2, color = "red", fontface = "bold") +
    labs(
      title = title,
      subtitle = sprintf("n = %d pixels, Mean: %.0f days, SD: %.0f days",
                         nrow(valid_recovery),
                         mean(valid_recovery$recovery_days),
                         sd(valid_recovery$recovery_days)),
      x = "Recovery Time (days)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}

#' Create Impact vs Recovery Scatter Plot
#'
#' Plots relationship between fire impact magnitude and recovery time
#'
#' @param combined_data Data frame with impact_magnitude and recovery_days columns
#' @param add_smoothing Add loess smoothing line (default = TRUE)
#' @return ggplot object
#' @export
create_impact_recovery_scatter <- function(combined_data, add_smoothing = TRUE) {
  
  # Filter valid data
  valid_data <- combined_data %>%
    filter(!is.na(impact_magnitude), 
           recovery_days > 0, 
           recovery_days < 9999)
  
  if (nrow(valid_data) == 0) {
    warning("No valid data for scatter plot")
    return(NULL)
  }
  
  p <- ggplot(valid_data, aes(x = impact_magnitude, y = recovery_days)) +
    geom_point(alpha = 0.4, color = "steelblue") +
    labs(
      title = "Fire Impact Magnitude vs Recovery Time",
      subtitle = sprintf("n = %d pixels, Correlation: %.3f",
                         nrow(valid_data),
                         cor(valid_data$impact_magnitude, 
                             valid_data$recovery_days, 
                             method = "spearman")),
      x = "Impact Magnitude (FPAR decline)",
      y = "Recovery Time (days)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  if (add_smoothing && nrow(valid_data) > 10) {
    p <- p + geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 1)
  }
  
  return(p)
}

#' Create Recovery Rate by Cluster Bar Plot
#'
#' Shows recovery success rate for each cluster
#'
#' @param summary_data Data frame with cluster and recovery rate columns
#' @return ggplot object
#' @export
create_recovery_rate_barplot <- function(summary_data) {
  
  # Calculate recovery rates
  cluster_summary <- summary_data %>%
    group_by(cluster) %>%
    summarise(
      total = n(),
      recovered = sum(recovered == TRUE, na.rm = TRUE),
      recovery_rate = recovered / total * 100,
      .groups = 'drop'
    ) %>%
    arrange(desc(recovery_rate))
  
  p <- ggplot(cluster_summary, aes(x = reorder(factor(cluster), recovery_rate), 
                                   y = recovery_rate)) +
    geom_bar(stat = "identity", fill = "forestgreen", alpha = 0.7) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red", alpha = 0.5) +
    coord_flip() +
    labs(
      title = "Recovery Success Rate by Cluster",
      subtitle = sprintf("Overall recovery rate: %.1f%%", 
                         mean(cluster_summary$recovery_rate)),
      x = "Cluster",
      y = "Recovery Rate (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  return(p)
}

#' Create Recovery Method Comparison Plot
#'
#' Compares recovery times by detection method
#'
#' @param recovery_data Data frame with recovery_method and recovery_days columns
#' @return ggplot object
#' @export
create_recovery_method_plot <- function(recovery_data) {
  
  # Filter valid data
  valid_data <- recovery_data %>%
    filter(recovery_days > 0, recovery_days < 9999,
           recovery_method %in% c("stage1_dynamic_catchup", "stage2_static_threshold"))
  
  if (nrow(valid_data) == 0) {
    warning("No valid data for method comparison")
    return(NULL)
  }
  
  # Rename methods for clarity
  valid_data <- valid_data %>%
    mutate(
      method_label = case_when(
        recovery_method == "stage1_dynamic_catchup" ~ "Stage 1: Dynamic Catchup",
        recovery_method == "stage2_static_threshold" ~ "Stage 2: Static Threshold",
        TRUE ~ recovery_method
      )
    )
  
  p <- ggplot(valid_data, aes(x = method_label, y = recovery_days, fill = method_label)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.5) +
    scale_fill_manual(values = c("Stage 1: Dynamic Catchup" = "#2ecc71", 
                                 "Stage 2: Static Threshold" = "#f39c12")) +
    labs(
      title = "Recovery Time by Detection Method",
      subtitle = sprintf("Stage 1: n=%d, Stage 2: n=%d",
                         sum(valid_data$recovery_method == "stage1_dynamic_catchup"),
                         sum(valid_data$recovery_method == "stage2_static_threshold")),
      x = "Detection Method",
      y = "Recovery Time (days)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  return(p)
}

#' Save Diagnostic Plots for Sample Pixels
#'
#' Creates and saves time series plots for sample pixels
#'
#' @param burnt_matrix Matrix of burnt pixel time series
#' @param unburnt_ts Unburnt reference time series
#' @param dates Vector of dates
#' @param impact_results Impact detection results
#' @param recovery_results Recovery estimation results
#' @param output_dir Output directory for plots
#' @param cluster_id Cluster identifier
#' @param fire_year Fire year
#' @param sample_size Number of pixels to plot (default = 3)
#' @export
save_sample_diagnostic_plots <- function(burnt_matrix, unburnt_ts, dates,
                                         impact_results, recovery_results,
                                         output_dir, cluster_id, fire_year,
                                         sample_size = 3) {
  
  # Create output directory if needed
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Sample pixels
  n_pixels <- min(nrow(burnt_matrix), sample_size)
  sample_indices <- sample(1:nrow(burnt_matrix), n_pixels)
  
  plots_list <- list()
  
  for (i in seq_along(sample_indices)) {
    idx <- sample_indices[i]
    
    # Get data for this pixel
    burnt_ts <- as.numeric(burnt_matrix[idx, ])
    
    impact_result <- list(
      impact_magnitude = impact_results$impact_magnitude[idx],
      impact_date_idx = impact_results$impact_date_idx[idx],
      pre_fire_mean = impact_results$pre_fire_mean[idx],
      impact_value = impact_results$impact_value[idx],
      unburnt_reference = impact_results$unburnt_reference[idx]
    )
    
    # Get recovery result for this pixel
    recovery_result <- list(
      recovery_days = recovery_results$recovery_days[idx],
      recovery_idx = recovery_results$recovery_idx[idx],
      recovered = recovery_results$recovered[idx],
      reference_mean = recovery_results$reference_mean[idx],
      modal_peak_position = recovery_results$modal_peak_position[idx],
      recovery_method = recovery_results$recovery_method[idx],
      reference_values = NA,  # Would need to pass this separately
      future_peak_indices = NA
    )
    
    # Create plot
    p <- create_timeseries_plot(
      burnt_ts = burnt_ts,
      unburnt_ts = unburnt_ts,
      dates = dates,
      fire_impact_result = impact_result,
      recovery_result = recovery_result,
      cluster_id = cluster_id,
      fire_year = fire_year
    )
    
    if (!is.null(p)) {
      plots_list[[i]] <- p
    }
  }
  
  # Save combined plot
  if (length(plots_list) > 0) {
    combined_plot <- gridExtra::arrangeGrob(grobs = plots_list, ncol = 1)
    
    output_file <- file.path(output_dir, 
                             sprintf("timeseries_cluster_%s_fire_%s.png", 
                                     cluster_id, fire_year))
    
    ggsave(output_file, combined_plot, 
           width = 14, height = 4 * length(plots_list), dpi = 300)
    
    cat(sprintf("Saved diagnostic plots: %s\n", output_file))
  }
  
  return(plots_list)
}