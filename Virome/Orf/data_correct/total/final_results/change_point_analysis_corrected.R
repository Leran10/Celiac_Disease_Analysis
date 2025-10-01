#!/usr/bin/env Rscript
# Change Point Analysis for Corrected Data
# Author: Claude AI
# Date: 2025-07-29

library(dplyr)
library(ggplot2)
library(changepoint)
library(gridExtra)
library(RColorBrewer)
library(stringr)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

cat("=== CHANGE POINT ANALYSIS FOR CORRECTED DATA ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading diversity data...\n")

# Load diversity data
diversity_data <- read.csv("diversity_data_full.csv")

cat("Loaded diversity data with", nrow(diversity_data), "samples\n")

# ============================================================================
# PART 2: PREPARE DATA FOR CHANGE POINT ANALYSIS
# ============================================================================

cat("\n2. Preparing data for change point analysis...\n")

# Calculate group means by timepoint for each metric
diversity_metrics <- c("richness", "shannon", "simpson", "evenness", "total_abundance", "dominance", "viral_load_cv")

group_means_data <- diversity_data %>%
  group_by(Dx.Status, onset_timeline_numeric) %>%
  summarise(
    across(all_of(diversity_metrics), ~mean(.x, na.rm = TRUE)),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(onset_timeline_numeric)

# Calculate differences between groups at each timepoint
timepoints <- sort(unique(group_means_data$onset_timeline_numeric))
change_point_data <- data.frame()

for(tp in timepoints) {
  tp_data <- group_means_data[group_means_data$onset_timeline_numeric == tp, ]
  
  if(nrow(tp_data) == 2) {  # Both CELIAC and CONTROL present
    celiac_data <- tp_data[tp_data$Dx.Status == "CELIAC", ]
    control_data <- tp_data[tp_data$Dx.Status == "CONTROL", ]
    
    row_data <- data.frame(timepoint = tp)
    
    for(metric in diversity_metrics) {
      # Calculate absolute difference between groups
      diff_val <- abs(control_data[[metric]] - celiac_data[[metric]])
      row_data[[paste0(metric, "_diff")]] <- diff_val
    }
    
    change_point_data <- rbind(change_point_data, row_data)
  }
}

cat("Prepared change point data for", nrow(change_point_data), "timepoints\n")

# ============================================================================
# PART 3: CHANGE POINT DETECTION
# ============================================================================

cat("\n3. Performing change point detection...\n")

change_point_results <- list()

for(metric in diversity_metrics) {
  cat("  Analyzing", metric, "...\n")
  
  diff_col <- paste0(metric, "_diff")
  
  if(diff_col %in% colnames(change_point_data)) {
    # Get the difference values
    diff_values <- change_point_data[[diff_col]]
    
    # Remove any NA or infinite values
    valid_indices <- which(!is.na(diff_values) & is.finite(diff_values))
    diff_values_clean <- diff_values[valid_indices]
    timepoints_clean <- change_point_data$timepoint[valid_indices]
    
    if(length(diff_values_clean) > 10) {
      
      # Method 1: PELT (Pruned Exact Linear Time)
      tryCatch({
        cpt_pelt <- cpt.mean(diff_values_clean, method="PELT", penalty="BIC")
        change_points_pelt <- cpts(cpt_pelt)
        if(length(change_points_pelt) > 0) {
          change_points_pelt_time <- timepoints_clean[change_points_pelt]
        } else {
          change_points_pelt_time <- numeric(0)
        }
      }, error = function(e) {
        change_points_pelt_time <- numeric(0)
      })
      
      # Method 2: Binary Segmentation
      tryCatch({
        cpt_binseg <- cpt.mean(diff_values_clean, method="BinSeg", Q=3)
        change_points_binseg <- cpts(cpt_binseg)
        if(length(change_points_binseg) > 0) {
          change_points_binseg_time <- timepoints_clean[change_points_binseg]
        } else {
          change_points_binseg_time <- numeric(0)
        }
      }, error = function(e) {
        change_points_binseg_time <- numeric(0)
      })
      
      # Method 3: Maximum difference approach
      tryCatch({
        # Simple approach: find maximum difference point
        max_diff_idx <- which.max(diff_values_clean)
        change_points_max_time <- timepoints_clean[max_diff_idx]
      }, error = function(e) {
        change_points_max_time <- numeric(0)
      })
      
      # Store results
      change_point_results[[metric]] <- list(
        timepoints = timepoints_clean,
        differences = diff_values_clean,
        pelt = change_points_pelt_time,
        binseg = change_points_binseg_time,
        max_method = change_points_max_time,
        max_diff_time = timepoints_clean[which.max(diff_values_clean)],
        max_diff_value = max(diff_values_clean)
      )
    }
  }
}

cat("Change point detection completed\n")

# ============================================================================
# PART 4: SUMMARIZE RESULTS
# ============================================================================

cat("\n4. Summarizing change point results...\n")

# Create summary table
change_point_summary <- data.frame()

for(metric in names(change_point_results)) {
  results <- change_point_results[[metric]]
  
  # Get consensus change point (most common or average)
  all_change_points <- c(results$pelt, results$binseg, results$max_method)
  
  if(length(all_change_points) > 0) {
    consensus_cp <- median(all_change_points, na.rm = TRUE)
  } else {
    consensus_cp <- results$max_diff_time
  }
  
  summary_row <- data.frame(
    metric = metric,
    consensus_change_point = consensus_cp,
    max_divergence_time = results$max_diff_time,
    max_divergence_value = results$max_diff_value,
    n_pelt_points = length(results$pelt),
    n_binseg_points = length(results$binseg),
    stringsAsFactors = FALSE
  )
  
  change_point_summary <- rbind(change_point_summary, summary_row)
}

# Calculate overall statistics
overall_mean_cp <- mean(change_point_summary$consensus_change_point, na.rm = TRUE)
overall_sd_cp <- sd(change_point_summary$consensus_change_point, na.rm = TRUE)
overall_mean_divergence <- mean(change_point_summary$max_divergence_time, na.rm = TRUE)

cat("Overall mean change point:", round(overall_mean_cp, 1), "months\n")
cat("Critical window (mean ± 1SD):", round(overall_mean_cp - overall_sd_cp, 1), "to", round(overall_mean_cp + overall_sd_cp, 1), "months\n")
cat("Mean maximum divergence time:", round(overall_mean_divergence, 1), "months\n")

# Save results
write.csv(change_point_summary, "change_point_summary_corrected.csv", row.names = FALSE)

# ============================================================================
# PART 5: CREATE VISUALIZATIONS
# ============================================================================

cat("\n5. Creating visualizations...\n")

# Function to create individual change point plots
create_change_point_plot <- function(metric, results, title_suffix = "") {
  
  plot_data <- data.frame(
    timepoint = results$timepoints,
    difference = results$differences
  )
  
  p <- ggplot(plot_data, aes(x = timepoint, y = difference)) +
    geom_line(size = 1, color = "#1f77b4", alpha = 0.8) +
    geom_point(size = 2, color = "#1f77b4", alpha = 0.8) +
    
    # Add change points
    {if(length(results$pelt) > 0) 
      geom_vline(xintercept = results$pelt, color = "#d62728", linetype = "dashed", alpha = 0.7)} +
    {if(length(results$binseg) > 0) 
      geom_vline(xintercept = results$binseg, color = "#ff7f0e", linetype = "dotted", alpha = 0.7)} +
    
    # Highlight maximum divergence
    geom_vline(xintercept = results$max_diff_time, color = "#2ca02c", size = 1.2, alpha = 0.8) +
    
    labs(
      title = paste("Change Point Analysis:", stringr::str_to_title(metric), title_suffix),
      subtitle = paste("Max divergence at", results$max_diff_time, "months"),
      x = "Months Relative to Onset",
      y = "Absolute Group Difference"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10)
    )
  
  return(p)
}

# Create individual plots
individual_plots <- list()
for(metric in names(change_point_results)) {
  individual_plots[[metric]] <- create_change_point_plot(metric, change_point_results[[metric]])
  
  # Save individual plot
  ggsave(
    filename = paste0("change_point_", metric, "_corrected.png"),
    plot = individual_plots[[metric]],
    width = 8, height = 6, dpi = 300
  )
}

# Create combined plot
if(length(individual_plots) >= 4) {
  combined_plot <- grid.arrange(
    individual_plots[[1]], individual_plots[[2]],
    individual_plots[[3]], individual_plots[[4]],
    ncol = 2, nrow = 2
  )
  
  ggsave("change_point_analysis_combined_corrected.png", combined_plot, 
         width = 14, height = 10, dpi = 300)
}

# Create divergence summary plot
create_divergence_summary <- function(change_point_summary) {
  
  # Prepare data for plotting
  plot_data <- change_point_summary %>%
    select(metric, consensus_change_point, max_divergence_time, max_divergence_value) %>%
    mutate(
      metric_clean = stringr::str_to_title(gsub("_", " ", metric))
    )
  
  p <- ggplot(plot_data, aes(x = reorder(metric_clean, -max_divergence_value))) +
    geom_col(aes(y = max_divergence_value), fill = "#1f77b4", alpha = 0.7) +
    geom_text(aes(y = max_divergence_value, 
                  label = paste("CP:", round(consensus_change_point, 1), "m")),
              vjust = -0.5, size = 3) +
    labs(
      title = "Change Point Summary: Maximum Group Divergence",
      subtitle = "With consensus change point times (CP) for corrected data",
      x = "Diversity Metric",
      y = "Maximum Divergence Value"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(p)
}

divergence_plot <- create_divergence_summary(change_point_summary)
ggsave("group_divergence_analysis_corrected.png", divergence_plot, 
       width = 10, height = 6, dpi = 300)

cat("All visualizations created\n")

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== CHANGE POINT ANALYSIS SUMMARY (CORRECTED DATA) ===\n")

cat("\nCHANGE POINT RESULTS:\n")
for(i in 1:nrow(change_point_summary)) {
  row <- change_point_summary[i, ]
  cat("- ", row$metric, ": CP = ", round(row$consensus_change_point, 1), 
      " months, Max divergence = ", round(row$max_divergence_value, 3), 
      " at ", round(row$max_divergence_time, 1), " months\n")
}

cat("\nOVERALL STATISTICS:\n")
cat("- Mean change point:", round(overall_mean_cp, 1), "months\n")
cat("- Critical window (mean ± 1SD):", round(overall_mean_cp - overall_sd_cp, 1), 
    "to", round(overall_mean_cp + overall_sd_cp, 1), "months\n")
cat("- Mean maximum divergence time:", round(overall_mean_divergence, 1), "months\n")

cat("\nGENERATED FILES:\n")
cat("- change_point_summary_corrected.csv: Summary results\n")
cat("- change_point_[metric]_corrected.png: Individual change point plots\n")
cat("- change_point_analysis_combined_corrected.png: Combined plot\n")
cat("- group_divergence_analysis_corrected.png: Divergence summary\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")