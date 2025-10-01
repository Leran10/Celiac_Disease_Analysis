#!/usr/bin/env Rscript
# Change Point Analysis with Group Trajectories - Corrected Data
# Author: Claude AI  
# Date: 2025-07-29

library(dplyr)
library(ggplot2)
library(changepoint)
library(gridExtra)
library(stringr)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

cat("=== CHANGE POINT TRAJECTORY ANALYSIS (CORRECTED DATA) ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD AND PREPARE DATA
# ============================================================================

cat("1. Loading and preparing trajectory data...\n")

# Load diversity data
diversity_data <- read.csv("diversity_data_full.csv")

# Calculate group means and standard errors by timepoint
diversity_metrics <- c("richness", "shannon", "simpson", "evenness", "total_abundance", "dominance", "viral_load_cv")

trajectory_data <- diversity_data %>%
  group_by(Dx.Status, onset_timeline_numeric) %>%
  summarise(
    across(all_of(diversity_metrics), 
           list(mean = ~mean(.x, na.rm = TRUE), 
                se = ~sd(.x, na.rm = TRUE)/sqrt(n())),
           .names = "{.col}_{.fn}"),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(onset_timeline_numeric)

cat("Prepared trajectory data with", nrow(trajectory_data), "group-timepoint combinations\n")

# ============================================================================
# PART 2: CHANGE POINT DETECTION FOR INDIVIDUAL GROUPS
# ============================================================================

cat("\n2. Performing change point detection for individual group trajectories...\n")

change_point_results <- list()

for(metric in diversity_metrics) {
  cat("  Analyzing", metric, "trajectories...\n")
  
  metric_mean <- paste0(metric, "_mean")
  
  # Get data for each group separately
  celiac_data <- trajectory_data[trajectory_data$Dx.Status == "CELIAC", ]
  control_data <- trajectory_data[trajectory_data$Dx.Status == "CONTROL", ]
  
  # Ensure both groups have data and are ordered by time
  celiac_data <- celiac_data[order(celiac_data$onset_timeline_numeric), ]
  control_data <- control_data[order(control_data$onset_timeline_numeric), ]
  
  if(nrow(celiac_data) > 5 && nrow(control_data) > 5) {
    
    # Change point detection for CELIAC group
    celiac_values <- celiac_data[[metric_mean]]
    celiac_times <- celiac_data$onset_timeline_numeric
    celiac_change_points <- numeric(0)
    
    if(length(celiac_values) > 5 && !any(is.na(celiac_values))) {
      tryCatch({
        # Use more conservative approach - limit to max 2 change points
        cpt_celiac <- cpt.mean(celiac_values, method="BinSeg", Q=2, penalty="SIC")
        if(length(cpts(cpt_celiac)) > 0) {
          celiac_cp_indices <- cpts(cpt_celiac)
          # Only keep change points that are not at the boundaries
          valid_indices <- celiac_cp_indices[celiac_cp_indices > 1 & celiac_cp_indices < length(celiac_values)]
          if(length(valid_indices) > 0) {
            celiac_change_points <- celiac_times[valid_indices]
          }
        }
      }, error = function(e) {
        # Fallback: find most significant change point
        if(length(celiac_values) > 3) {
          diffs <- abs(diff(celiac_values))
          if(max(diffs) > sd(celiac_values, na.rm = TRUE)) {  # Only if change is substantial
            max_change_idx <- which.max(diffs) + 1
            if(max_change_idx > 1 && max_change_idx < length(celiac_values)) {
              celiac_change_points <- celiac_times[max_change_idx]
            }
          }
        }
      })
    }
    
    # Change point detection for CONTROL group  
    control_values <- control_data[[metric_mean]]
    control_times <- control_data$onset_timeline_numeric
    control_change_points <- numeric(0)
    
    if(length(control_values) > 5 && !any(is.na(control_values))) {
      tryCatch({
        # Use more conservative approach - limit to max 2 change points
        cpt_control <- cpt.mean(control_values, method="BinSeg", Q=2, penalty="SIC")
        if(length(cpts(cpt_control)) > 0) {
          control_cp_indices <- cpts(cpt_control)
          # Only keep change points that are not at the boundaries
          valid_indices <- control_cp_indices[control_cp_indices > 1 & control_cp_indices < length(control_values)]
          if(length(valid_indices) > 0) {
            control_change_points <- control_times[valid_indices]
          }
        }
      }, error = function(e) {
        # Fallback: find most significant change point
        if(length(control_values) > 3) {
          diffs <- abs(diff(control_values))
          if(max(diffs) > sd(control_values, na.rm = TRUE)) {  # Only if change is substantial
            max_change_idx <- which.max(diffs) + 1
            if(max_change_idx > 1 && max_change_idx < length(control_values)) {
              control_change_points <- control_times[max_change_idx]
            }
          }
        }
      })
    }
    
    # Calculate maximum divergence point
    common_times <- intersect(celiac_times, control_times)
    max_div_time <- NA
    max_div_value <- 0
    
    if(length(common_times) > 1) {
      divergences <- sapply(common_times, function(t) {
        c_val <- celiac_data[celiac_data$onset_timeline_numeric == t, metric_mean]
        ct_val <- control_data[control_data$onset_timeline_numeric == t, metric_mean]
        if(length(c_val) > 0 && length(ct_val) > 0) {
          abs(c_val - ct_val)
        } else {
          0
        }
      })
      
      max_div_idx <- which.max(divergences)
      max_div_time <- common_times[max_div_idx]
      max_div_value <- divergences[max_div_idx]
    }
    
    # Store results
    change_point_results[[metric]] <- list(
      celiac_data = celiac_data,
      control_data = control_data,
      celiac_change_points = celiac_change_points,
      control_change_points = control_change_points,
      max_divergence_time = max_div_time,
      max_divergence_value = max_div_value,
      metric_mean = metric_mean,
      metric_se = paste0(metric, "_se")
    )
  }
}

cat("Change point detection completed\n")

# ============================================================================
# PART 3: CREATE TRAJECTORY PLOTS WITH CHANGE POINTS
# ============================================================================

cat("\n3. Creating trajectory plots with change points...\n")

create_trajectory_change_point_plot <- function(metric, results) {
  
  # Combine data for plotting
  plot_data <- rbind(
    results$celiac_data %>% select(Dx.Status, onset_timeline_numeric, 
                                   mean = all_of(results$metric_mean),
                                   se = all_of(results$metric_se)),
    results$control_data %>% select(Dx.Status, onset_timeline_numeric,
                                    mean = all_of(results$metric_mean), 
                                    se = all_of(results$metric_se))
  )
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = onset_timeline_numeric, y = mean, color = Dx.Status)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                  width = 2, alpha = 0.6) +
    
    # Add change point lines for CELIAC (red dashed)
    {if(length(results$celiac_change_points) > 0) {
      geom_vline(xintercept = results$celiac_change_points, 
                 color = "#d62728", linetype = "dashed", size = 1, alpha = 0.8)
    }} +
    
    # Add change point lines for CONTROL (blue dashed)  
    {if(length(results$control_change_points) > 0) {
      geom_vline(xintercept = results$control_change_points,
                 color = "#1f77b4", linetype = "dashed", size = 1, alpha = 0.8)
    }} +
    
    # Add maximum divergence line (black dotted)
    {if(!is.na(results$max_divergence_time)) {
      geom_vline(xintercept = results$max_divergence_time,
                 color = "black", linetype = "dotted", size = 1, alpha = 0.8)
    }} +
    
    # Colors and styling  
    scale_color_manual(values = c("CELIAC" = "#d62728", "CONTROL" = "#1f77b4")) +
    
    labs(
      title = paste("Change Point Analysis:", str_to_upper(metric)),
      subtitle = "Viral diversity trajectories with detected change points",
      x = "Time to Onset (months)",
      y = paste(str_to_title(gsub("_", " ", metric))),
      color = "Disease Status"
    ) +
    
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  # Add annotations for change points
  y_max <- max(plot_data$mean + plot_data$se, na.rm = TRUE)
  y_range <- max(plot_data$mean, na.rm = TRUE) - min(plot_data$mean, na.rm = TRUE)
  annotation_height <- y_max + 0.1 * y_range
  
  # Annotate CELIAC change points
  if(length(results$celiac_change_points) > 0) {
    for(i in seq_along(results$celiac_change_points)) {
      cp <- results$celiac_change_points[i]
      p <- p + annotate("text", x = cp, y = annotation_height,
                       label = paste("CELIAC CP:", cp),
                       color = "#d62728", size = 3, hjust = 0.5)
    }
  }
  
  # Annotate CONTROL change points  
  if(length(results$control_change_points) > 0) {
    for(i in seq_along(results$control_change_points)) {
      cp <- results$control_change_points[i]
      p <- p + annotate("text", x = cp, y = annotation_height * 0.95,
                       label = paste("CONTROL CP:", cp), 
                       color = "#1f77b4", size = 3, hjust = 0.5)
    }
  }
  
  # Annotate maximum divergence
  if(!is.na(results$max_divergence_time)) {
    p <- p + annotate("text", x = results$max_divergence_time, y = annotation_height * 0.9,
                     label = paste("Max Div:", results$max_divergence_time),
                     color = "black", size = 3, hjust = 0.5)
  }
  
  return(p)
}

# Create individual plots
individual_plots <- list()
for(metric in names(change_point_results)) {
  cat("  Creating plot for", metric, "...\n")
  
  individual_plots[[metric]] <- create_trajectory_change_point_plot(metric, change_point_results[[metric]])
  
  # Save individual plot
  ggsave(
    filename = paste0("change_point_trajectory_", metric, "_corrected.png"),
    plot = individual_plots[[metric]],
    width = 10, height = 7, dpi = 300
  )
}

# Create combined plot (2x2 grid of top 4 metrics)
if(length(individual_plots) >= 4) {
  top_metrics <- c("richness", "shannon", "simpson", "evenness")
  available_metrics <- intersect(top_metrics, names(individual_plots))
  
  if(length(available_metrics) >= 4) {
    combined_plot <- grid.arrange(
      individual_plots[[available_metrics[1]]], individual_plots[[available_metrics[2]]],
      individual_plots[[available_metrics[3]]], individual_plots[[available_metrics[4]]],
      ncol = 2, nrow = 2
    )
    
    ggsave("change_point_trajectory_combined_corrected.png", combined_plot,
           width = 16, height = 12, dpi = 300)
  }
}

# ============================================================================
# PART 4: SUMMARY TABLE
# ============================================================================

cat("\n4. Creating summary table...\n")

# Create summary of change points
cp_summary <- data.frame()

for(metric in names(change_point_results)) {
  results <- change_point_results[[metric]]
  
  summary_row <- data.frame(
    metric = metric,
    celiac_change_points = paste(results$celiac_change_points, collapse = ", "),
    control_change_points = paste(results$control_change_points, collapse = ", "),
    max_divergence_time = ifelse(is.na(results$max_divergence_time), "NA", as.character(results$max_divergence_time)),
    max_divergence_value = ifelse(is.numeric(results$max_divergence_value), 
                                  round(results$max_divergence_value, 3), 
                                  as.character(results$max_divergence_value)),
    stringsAsFactors = FALSE
  )
  
  cp_summary <- rbind(cp_summary, summary_row)
}

write.csv(cp_summary, "change_point_trajectory_summary_corrected.csv", row.names = FALSE)

cat("\n=== CHANGE POINT TRAJECTORY ANALYSIS COMPLETED ===\n")

cat("\nSUMMARY OF RESULTS:\n")
for(i in 1:nrow(cp_summary)) {
  row <- cp_summary[i, ]
  cat("- ", row$metric, ":\n")
  cat("  CELIAC change points: ", row$celiac_change_points, "\n")
  cat("  CONTROL change points: ", row$control_change_points, "\n") 
  cat("  Max divergence: ", row$max_divergence_value, " at ", row$max_divergence_time, " months\n\n")
}

cat("GENERATED FILES:\n")
cat("- change_point_trajectory_[metric]_corrected.png: Individual trajectory plots\n")
cat("- change_point_trajectory_combined_corrected.png: Combined 2x2 plot\n")
cat("- change_point_trajectory_summary_corrected.csv: Summary table\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")