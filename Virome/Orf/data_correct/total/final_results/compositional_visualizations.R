#!/usr/bin/env Rscript
# Compositional Analysis Visualizations and Summary
# Author: Claude AI
# Date: 2025-07-25

library(dplyr)
library(vegan)
library(ggplot2)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(ggrepel)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

cat("=== COMPOSITIONAL ANALYSIS VISUALIZATIONS ===\n")
cat("Starting visualization at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD ALL RESULTS
# ============================================================================

cat("1. Loading all results...\n")

# Load all results
diversity_data_full <- read.csv("diversity_data_full.csv")
slope_data <- read.csv("slope_data.csv")
stability_data <- read.csv("stability_data.csv")
turnover_data <- read.csv("turnover_data.csv")

trajectory_results <- read.csv("diversity_trajectory_results.csv")
slope_results <- read.csv("slope_analysis_results.csv")
stability_results <- read.csv("stability_analysis_results.csv")
turnover_results <- read.csv("turnover_analysis_results.csv")

cat("All results loaded successfully\n")

# ============================================================================
# PART 2: CREATE COMPREHENSIVE RESULTS SUMMARY
# ============================================================================

cat("\n2. Creating comprehensive results summary...\n")

# Combine all results into comprehensive summary
create_comprehensive_summary <- function() {
  
  # Trajectory results
  traj_summary <- trajectory_results %>%
    select(metric, logFC, P.Value, adj.P.Val, n_samples, n_patients) %>%
    mutate(analysis_type = "Trajectory_Interaction", effect_size = logFC)
  
  # Slope results
  slope_summary <- slope_results %>%
    select(metric, logFC, P.Value, adj.P.Val, n_patients, celiac_mean_slope, control_mean_slope) %>%
    mutate(analysis_type = "Slope_Difference", effect_size = logFC) %>%
    rename(n_samples = n_patients)
  
  # Stability results
  stability_summary <- stability_results %>%
    select(metric, logFC, P.Value, adj.P.Val, n_patients, celiac_mean_stability, control_mean_stability) %>%
    mutate(analysis_type = "Stability_Difference", effect_size = logFC) %>%
    rename(n_samples = n_patients)
  
  # Turnover results
  turnover_summary <- turnover_results %>%
    select(metric, logFC, P.Value, adj.P.Val, n_patients, celiac_mean_turnover, control_mean_turnover) %>%
    mutate(analysis_type = "Turnover_Difference", effect_size = logFC) %>%
    rename(n_samples = n_patients)
  
  # Combine all
  comprehensive_results <- bind_rows(
    traj_summary,
    slope_summary,
    stability_summary,
    turnover_summary
  )
  
  return(comprehensive_results)
}

comprehensive_results <- create_comprehensive_summary()
write.csv(comprehensive_results, "comprehensive_results_summary.csv", row.names = FALSE)

cat("Comprehensive results summary created\n")

# ============================================================================
# PART 3: DIVERSITY TRAJECTORY HEATMAP
# ============================================================================

cat("\n3. Creating diversity trajectory heatmap...\n")

# Create heatmap of trajectory results
create_trajectory_heatmap <- function(results) {
  
  # Prepare data for heatmap
  heatmap_data <- results %>%
    select(metric, logFC, adj.P.Val) %>%
    mutate(
      neg_log_pval = -log10(adj.P.Val),
      significance = ifelse(adj.P.Val < 0.05, "Significant", "Not Significant")
    )
  
  # Create matrix for heatmap
  heatmap_matrix <- matrix(heatmap_data$logFC, 
                          nrow = 1, 
                          dimnames = list("Dx.Status x Timeline", heatmap_data$metric))
  
  # Create heatmap
  png("diversity_trajectory_heatmap.png", width = 1200, height = 400, res = 150)
  
  pheatmap(heatmap_matrix,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           scale = "none",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "Diversity Trajectory Analysis: Dx.Status × Timeline Interactions\n(Effect sizes from limma analysis)",
           fontsize = 10,
           cellwidth = 80,
           cellheight = 40,
           display_numbers = TRUE,
           number_format = "%.3f")
  
  dev.off()
  
  cat("Trajectory heatmap saved\n")
}

create_trajectory_heatmap(trajectory_results)

# ============================================================================
# PART 4: SLOPE ANALYSIS VISUALIZATION
# ============================================================================

cat("\n4. Creating slope analysis visualizations...\n")

# Bar plot of slope differences
create_slope_barplot <- function(results) {
  
  # Prepare data
  plot_data <- results %>%
    mutate(
      significance = ifelse(adj.P.Val < 0.05, "Significant", "Not Significant"),
      metric_clean = gsub("_", " ", stringr::str_to_title(metric))
    )
  
  p <- ggplot(plot_data, aes(x = reorder(metric_clean, -abs(logFC)), y = logFC)) +
    geom_col(aes(fill = significance), alpha = 0.8) +
    geom_text(aes(label = paste("p =", round(adj.P.Val, 4))), 
              vjust = ifelse(plot_data$logFC > 0, -0.5, 1.5), size = 3) +
    scale_fill_manual(values = c("Significant" = "#E31A1C", "Not Significant" = "grey70")) +
    labs(
      title = "Slope Analysis: Disease Group Differences",
      subtitle = "Effect sizes from patient-level slope comparisons (CELIAC vs CONTROL)",
      x = "Diversity Metric",
      y = "Effect Size (logFC)",
      fill = "Significance"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave("slope_analysis_results.png", p, width = 10, height = 6, dpi = 300)
  
  cat("Slope bar plot saved\n")
}

create_slope_barplot(slope_results)

# Volcano plot of slope analysis
create_slope_volcano <- function(results) {
  
  plot_data <- results %>%
    mutate(
      neg_log_pval = -log10(adj.P.Val),
      significance = case_when(
        adj.P.Val < 0.01 ~ "Highly Significant (p < 0.01)",
        adj.P.Val < 0.05 ~ "Significant (p < 0.05)",
        TRUE ~ "Not Significant"
      )
    )
  
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log_pval)) +
    geom_point(aes(color = significance), size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = metric), size = 3, max.overlaps = 10) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c(
      "Highly Significant (p < 0.01)" = "#B2182B",
      "Significant (p < 0.05)" = "#E31A1C",
      "Not Significant" = "grey70"
    )) +
    labs(
      title = "Volcano Plot: Slope Analysis Results",
      subtitle = "Patient-level slope differences between CELIAC and CONTROL groups",
      x = "Effect Size (logFC)",
      y = "-log10(adjusted P-value)",
      color = "Significance"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave("slope_volcano_plot.png", p, width = 10, height = 8, dpi = 300)
  
  cat("Slope volcano plot saved\n")
}

create_slope_volcano(slope_results)

# ============================================================================
# PART 5: DIVERSITY TRAJECTORIES BY GROUP
# ============================================================================

cat("\n5. Creating diversity trajectories by group...\n")

# Plot mean trajectories by disease status
create_trajectory_plots <- function(diversity_data) {
  
  # Calculate group means by timepoint
  trajectory_data <- diversity_data %>%
    group_by(Dx.Status, onset_timeline_numeric) %>%
    summarise(
      across(c(richness, shannon, simpson, evenness, total_abundance, dominance, viral_load_cv),
             list(mean = ~mean(.x, na.rm = TRUE), se = ~sd(.x, na.rm = TRUE)/sqrt(n())),
             .names = "{.col}_{.fn}"),
      n_samples = n(),
      .groups = "drop"
    )
  
  # Create plots for each metric
  metrics <- c("richness", "shannon", "simpson", "evenness")
  
  plot_list <- list()
  
  for(metric in metrics) {
    mean_col <- paste0(metric, "_mean")
    se_col <- paste0(metric, "_se")
    
    p <- ggplot(trajectory_data, aes(x = onset_timeline_numeric, y = .data[[mean_col]], color = Dx.Status)) +
      geom_line(size = 1, alpha = 0.8) +
      geom_point(size = 2) +
      geom_ribbon(aes(ymin = .data[[mean_col]] - .data[[se_col]], 
                      ymax = .data[[mean_col]] + .data[[se_col]],
                      fill = Dx.Status), alpha = 0.2, color = NA) +
      scale_color_manual(values = c("CELIAC" = "#E31A1C", "CONTROL" = "#1F78B4")) +
      scale_fill_manual(values = c("CELIAC" = "#E31A1C", "CONTROL" = "#1F78B4")) +
      labs(
        title = paste("Trajectory:", stringr::str_to_title(metric)),
        x = "Months Relative to Onset",
        y = stringr::str_to_title(metric),
        color = "Disease Status",
        fill = "Disease Status"
      ) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot_list[[metric]] <- p
  }
  
  # Combine plots
  combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)
  
  ggsave("diversity_trajectories_by_group.png", combined_plot, width = 12, height = 10, dpi = 300)
  
  cat("Trajectory plots saved\n")
}

create_trajectory_plots(diversity_data_full)

# ============================================================================
# PART 6: SUMMARY STATISTICS VISUALIZATION
# ============================================================================

cat("\n6. Creating summary statistics visualization...\n")

create_summary_stats_plot <- function() {
  
  # Count significant results by analysis type
  sig_counts <- comprehensive_results %>%
    group_by(analysis_type) %>%
    summarise(
      total_tests = n(),
      significant_05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
      significant_01 = sum(adj.P.Val < 0.01, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      prop_sig_05 = significant_05 / total_tests,
      prop_sig_01 = significant_01 / total_tests
    )
  
  # Create bar plot
  plot_data <- sig_counts %>%
    select(analysis_type, significant_05, significant_01) %>%
    reshape2::melt(id.vars = "analysis_type", variable.name = "threshold", value.name = "count") %>%
    mutate(
      threshold = ifelse(threshold == "significant_05", "p < 0.05", "p < 0.01"),
      analysis_clean = gsub("_", " ", analysis_type)
    )
  
  p <- ggplot(plot_data, aes(x = analysis_clean, y = count, fill = threshold)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.5) +
    scale_fill_manual(values = c("p < 0.05" = "#FB9A99", "p < 0.01" = "#E31A1C")) +
    labs(
      title = "Compositional Analysis Summary",
      subtitle = "Number of significant results by analysis type",
      x = "Analysis Type",
      y = "Number of Significant Results",
      fill = "Significance Threshold"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave("analysis_summary_statistics.png", p, width = 10, height = 6, dpi = 300)
  
  cat("Summary statistics plot saved\n")
}

create_summary_stats_plot()

# ============================================================================
# PART 7: CREATE FINAL SUMMARY REPORT
# ============================================================================

cat("\n7. Creating final summary report...\n")

# Print comprehensive summary
cat("\n=== COMPOSITIONAL ANALYSIS RESULTS SUMMARY ===\n")

cat("\nDATASET OVERVIEW:\n")
cat("- Total ORFs:", nrow(read.csv("total_orf.abundance.clean.csv")), "\n")
cat("- Total samples:", length(unique(diversity_data_full$sample_id)), "\n")
cat("- Total patients:", length(unique(diversity_data_full$patientID)), "\n")
cat("- CELIAC samples:", sum(diversity_data_full$Dx.Status == "CELIAC"), "\n")
cat("- CONTROL samples:", sum(diversity_data_full$Dx.Status == "CONTROL"), "\n")

cat("\nTRAJECTORY ANALYSIS RESULTS:\n")
sig_traj <- trajectory_results[trajectory_results$adj.P.Val < 0.05, ]
if(nrow(sig_traj) > 0) {
  cat("Significant Dx.Status × Timeline interactions (adj.p < 0.05):\n")
  for(i in 1:nrow(sig_traj)) {
    cat("- ", sig_traj$metric[i], ": effect =", round(sig_traj$logFC[i], 4), 
        ", p =", format(sig_traj$adj.P.Val[i], scientific = TRUE), "\n")
  }
} else {
  cat("No significant trajectory interactions found\n")
}

cat("\nSLOPE ANALYSIS RESULTS:\n")
sig_slope <- slope_results[slope_results$adj.P.Val < 0.05, ]
if(nrow(sig_slope) > 0) {
  cat("Significant slope differences (adj.p < 0.05):\n")
  for(i in 1:nrow(sig_slope)) {
    cat("- ", sig_slope$metric[i], ": effect =", round(sig_slope$logFC[i], 4), 
        ", p =", format(sig_slope$adj.P.Val[i], scientific = TRUE), "\n")
  }
} else {
  cat("No significant slope differences found\n")
}

cat("\nSTABILITY ANALYSIS RESULTS:\n")
sig_stab <- stability_results[stability_results$adj.P.Val < 0.05, ]
if(nrow(sig_stab) > 0) {
  cat("Significant stability differences (adj.p < 0.05):\n")
  for(i in 1:nrow(sig_stab)) {
    cat("- ", sig_stab$metric[i], ": effect =", round(sig_stab$logFC[i], 4), 
        ", p =", format(sig_stab$adj.P.Val[i], scientific = TRUE), "\n")
  }
} else {
  cat("No significant stability differences found\n")
}

cat("\nTURNOVER ANALYSIS RESULTS:\n")
sig_turn <- turnover_results[turnover_results$adj.P.Val < 0.05, ]
if(nrow(sig_turn) > 0) {
  cat("Significant turnover differences (adj.p < 0.05):\n")
  for(i in 1:nrow(sig_turn)) {
    cat("- ", sig_turn$metric[i], ": effect =", round(sig_turn$logFC[i], 4), 
        ", p =", format(sig_turn$adj.P.Val[i], scientific = TRUE), "\n")
  }
} else {
  cat("No significant turnover differences found\n")
}

cat("\n=== COMPOSITIONAL ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- comprehensive_results_summary.csv: All analysis results combined\n")
cat("- diversity_trajectory_heatmap.png: Heatmap of trajectory results\n")
cat("- slope_analysis_results.png: Bar plot of slope differences\n")
cat("- slope_volcano_plot.png: Volcano plot of slope analysis\n")
cat("- diversity_trajectories_by_group.png: Mean trajectories by disease status\n")
cat("- analysis_summary_statistics.png: Summary statistics plots\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")