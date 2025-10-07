# Comprehensive ORF PA Analysis with Enhanced Visualizations
library(dplyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(tibble)
library(tidyr)
library(cowplot)

cat("Creating comprehensive ORF PA analysis and visualizations...\n")

# Load PA results
total_pa_results <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/total_PA_model1_glmmTMB_logit.csv")
US_pa_results <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/US_PA_model1_glmmTMB_logit.csv")
total_pa_timepoints <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/total_PA_model1_glmmTMB_logit_timepoint_specific_results.csv")
US_pa_timepoints <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/US_PA_model1_glmmTMB_logit_timepoint_specific_results.csv")

cat("Loaded PA results - Total:", nrow(total_pa_results), "US:", nrow(US_pa_results), "\n")
cat("Timepoint results - Total:", nrow(total_pa_timepoints), "US:", nrow(US_pa_timepoints), "\n")

# Create output directory
dir.create("../Orf_Contig_Phrog_compositional/orf_figures", showWarnings = FALSE, recursive = TRUE)

# 1. Enhanced PA Effect Size Distribution Plot
create_pa_effect_distribution <- function(pa_results, cohort_name) {
  # Get main effects
  main_effects <- pa_results %>%
    filter(term == "Dx.StatusCELIAC") %>%
    filter(!is.na(estimate))
  
  if(nrow(main_effects) == 0) return(NULL)
  
  # Create distribution plot
  p1 <- ggplot(main_effects, aes(x = estimate)) +
    geom_histogram(bins = 50, fill = "#28a745", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    labs(
      title = paste(cohort_name, "Cohort - PA Effect Size Distribution"),
      subtitle = "Log Odds Ratios: CELIAC vs CONTROL",
      x = "Log Odds Ratio (Effect Size)",
      y = "Number of ORFs"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  
  # Add significance annotation
  n_sig <- sum(main_effects$p.value < 0.05, na.rm = TRUE)
  n_total <- nrow(main_effects)
  
  p1 <- p1 + annotate("text", x = Inf, y = Inf, 
                      label = paste("Significant:", n_sig, "/", n_total, "ORFs"),
                      hjust = 1.1, vjust = 1.5, size = 4, color = "red")
  
  return(p1)
}

# 2. PA Temporal Progression Analysis
create_pa_temporal_progression <- function(timepoint_data, cohort_name) {
  if(nrow(timepoint_data) == 0) return(NULL)
  
  # Get significant timepoint effects
  sig_timepoints <- timepoint_data %>%
    filter(!is.na(p.value) & p.value < 0.05) %>%
    group_by(timepoint) %>%
    summarise(
      n_significant = n(),
      mean_effect = mean(estimate, na.rm = TRUE),
      median_effect = median(estimate, na.rm = TRUE),
      .groups = "drop"
    )
  
  if(nrow(sig_timepoints) == 0) return(NULL)
  
  # Order timepoints chronologically
  timepoint_order <- c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0")
  sig_timepoints$timepoint <- factor(sig_timepoints$timepoint, levels = timepoint_order)
  sig_timepoints <- sig_timepoints %>% filter(timepoint %in% timepoint_order)
  
  # Create progression plot
  p <- ggplot(sig_timepoints, aes(x = timepoint)) +
    geom_col(aes(y = n_significant), fill = "#28a745", alpha = 0.8) +
    geom_line(aes(y = mean_effect * 2, group = 1), color = "red", size = 1.2) +
    geom_point(aes(y = mean_effect * 2), color = "red", size = 3) +
    scale_y_continuous(
      name = "Number of Significant ORFs",
      sec.axis = sec_axis(~ . / 2, name = "Mean Effect Size (Log Odds)")
    ) +
    labs(
      title = paste(cohort_name, "Cohort - PA Temporal Progression"),
      subtitle = "Significant ORFs and Effect Sizes Across Disease Timeline",
      x = "Timepoint (months before diagnosis)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(color = "#28a745"),
      axis.title.y.right = element_text(color = "red")
    )
  
  return(p)
}

# 3. PA Significance Heatmap (Enhanced)
create_pa_significance_heatmap <- function(timepoint_data, cohort_name) {
  if(nrow(timepoint_data) == 0) return(NULL)
  
  # Get significant results with p-values
  sig_data <- timepoint_data %>%
    filter(!is.na(p.value) & p.value < 0.05) %>%
    mutate(
      neg_log_p = -log10(p.value),
      signed_significance = ifelse(estimate > 0, neg_log_p, -neg_log_p)
    )
  
  if(nrow(sig_data) == 0) return(NULL)
  
  # Prepare data for heatmap
  timepoint_order <- c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0")
  
  heatmap_data <- sig_data %>%
    select(gene, timepoint, signed_significance) %>%
    mutate(timepoint = factor(timepoint, levels = timepoint_order)) %>%
    filter(timepoint %in% timepoint_order) %>%
    pivot_wider(names_from = timepoint, values_from = signed_significance, values_fill = 0) %>%
    column_to_rownames("gene")
  
  if(nrow(heatmap_data) == 0 || ncol(heatmap_data) == 0) return(NULL)
  
  # Create enhanced heatmap
  pdf(paste0("../Orf_Contig_Phrog_compositional/orf_figures/", cohort_name, "_PA_significance_heatmap.pdf"), 
      width = 12, height = max(8, nrow(heatmap_data) * 0.3))
  
  pheatmap(
    as.matrix(heatmap_data),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",
    cluster_cols = FALSE,  # Keep chronological order
    cluster_rows = TRUE,
    fontsize = 10,
    main = paste(cohort_name, "Cohort - PA Significance Heatmap\nSigned -log10(p-value): Red=Positive Effect, Blue=Negative Effect"),
    na_col = "grey90",
    breaks = seq(-max(abs(heatmap_data), na.rm = TRUE), max(abs(heatmap_data), na.rm = TRUE), length.out = 101)
  )
  
  dev.off()
  
  cat("PA significance heatmap saved for", cohort_name, ":", nrow(heatmap_data), "genes\n")
  return(heatmap_data)
}

# 4. Combined PA Analysis Summary
create_combined_pa_summary <- function() {
  # Combine timepoint summaries
  total_summary <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/total_PA_timepoint_summary.csv") %>%
    mutate(cohort = "Total")
  US_summary <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/US_PA_timepoint_summary.csv") %>%
    mutate(cohort = "US")
  
  combined_summary <- bind_rows(total_summary, US_summary)
  
  # Create comprehensive summary plot
  p1 <- ggplot(combined_summary, aes(x = timepoint, y = significant, fill = cohort)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Total" = "#2E8B57", "US" = "#4682B4")) +
    labs(
      title = "ORF PA Analysis - Significant Results by Timepoint",
      subtitle = "Comparison Across Cohorts",
      x = "Timepoint (months before diagnosis)",
      y = "Number of Significant ORFs",
      fill = "Cohort"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  p2 <- ggplot(combined_summary, aes(x = timepoint, y = mean_effect, color = cohort, group = cohort)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Total" = "#2E8B57", "US" = "#4682B4")) +
    labs(
      title = "Mean Effect Sizes Across Timepoints",
      x = "Timepoint (months before diagnosis)",
      y = "Mean Log Odds Ratio",
      color = "Cohort"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Combine plots
  combined_plot <- plot_grid(p1, p2, ncol = 1, align = "v")
  ggsave("../Orf_Contig_Phrog_compositional/orf_figures/combined_PA_comprehensive_summary.pdf", 
         combined_plot, width = 12, height = 10, device = "pdf")
  
  cat("Combined PA summary plot saved\n")
  return(combined_plot)
}

# 5. PA Model Diagnostic Plots
create_pa_diagnostic_plots <- function(pa_results, cohort_name) {
  if(nrow(pa_results) == 0) return(NULL)
  
  # Convergence analysis
  convergence_summary <- pa_results %>%
    group_by(orf_id) %>%
    summarise(
      n_terms = n(),
      n_converged = sum(!is.na(p.value)),
      convergence_rate = n_converged / n_terms,
      .groups = "drop"
    )
  
  p1 <- ggplot(convergence_summary, aes(x = convergence_rate)) +
    geom_histogram(bins = 20, fill = "#28a745", alpha = 0.7, color = "white") +
    labs(
      title = paste(cohort_name, "- PA Model Convergence"),
      x = "Convergence Rate (fraction of terms with valid p-values)",
      y = "Number of ORFs"
    ) +
    theme_bw()
  
  # Effect size vs significance
  main_effects <- pa_results %>%
    filter(term == "Dx.StatusCELIAC") %>%
    filter(!is.na(estimate) & !is.na(p.value))
  
  if(nrow(main_effects) > 0) {
    p2 <- ggplot(main_effects, aes(x = estimate, y = -log10(p.value))) +
      geom_point(aes(color = p.value < 0.05), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), 
                         name = "Significant") +
      labs(
        title = paste(cohort_name, "- PA Effect Size vs Significance"),
        x = "Log Odds Ratio",
        y = "-log10(p-value)"
      ) +
      theme_bw()
  } else {
    p2 <- ggplot() + theme_void() + labs(title = "No main effects data")
  }
  
  # Combine diagnostic plots
  diagnostic_plot <- plot_grid(p1, p2, ncol = 2)
  ggsave(paste0("../Orf_Contig_Phrog_compositional/orf_figures/", cohort_name, "_PA_diagnostic_plots.pdf"), 
         diagnostic_plot, width = 14, height = 6, device = "pdf")
  
  cat("PA diagnostic plots saved for", cohort_name, "\n")
  return(diagnostic_plot)
}

# Generate all PA visualizations
cat("\n=== GENERATING COMPREHENSIVE PA VISUALIZATIONS ===\n")

# 1. Effect distributions
total_dist <- create_pa_effect_distribution(total_pa_results, "Total")
US_dist <- create_pa_effect_distribution(US_pa_results, "US")

if(!is.null(total_dist)) {
  ggsave("../Orf_Contig_Phrog_compositional/orf_figures/total_PA_effect_distribution.pdf", 
         total_dist, width = 10, height = 7, device = "pdf")
}
if(!is.null(US_dist)) {
  ggsave("../Orf_Contig_Phrog_compositional/orf_figures/US_PA_effect_distribution.pdf", 
         US_dist, width = 10, height = 7, device = "pdf")
}

# 2. Temporal progressions
total_temporal <- create_pa_temporal_progression(total_pa_timepoints, "Total")
US_temporal <- create_pa_temporal_progression(US_pa_timepoints, "US")

if(!is.null(total_temporal)) {
  ggsave("../Orf_Contig_Phrog_compositional/orf_figures/total_PA_temporal_progression.pdf", 
         total_temporal, width = 10, height = 7, device = "pdf")
}
if(!is.null(US_temporal)) {
  ggsave("../Orf_Contig_Phrog_compositional/orf_figures/US_PA_temporal_progression.pdf", 
         US_temporal, width = 10, height = 7, device = "pdf")
}

# 3. Significance heatmaps
total_sig_heatmap <- create_pa_significance_heatmap(total_pa_timepoints, "total")
US_sig_heatmap <- create_pa_significance_heatmap(US_pa_timepoints, "US")

# 4. Combined summary
combined_summary <- create_combined_pa_summary()

# 5. Diagnostic plots
total_diagnostic <- create_pa_diagnostic_plots(total_pa_results, "total")
US_diagnostic <- create_pa_diagnostic_plots(US_pa_results, "US")

# Create PA analysis summary statistics
cat("\n=== PA ANALYSIS SUMMARY STATISTICS ===\n")

# Total cohort stats
total_main <- total_pa_results %>% filter(term == "Dx.StatusCELIAC")
total_sig_main <- sum(total_main$p.value < 0.05, na.rm = TRUE)
total_sig_timepoints <- sum(total_pa_timepoints$p.value < 0.05, na.rm = TRUE)

cat("TOTAL COHORT PA ANALYSIS:\n")
cat("- Total ORFs analyzed:", length(unique(total_pa_results$orf_id)), "\n")
cat("- Main effects significant:", total_sig_main, "\n")
cat("- Timepoint effects significant:", total_sig_timepoints, "\n")
cat("- Effect size range:", round(min(total_main$estimate, na.rm=TRUE), 2), "to", round(max(total_main$estimate, na.rm=TRUE), 2), "\n")

# US cohort stats
US_main <- US_pa_results %>% filter(term == "Dx.StatusCELIAC")
US_sig_main <- sum(US_main$p.value < 0.05, na.rm = TRUE)
US_sig_timepoints <- sum(US_pa_timepoints$p.value < 0.05, na.rm = TRUE)

cat("\nUS COHORT PA ANALYSIS:\n")
cat("- Total ORFs analyzed:", length(unique(US_pa_results$orf_id)), "\n")
cat("- Main effects significant:", US_sig_main, "\n")
cat("- Timepoint effects significant:", US_sig_timepoints, "\n")
cat("- Effect size range:", round(min(US_main$estimate, na.rm=TRUE), 2), "to", round(max(US_main$estimate, na.rm=TRUE), 2), "\n")

# Save summary stats
pa_summary_stats <- data.frame(
  cohort = c("Total", "US"),
  orfs_analyzed = c(length(unique(total_pa_results$orf_id)), length(unique(US_pa_results$orf_id))),
  main_effects_significant = c(total_sig_main, US_sig_main),
  timepoint_effects_significant = c(total_sig_timepoints, US_sig_timepoints),
  effect_range_min = c(min(total_main$estimate, na.rm=TRUE), min(US_main$estimate, na.rm=TRUE)),
  effect_range_max = c(max(total_main$estimate, na.rm=TRUE), max(US_main$estimate, na.rm=TRUE))
)

write.csv(pa_summary_stats, "../Orf_Contig_Phrog_compositional/orf_results/PA_analysis_summary_stats.csv", row.names = FALSE)

cat("\n=== COMPREHENSIVE PA ANALYSIS COMPLETED ===\n")
cat("New PA figures generated:\n")
cat("- Effect distribution plots (2)\n")
cat("- Temporal progression plots (2)\n") 
cat("- Significance heatmaps (2)\n")
cat("- Combined summary plot (1)\n")
cat("- Diagnostic plots (2)\n")
cat("- Summary statistics saved\n")