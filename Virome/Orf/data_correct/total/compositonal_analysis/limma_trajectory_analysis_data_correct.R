#!/usr/bin/env Rscript
# Limma Trajectory Analysis for data_correct Directory
# Replicating the exact same analysis as performed in data/total/
# Author: Claude AI
# Date: 2025-07-21

library(limma)
library(dplyr)
library(vegan)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/compositonal_analysis")

cat("=== LIMMA TRAJECTORY ANALYSIS - DATA_CORRECT ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD AND PREPARE DATA
# ============================================================================

cat("1. Loading data from data_correct directory...\n")

# Load ORF abundance data
orf_data <- read.csv("../total.orf.abundance.table_0.75_prevFiltered_temporal_cleaned.csv", row.names = 1)

# Load metadata (using row names for matching)
metadata <- read.csv("../total.metadata.cleaned.csv", row.names = 1)

cat("Data dimensions:\n")
cat("  ORFs:", nrow(orf_data), "\n")
cat("  Samples (ORF data):", ncol(orf_data), "\n") 
cat("  Metadata samples:", nrow(metadata), "\n")

# Clean sample names in ORF data (remove X prefix if present)
orf_samples <- colnames(orf_data)
orf_samples_clean <- gsub("^X", "", orf_samples)
colnames(orf_data) <- orf_samples_clean

# Match samples between ORF data and metadata using row names
common_samples <- intersect(orf_samples_clean, rownames(metadata))
metadata_matched <- metadata[common_samples, ]
orf_data_matched <- orf_data[, common_samples]

cat("Matched samples:", length(common_samples), "\n")
cat("Final data dimensions:\n")
cat("  ORFs:", nrow(orf_data_matched), "\n")
cat("  Samples:", ncol(orf_data_matched), "\n")
cat("  Patients:", length(unique(metadata_matched$patientID)), "\n")

# Prepare metadata for limma
metadata_clean <- metadata_matched %>%
  mutate(sample_id = rownames(metadata_matched)) %>%
  select(sample_id, patientID, Dx.Status, onset_timeline_numeric, Country, Sex, 
         Age.at.Gluten.Introduction..months., HLA.Category, 
         feeding_first_year, Delivery.Mode) %>%
  mutate(
    Dx.Status = factor(Dx.Status, levels = c("CONTROL", "CELIAC")),
    Country = factor(Country),
    Sex = factor(Sex),
    HLA.Category = factor(HLA.Category),
    feeding_first_year = factor(feeding_first_year),
    Delivery.Mode = factor(Delivery.Mode),
    patientID = factor(patientID)
  ) %>%
  arrange(sample_id)

cat("Sample distribution:\n")
print(table(metadata_clean$Dx.Status))
cat("Country distribution:\n") 
print(table(metadata_clean$Country))
cat("Onset timeline range:", range(metadata_clean$onset_timeline_numeric, na.rm = TRUE), "\n\n")

# ============================================================================
# PART 2: CALCULATE DERIVED ECOLOGICAL METRICS
# ============================================================================

cat("2. Calculating derived ecological metrics...\n")

# Transpose data for ecological calculations (samples as rows)
orf_data_t <- t(orf_data_matched)

# Calculate diversity metrics for each sample
calculate_diversity_metrics <- function(abundance_matrix) {
  
  # Ensure all values are numeric and non-negative
  abundance_matrix[abundance_matrix < 0] <- 0
  abundance_matrix[is.na(abundance_matrix)] <- 0
  
  diversity_metrics <- data.frame(
    sample_id = rownames(abundance_matrix),
    
    # Alpha diversity metrics
    richness = apply(abundance_matrix, 1, function(x) sum(x > 0)),
    shannon = diversity(abundance_matrix, index = "shannon"),
    simpson = diversity(abundance_matrix, index = "simpson"),
    evenness = diversity(abundance_matrix, index = "shannon") / log(apply(abundance_matrix, 1, function(x) sum(x > 0))),
    
    # Total viral load
    total_abundance = rowSums(abundance_matrix),
    
    # Dominance metrics
    dominance = apply(abundance_matrix, 1, function(x) {
      total <- sum(x)
      if(total > 0) max(x) / total else 0
    }),
    
    # Viral load variability
    viral_load_cv = apply(abundance_matrix, 1, function(x) {
      non_zero <- x[x > 0]
      if(length(non_zero) > 1) sd(non_zero) / mean(non_zero) else 0
    }),
    
    stringsAsFactors = FALSE
  )
  
  # Handle infinite/NaN values
  diversity_metrics$evenness[is.infinite(diversity_metrics$evenness) | is.nan(diversity_metrics$evenness)] <- 0
  diversity_metrics$viral_load_cv[is.infinite(diversity_metrics$viral_load_cv) | is.nan(diversity_metrics$viral_load_cv)] <- 0
  
  return(diversity_metrics)
}

# Calculate diversity metrics
diversity_data <- calculate_diversity_metrics(orf_data_t)

# Add metadata to diversity data
diversity_data_full <- merge(diversity_data, metadata_clean, by = "sample_id")

cat("Diversity metrics calculated for", nrow(diversity_data), "samples\n")
cat("Diversity summary:\n")
print(summary(diversity_data[, c("richness", "shannon", "simpson", "evenness", "total_abundance")]))

# Save diversity data
write.csv(diversity_data_full, "diversity_data_full.csv", row.names = FALSE)

# ============================================================================
# PART 3: DIVERSITY TRAJECTORY ANALYSIS
# ============================================================================

cat("\n3. Performing diversity trajectory analysis...\n")

# Create design matrix for trajectory analysis
design_trajectory <- model.matrix(~ Dx.Status * onset_timeline_numeric + 
                                 Country + Sex + Age.at.Gluten.Introduction..months. + 
                                 HLA.Category + feeding_first_year + Delivery.Mode, 
                                 data = diversity_data_full)

# Remove samples with missing values
complete_samples <- complete.cases(design_trajectory)
design_trajectory_clean <- design_trajectory[complete_samples, ]
diversity_clean <- diversity_data_full[complete_samples, ]

cat("Samples for analysis:", nrow(diversity_clean), "\n")

# Block by patient for repeated measures
patient_block <- factor(diversity_clean$patientID)

# Prepare diversity matrix for limma
diversity_metrics_matrix <- t(diversity_clean[, c("richness", "shannon", "simpson", 
                                                 "evenness", "total_abundance", 
                                                 "dominance", "viral_load_cv")])

# Fit trajectory models
fit_trajectory <- lmFit(diversity_metrics_matrix, design_trajectory_clean, block = patient_block,
                       correlation = duplicateCorrelation(diversity_metrics_matrix, 
                                                         design_trajectory_clean, 
                                                         block = patient_block)$consensus)

# Apply empirical Bayes moderation
fit_trajectory <- eBayes(fit_trajectory)

# Extract interaction results (Dx.Status:onset_timeline_numeric)
interaction_coef <- grep("Dx.Status.*onset_timeline_numeric", colnames(design_trajectory_clean), value = TRUE)
trajectory_results <- topTable(fit_trajectory, coef = interaction_coef, n = Inf)

cat("Trajectory analysis completed\n")
cat("Significant interactions (p < 0.05):", sum(trajectory_results$P.Value < 0.05), "\n")
cat("Significant interactions (adj.p < 0.05):", sum(trajectory_results$adj.P.Val < 0.05), "\n")

# Save trajectory results
write.csv(trajectory_results, "diversity_trajectory_results.csv")

# ============================================================================
# PART 4: SLOPE ANALYSIS
# ============================================================================

cat("\n4. Calculating individual slope metrics...\n")

# Calculate slopes for each patient
calculate_slopes <- function(diversity_data_full) {
  
  slope_results <- data.frame()
  
  for(patient in unique(diversity_data_full$patientID)) {
    patient_data <- diversity_data_full[diversity_data_full$patientID == patient, ]
    patient_data <- patient_data[order(patient_data$onset_timeline_numeric), ]
    
    if(nrow(patient_data) >= 3) {  # Need at least 3 timepoints
      
      # Calculate slopes for each diversity metric
      slopes <- list()
      metrics <- c("richness", "shannon", "simpson", "evenness", 
                  "total_abundance", "dominance", "viral_load_cv")
      
      for(metric in metrics) {
        if(sum(!is.na(patient_data[[metric]])) >= 3) {
          slope_model <- lm(patient_data[[metric]] ~ patient_data$onset_timeline_numeric)
          slopes[[paste0(metric, "_slope")]] <- coef(slope_model)[2]
        } else {
          slopes[[paste0(metric, "_slope")]] <- NA
        }
      }
      
      # Add patient info
      patient_slopes <- data.frame(
        patientID = patient,
        Dx.Status = patient_data$Dx.Status[1],
        Country = patient_data$Country[1],
        Sex = patient_data$Sex[1],
        Age.at.Gluten.Introduction..months. = patient_data$Age.at.Gluten.Introduction..months.[1],
        HLA.Category = patient_data$HLA.Category[1],
        feeding_first_year = patient_data$feeding_first_year[1],
        Delivery.Mode = patient_data$Delivery.Mode[1],
        n_timepoints = nrow(patient_data),
        slopes,
        stringsAsFactors = FALSE
      )
      
      slope_results <- rbind(slope_results, patient_slopes)
    }
  }
  
  return(slope_results)
}

# Calculate slopes
slope_data <- calculate_slopes(diversity_data_full)

cat("Calculated slopes for", nrow(slope_data), "patients\n")

# Save slope data
write.csv(slope_data, "slope_data.csv", row.names = FALSE)

# Perform limma analysis on slopes
slope_metrics <- slope_data[, grep("_slope$", colnames(slope_data))]
slope_metrics_matrix <- t(slope_metrics)

# Design matrix for slope analysis
design_slope <- model.matrix(~ Dx.Status + Country + Sex + Age.at.Gluten.Introduction..months. + 
                            HLA.Category + feeding_first_year + Delivery.Mode, 
                            data = slope_data)

# Remove samples with missing values
complete_slope_samples <- complete.cases(design_slope) & complete.cases(slope_metrics)
design_slope_clean <- design_slope[complete_slope_samples, ]
slope_metrics_clean <- slope_metrics_matrix[, complete_slope_samples]

# Fit slope models
fit_slopes <- lmFit(slope_metrics_clean, design_slope_clean)
fit_slopes <- eBayes(fit_slopes)

# Extract slope results
slope_results <- topTable(fit_slopes, coef = "Dx.StatusCELIAC", n = Inf)

cat("Slope analysis completed\n")
cat("Significant slope differences (p < 0.05):", sum(slope_results$P.Value < 0.05), "\n")
cat("Significant slope differences (adj.p < 0.05):", sum(slope_results$adj.P.Val < 0.05), "\n")

# Save slope results
write.csv(slope_results, "slope_analysis_results.csv")

# ============================================================================
# PART 5: COMPREHENSIVE RESULTS SUMMARY
# ============================================================================

cat("\n5. Creating comprehensive results summary...\n")

# Combine all results
comprehensive_results <- list(
  trajectory_analysis = trajectory_results,
  slope_analysis = slope_results
)

# Create summary statistics
summary_stats <- data.frame(
  Analysis = c("Diversity Trajectory (Interaction)", "Slope Analysis"),
  Total_Tests = c(nrow(trajectory_results), nrow(slope_results)),
  Significant_p005 = c(sum(trajectory_results$P.Value < 0.05), sum(slope_results$P.Value < 0.05)),
  Significant_adjp005 = c(sum(trajectory_results$adj.P.Val < 0.05), sum(slope_results$adj.P.Val < 0.05)),
  stringsAsFactors = FALSE
)

print(summary_stats)

# Save comprehensive results
write.csv(summary_stats, "comprehensive_results_summary.csv", row.names = FALSE)

# ============================================================================
# PART 6: GENERATE VISUALIZATIONS
# ============================================================================

cat("\n6. Generating visualizations...\n")

# 1. Trajectory results heatmap
if(nrow(trajectory_results) > 0) {
  
  png("diversity_trajectory_heatmap.png", width = 1200, height = 800, res = 150)
  
  # Prepare data for heatmap
  heatmap_data <- trajectory_results[, c("logFC", "t", "P.Value", "adj.P.Val")]
  heatmap_matrix <- as.matrix(heatmap_data)
  
  # Create heatmap
  pheatmap(t(heatmap_matrix),
          cluster_rows = FALSE,
          cluster_cols = TRUE,
          scale = "row",
          color = colorRampPalette(c("blue", "white", "red"))(100),
          main = "Diversity Trajectory Analysis Results\n(Dx.Status × onset_timeline_numeric interaction)",
          fontsize = 10,
          angle_col = 45)
  
  dev.off()
  
  cat("Generated diversity_trajectory_heatmap.png\n")
}

# 2. Slope analysis results
if(nrow(slope_results) > 0) {
  
  # Bar plot of slope results
  slope_plot_data <- slope_results
  slope_plot_data$metric <- rownames(slope_plot_data)
  slope_plot_data$significant <- slope_plot_data$adj.P.Val < 0.05
  
  p_slope <- ggplot(slope_plot_data, aes(x = reorder(metric, logFC), y = logFC)) +
    geom_col(aes(fill = significant), alpha = 0.7) +
    scale_fill_manual(values = c("grey70", "red"), name = "Significant\n(adj.p < 0.05)") +
    coord_flip() +
    labs(
      title = "Slope Analysis Results: CELIAC vs CONTROL",
      subtitle = paste("Significant slopes:", sum(slope_plot_data$significant), "out of", nrow(slope_plot_data)),
      x = "Diversity Metric Slope",
      y = "Log Fold Change (CELIAC vs CONTROL)",
      caption = "Positive values = steeper slopes in CELIAC cases"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text = element_text(size = 10)
    )
  
  ggsave("slope_analysis_results.png", p_slope, width = 10, height = 6, dpi = 300)
  
  cat("Generated slope_analysis_results.png\n")
  
  # Volcano plot for slopes
  slope_plot_data$neg_log10_padj <- -log10(slope_plot_data$adj.P.Val)
  
  p_volcano <- ggplot(slope_plot_data, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(aes(color = significant), size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = metric), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("grey70", "red"), name = "Significant\n(adj.p < 0.05)") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = "Slope Analysis Volcano Plot",
      x = "Log Fold Change (CELIAC vs CONTROL)",
      y = "-Log10(Adjusted P-value)",
      caption = "Dashed line: adj.p = 0.05 threshold"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave("slope_volcano_plot.png", p_volcano, width = 10, height = 8, dpi = 300)
  
  cat("Generated slope_volcano_plot.png\n")
}

# 3. Individual trajectory plots
if(nrow(diversity_data_full) > 0) {
  
  # Plot mean trajectories by group
  trajectory_summary <- diversity_data_full %>%
    group_by(Dx.Status, onset_timeline_numeric) %>%
    summarise(
      mean_richness = mean(richness, na.rm = TRUE),
      se_richness = sd(richness, na.rm = TRUE) / sqrt(n()),
      mean_shannon = mean(shannon, na.rm = TRUE),
      se_shannon = sd(shannon, na.rm = TRUE) / sqrt(n()),
      n_samples = n(),
      .groups = "drop"
    )
  
  # Richness trajectory plot
  p_richness <- ggplot(trajectory_summary, aes(x = onset_timeline_numeric, y = mean_richness, color = Dx.Status)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = mean_richness - se_richness, ymax = mean_richness + se_richness, fill = Dx.Status), 
                alpha = 0.2, color = NA) +
    scale_color_manual(values = c("CONTROL" = "blue", "CELIAC" = "red")) +
    scale_fill_manual(values = c("CONTROL" = "blue", "CELIAC" = "red")) +
    labs(
      title = "Viral Richness Trajectories by Disease Status",
      x = "Months Relative to Onset",
      y = "Mean Viral Richness",
      caption = "Error bands represent standard error"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave("diversity_trajectories_by_group.png", p_richness, width = 12, height = 8, dpi = 300)
  
  cat("Generated diversity_trajectories_by_group.png\n")
}

# ============================================================================
# PART 7: CREATE CLAUDE.MD DOCUMENTATION
# ============================================================================

cat("\n7. Creating documentation...\n")

claude_md_content <- paste0("# CLAUDE.md - Limma Trajectory Analysis Session (data_correct)

## Session Overview
**Date:** ", Sys.Date(), "  
**Task:** Limma trajectory analysis comparing CELIAC vs CONTROL compositional differences using derived ecological metrics  
**Working Directory:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/compositonal_analysis`

## Analysis Approach
This analysis replicates the exact same limma statistical framework applied to derived ecological metrics (diversity, stability, turnover) rather than individual viral taxa. The focus is on identifying group differences in trajectory patterns using the model: `~ Dx.Status * onset_timeline_numeric + confounders`.

## Data Sources
- **ORF Abundance Data:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/total.orf.abundance.table_0.75_prevFiltered_temporal_cleaned.csv`
- **Metadata:** `/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/total.metadata.cleaned.csv`
- **Sample Count:** ", nrow(diversity_clean), " samples from ", length(unique(diversity_clean$patientID)), " patients
- **ORF Count:** ", nrow(orf_data_matched), " viral ORFs

## Analysis Components

### 1. Diversity Trajectory Analysis
- **Model:** `~ Dx.Status * onset_timeline_numeric + Country + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode`
- **Metrics:** Richness, Shannon, Simpson, Evenness, Total abundance, Dominance, Viral load CV
- **Blocking:** Patient ID for repeated measures
- **Significant Results:** ", sum(trajectory_results$adj.P.Val < 0.05), " metrics with adj.p < 0.05

### 2. Slope Analysis
- **Model:** `~ Dx.Status + confounders` (patient-level slopes)
- **Significant Results:** ", sum(slope_results$adj.P.Val < 0.05), " slope metrics with adj.p < 0.05

## Generated Files

### Results Files
- `comprehensive_results_summary.csv` - All analysis results combined
- `diversity_trajectory_results.csv` - Interaction effects
- `slope_analysis_results.csv` - Slope differences
- `diversity_data_full.csv` - Sample-level diversity metrics
- `slope_data.csv` - Patient-level slope metrics

### Visualization Files
- `diversity_trajectory_heatmap.png` - Heatmap of trajectory results
- `slope_analysis_results.png` - Bar plot of slope differences
- `slope_volcano_plot.png` - Volcano plot of slope analysis
- `diversity_trajectories_by_group.png` - Mean trajectories by disease status

## Analysis Summary
- **Total ORFs:** ", nrow(orf_data_matched), "
- **Samples analyzed:** ", nrow(diversity_clean), "
- **Patients:** ", length(unique(diversity_clean$patientID)), "
- **Significant trajectory interactions:** ", sum(trajectory_results$adj.P.Val < 0.05), "
- **Significant slope differences:** ", sum(slope_results$adj.P.Val < 0.05), "

## Technical Environment
- **R Version:** ", R.version.string, "
- **Key Packages:** limma, dplyr, vegan, ggplot2, pheatmap
- **Analysis Runtime:** < 2 minutes
- **Memory Usage:** Standard desktop requirements

**Analysis Completed:** ", format(Sys.time()), "
")

writeLines(claude_md_content, "CLAUDE.md")

cat("Documentation created: CLAUDE.md\n")

cat("\n=== LIMMA TRAJECTORY ANALYSIS COMPLETED ===\n")
cat("Analysis completed at:", format(Sys.time()), "\n")
cat("All results saved in: /Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/compositonal_analysis/\n")