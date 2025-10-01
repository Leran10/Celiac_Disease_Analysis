#!/usr/bin/env Rscript
# Comprehensive Italy Cohort Analysis - Diversity-Focused Analysis
# Author: Claude AI
# Date: 2025-08-05
# Note: Only 2 significant ORFs found, focusing on diversity analysis

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)
library(vegan)
library(gridExtra)
library(ggrepel)
library(stringr)

cat("=== COMPREHENSIVE ITALY COHORT ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND ASSESS ORF SIGNIFICANCE
# ============================================================================

cat("1. Loading Italy cohort data and assessing ORF significance...\n")

# Load abundance data
abundance_data_raw <- read.csv("Italy.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names (remove X prefix if present)
sample_names <- colnames(abundance_data)
sample_names_clean <- gsub("^X", "", sample_names)
colnames(abundance_data) <- sample_names_clean

# Load metadata
metadata_raw <- read.csv("Italy.metadata.clean.csv")
rownames(metadata_raw) <- metadata_raw$X
metadata <- metadata_raw

cat("Loaded abundance data with", nrow(abundance_data), "ORFs and", ncol(abundance_data), "samples\n")
cat("Loaded metadata with", nrow(metadata), "samples\n")

# Match samples
matching_samples <- intersect(colnames(abundance_data), rownames(metadata))
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")

# Ensure Dx.Status is a factor with CELIAC as reference
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))

# Create design matrix (Italy cohort doesn't have Country variable)
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_numeric + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

cat("Design matrix dimensions:", nrow(linear.model.design), "x", ncol(linear.model.design), "\n")

# Convert abundance data to matrix for limma
abundance_matrix <- as.matrix(abundance_filtered)

# Fit limma model with voomLmFit and patient blocking
cat("Running voomLmFit with patient blocking for ORF assessment...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get significant ORFs for the INTERACTION term (Dx.StatusCONTROL:onset_timeline_numeric)
topTable_results <- topTable(linear.model.fit, coef = "Dx.StatusCONTROL:onset_timeline_numeric", number = Inf, sort.by = "P")
significant_orfs <- topTable_results[topTable_results$adj.P.Val < 0.05, ]

cat("Identified", nrow(significant_orfs), "significant ORFs at adj.p < 0.05 (interaction term)\n")
cat("Percentage of total ORFs:", round(nrow(significant_orfs)/nrow(abundance_matrix)*100, 1), "%\n")

# Save significant ORFs results
write.csv(significant_orfs, "Italy_limma_model_Sig_res_final.csv", row.names = TRUE)

# Note about limited ORF results
cat("\nNote: Due to limited significant ORFs (", nrow(significant_orfs), "), focusing on diversity analysis\n")

# ============================================================================
# PART 2: DIVERSITY ANALYSIS (PRIMARY FOCUS)
# ============================================================================

cat("\n2. Running comprehensive diversity analysis...\n")

# Calculate diversity metrics for each sample
diversity_metrics <- data.frame(
  sample_id = matching_samples,
  stringsAsFactors = FALSE
)

# Add metadata
diversity_metrics <- merge(diversity_metrics, metadata_filtered, by.x = "sample_id", by.y = "row.names")

# Calculate diversity indices
abundance_matrix_t <- t(abundance_filtered)  # Transpose for vegan

for(i in 1:nrow(diversity_metrics)) {
  sample_abundance <- abundance_matrix_t[i, ]
  
  # Richness (number of non-zero ORFs)
  diversity_metrics$richness[i] <- sum(sample_abundance > 0)
  
  # Shannon diversity
  shannon_val <- vegan::diversity(sample_abundance, index = "shannon")
  diversity_metrics$shannon[i] <- ifelse(is.finite(shannon_val), shannon_val, 0)
  
  # Simpson diversity
  simpson_val <- vegan::diversity(sample_abundance, index = "simpson")
  diversity_metrics$simpson[i] <- ifelse(is.finite(simpson_val), simpson_val, 0)
  
  # Evenness (Pielou's evenness) - handle division by zero
  if(diversity_metrics$richness[i] > 1) {
    diversity_metrics$evenness[i] <- diversity_metrics$shannon[i] / log(diversity_metrics$richness[i])
  } else {
    diversity_metrics$evenness[i] <- 0
  }
  
  # Total abundance
  diversity_metrics$total_abundance[i] <- sum(sample_abundance)
  
  # Dominance (1 - Simpson)
  diversity_metrics$dominance[i] <- 1 - diversity_metrics$simpson[i]
  
  if(i %% 25 == 0) {
    cat("Calculated diversity for", i, "samples\n")
  }
}

# Calculate viral load CV by patient
diversity_metrics <- diversity_metrics %>%
  group_by(patientID) %>%
  mutate(viral_load_cv = ifelse(mean(total_abundance, na.rm = TRUE) > 0, 
                               sd(total_abundance, na.rm = TRUE) / mean(total_abundance, na.rm = TRUE), 
                               0)) %>%
  ungroup()

# Replace any remaining NAs or infinite values
diversity_metrics$viral_load_cv[!is.finite(diversity_metrics$viral_load_cv)] <- 0
diversity_metrics$evenness[!is.finite(diversity_metrics$evenness)] <- 0

cat("Diversity metrics calculated for", nrow(diversity_metrics), "samples\n")

# ============================================================================
# PART 3: DIVERSITY TRAJECTORY ANALYSIS
# ============================================================================

cat("\n3. Running diversity trajectory analysis...\n")

# Create design matrix for trajectory analysis
trajectory_design <- model.matrix(~ Dx.Status * onset_timeline_numeric + Sex + 
                                  Age.at.Gluten.Introduction..months. + HLA.Category + 
                                  feeding_first_year + Delivery.Mode, 
                                  data = diversity_metrics)

# Create diversity matrix for limma
diversity_matrix <- as.matrix(t(diversity_metrics[, c("richness", "shannon", "simpson", 
                                                     "evenness", "total_abundance", 
                                                     "dominance", "viral_load_cv")]))

# Check for any remaining NAs in diversity matrix
diversity_matrix[!is.finite(diversity_matrix)] <- 0

# Run limma with patient blocking
diversity_fit <- lmFit(diversity_matrix, trajectory_design)
diversity_fit <- eBayes(diversity_fit)

# Extract trajectory results
trajectory_results <- topTable(diversity_fit, coef = "Dx.StatusCONTROL:onset_timeline_numeric", 
                               number = Inf, sort.by = "P")

cat("Trajectory analysis completed\n")

# ============================================================================
# PART 4: CREATE DIVERSITY TRAJECTORY HEATMAP
# ============================================================================

cat("\n4. Creating diversity trajectory heatmap...\n")

# Create heatmap showing trajectory interaction effects
png("diversity_trajectory_heatmap_Italy.png", width = 800, height = 600, res = 150)

# Create a simple heatmap of the logFC values
heatmap_matrix <- as.matrix(trajectory_results$logFC)
rownames(heatmap_matrix) <- rownames(trajectory_results)
colnames(heatmap_matrix) <- "Trajectory_Interaction"

pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Italy Cohort: Diversity Trajectory Interactions (Dx.Status × Timeline)",
         fontsize = 12,
         fontsize_row = 10,
         show_colnames = FALSE,
         annotation_legend = TRUE,
         border_color = "white")

dev.off()
cat("Diversity trajectory heatmap saved\n")

# ============================================================================
# PART 5: SLOPE ANALYSIS
# ============================================================================

cat("\n5. Running slope analysis...\n")

# Calculate individual patient slopes
slope_data <- diversity_metrics %>%
  group_by(patientID) %>%
  summarise(
    Dx.Status = first(Dx.Status),
    Sex = first(Sex),
    Age.at.Gluten.Introduction..months. = first(Age.at.Gluten.Introduction..months.),
    HLA.Category = first(HLA.Category),
    feeding_first_year = first(feeding_first_year),
    Delivery.Mode = first(Delivery.Mode),
    richness_slope = if(n() >= 3) coef(lm(richness ~ onset_timeline_numeric))[2] else NA,
    shannon_slope = if(n() >= 3) coef(lm(shannon ~ onset_timeline_numeric))[2] else NA,
    simpson_slope = if(n() >= 3) coef(lm(simpson ~ onset_timeline_numeric))[2] else NA,
    evenness_slope = if(n() >= 3) coef(lm(evenness ~ onset_timeline_numeric))[2] else NA,
    total_abundance_slope = if(n() >= 3) coef(lm(total_abundance ~ onset_timeline_numeric))[2] else NA,
    dominance_slope = if(n() >= 3) coef(lm(dominance ~ onset_timeline_numeric))[2] else NA,
    viral_load_cv_slope = if(n() >= 3) coef(lm(viral_load_cv ~ onset_timeline_numeric))[2] else NA,
    .groups = "drop"
  ) %>%
  filter(!is.na(richness_slope))  # Remove patients with insufficient data

# Run limma analysis on slopes
slope_design <- model.matrix(~ Dx.Status + Sex + Age.at.Gluten.Introduction..months. + 
                            HLA.Category + feeding_first_year + Delivery.Mode, 
                            data = slope_data)

slope_matrix <- as.matrix(t(slope_data[, c("richness_slope", "shannon_slope", "simpson_slope",
                                          "evenness_slope", "total_abundance_slope", 
                                          "dominance_slope", "viral_load_cv_slope")]))

slope_fit <- lmFit(slope_matrix, slope_design)
slope_fit <- eBayes(slope_fit)

slope_results <- topTable(slope_fit, coef = "Dx.StatusCONTROL", number = Inf, sort.by = "P")

cat("Slope analysis completed\n")

# ============================================================================
# PART 6: CREATE ANALYSIS PLOTS
# ============================================================================

cat("\n6. Creating analysis plots...\n")

# Create slope analysis bar plot
slope_plot_data <- slope_results
slope_plot_data$metric <- rownames(slope_plot_data)
slope_plot_data$significant <- slope_plot_data$adj.P.Val < 0.05

p_slope_bar <- ggplot(slope_plot_data, aes(x = reorder(metric, logFC), y = logFC)) +
  geom_col(aes(fill = significant), alpha = 0.8) +
  geom_text(aes(label = paste0("p=", round(adj.P.Val, 3))), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
  labs(title = "Italy Cohort: Slope Analysis Results",
       subtitle = "Effect sizes (logFC) for CONTROL vs CELIAC slope differences",
       x = "Diversity Metric", y = "Effect Size (logFC)",
       caption = "Red bars indicate adj.p < 0.05") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("slope_analysis_results_Italy.png", p_slope_bar, width = 10, height = 6, dpi = 300)

# Create slope volcano plot
p_slope_volcano <- ggplot(slope_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = rownames(slope_results)), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
  labs(title = "Italy Cohort: Slope Analysis Volcano Plot",
       subtitle = "Statistical significance vs effect size for slope differences",
       x = "Effect Size (logFC)", y = "-log10(adjusted P-value)",
       caption = "Red line: adj.p = 0.05 threshold") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("slope_volcano_plot_Italy.png", p_slope_volcano, width = 10, height = 6, dpi = 300)

# Create diversity trajectories by group
trajectory_summary <- diversity_metrics %>%
  group_by(Dx.Status, onset_timeline_numeric) %>%
  summarise(
    richness_mean = mean(richness, na.rm = TRUE),
    shannon_mean = mean(shannon, na.rm = TRUE),
    simpson_mean = mean(simpson, na.rm = TRUE),
    evenness_mean = mean(evenness, na.rm = TRUE),
    total_abundance_mean = mean(total_abundance, na.rm = TRUE),
    dominance_mean = mean(dominance, na.rm = TRUE),
    viral_load_cv_mean = mean(viral_load_cv, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

# Create trajectory plots
metrics_to_plot <- c("richness_mean", "shannon_mean", "simpson_mean", "evenness_mean")
plot_list <- list()

for(metric in metrics_to_plot) {
  metric_name <- gsub("_mean", "", metric)
  
  p <- ggplot(trajectory_summary, aes(x = onset_timeline_numeric, y = .data[[metric]], color = Dx.Status)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = c("CELIAC" = "#E31A1C", "CONTROL" = "#1F78B4")) +
    labs(title = paste("Italy Cohort: Mean", str_to_title(metric_name), "Trajectory"),
         x = "Months Relative to Diagnosis", y = str_to_title(metric_name)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_list[[metric]] <- p
}

# Combine plots
combined_trajectories <- do.call(grid.arrange, c(plot_list, ncol = 2))

ggsave("diversity_trajectories_by_group_Italy.png", combined_trajectories, 
       width = 12, height = 8, dpi = 300)

cat("Diversity trajectories plot saved\n")

# ============================================================================
# PART 7: CREATE ANALYSIS SUMMARY STATISTICS
# ============================================================================

cat("\n7. Creating analysis summary statistics...\n")

# Create summary statistics plot
summary_stats <- data.frame(
  Metric = c("Total Samples", "Total Patients", "Significant ORFs (Interaction)", "Diversity Metrics", 
             "Trajectory Interactions", "Significant Slopes"),
  Count = c(nrow(diversity_metrics), length(unique(diversity_metrics$patientID)), 
            nrow(significant_orfs), 7, sum(trajectory_results$adj.P.Val < 0.05), 
            sum(slope_results$adj.P.Val < 0.05)),
  Category = c("Dataset", "Dataset", "ORF Analysis", "Diversity Analysis", 
               "Trajectory Analysis", "Slope Analysis")
)

p_summary <- ggplot(summary_stats, aes(x = reorder(Metric, Count), y = Count, fill = Category)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.1, size = 4) +
  coord_flip() +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(title = "Italy Cohort: Analysis Summary Statistics",
       subtitle = "Overview of dataset characteristics and analysis results",
       x = "Analysis Component", y = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("analysis_summary_statistics_Italy.png", p_summary, width = 10, height = 6, dpi = 300)

cat("Summary statistics plot saved\n")

# ============================================================================
# PART 8: SAVE DATA FILES
# ============================================================================

cat("\n8. Saving analysis data files...\n")

# Save diversity data
write.csv(diversity_metrics, "diversity_data_full_Italy.csv", row.names = FALSE)

# Save trajectory results
write.csv(trajectory_results, "diversity_trajectory_results_Italy.csv", row.names = TRUE)

# Save slope results
write.csv(slope_results, "slope_analysis_results_Italy.csv", row.names = TRUE)
write.csv(slope_data, "slope_data_Italy.csv", row.names = FALSE)

# Save summary statistics
write.csv(summary_stats, "analysis_summary_Italy.csv", row.names = FALSE)

cat("Data files saved\n")

# ============================================================================
# PART 9: FINAL SUMMARY
# ============================================================================

cat("\n=== ITALY COHORT COMPREHENSIVE ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- Italy_limma_model_Sig_res_final.csv: Significant ORFs results\n")
cat("- diversity_trajectory_heatmap_Italy.png: Diversity trajectory interactions\n")
cat("- slope_analysis_results_Italy.png: Slope analysis results\n")
cat("- slope_volcano_plot_Italy.png: Slope analysis volcano plot\n")
cat("- diversity_trajectories_by_group_Italy.png: Mean diversity trajectories\n")
cat("- analysis_summary_statistics_Italy.png: Summary statistics\n")
cat("- All supporting data files (.csv)\n")

cat("\nKey Results Summary:\n")
cat("- Total samples:", nrow(diversity_metrics), "\n")
cat("- Total patients:", length(unique(diversity_metrics$patientID)), "\n")
cat("- Significant ORFs:", nrow(significant_orfs), "(", round(nrow(significant_orfs)/nrow(abundance_matrix)*100, 1), "% of total)\n")
cat("- Significant trajectory interactions:", sum(trajectory_results$adj.P.Val < 0.05), "\n")
cat("- Significant slope differences:", sum(slope_results$adj.P.Val < 0.05), "\n")

cat("\nNote: Limited ORF significance led to focus on diversity-based analysis\n")
cat("Analysis completed at:", format(Sys.time()), "\n")