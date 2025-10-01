#!/usr/bin/env Rscript
# Investigate the contradiction between limma results and temporal patterns
# Why do some ORFs show red (higher in CELIAC) when limma says positive logFC (higher in CONTROL)?

library(dplyr)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/compositonal_analysis")

cat("=== INVESTIGATING LIMMA vs TEMPORAL CONTRADICTION ===\n")

# Load limma results
limma_results <- read.csv("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/total_limma_results.csv", row.names = 1)

# Load temporal matrix
temporal_matrix <- read.csv("confounder_adjusted_temporal_matrix.csv", row.names = 1)

# Focus on the specific ORF that shows red (virus_comp_315_cycle_1_1)
example_orf <- "virus_comp_315_cycle_1_1"

cat("Investigating ORF:", example_orf, "\n")
cat("Limma logFC:", limma_results[example_orf, "logFC"], "\n")
cat("Limma interpretation: Positive logFC = Higher in CONTROL, Lower in CELIAC\n")
cat("Temporal differences (CELIAC - CONTROL):\n")
print(temporal_matrix[example_orf, ])

cat("\n=== THE ISSUE ===\n")
cat("Temporal matrix shows CELIAC - CONTROL differences\n")
cat("Positive values = Higher in CELIAC\n")
cat("But limma logFC is CONTROL - CELIAC comparison!\n")

cat("\n=== SOLUTION ===\n")
cat("We need to flip the sign of temporal differences to match limma interpretation\n")
cat("OR clearly specify what the temporal heatmap represents\n")

# Check the mean temporal difference vs limma logFC for validation
temporal_means <- rowMeans(temporal_matrix, na.rm = TRUE)
matched_orfs <- intersect(rownames(temporal_matrix), rownames(limma_results))

comparison_data <- data.frame(
  orf_id = matched_orfs,
  limma_logfc = limma_results[matched_orfs, "logFC"],
  temporal_mean = temporal_means[matched_orfs],
  expected_sign = ifelse(limma_results[matched_orfs, "logFC"] > 0, "negative", "positive")
)

cat("\nComparison for first 10 ORFs:\n")
print(head(comparison_data, 10))

# Check correlation
correlation <- cor(comparison_data$limma_logfc, comparison_data$temporal_mean, use = "complete.obs")
cat("\nCorrelation between limma logFC and temporal mean:", correlation, "\n")
cat("Expected: NEGATIVE correlation (opposite signs)\n")

# Count sign agreements
sign_agreement <- sum(sign(comparison_data$limma_logfc) == -sign(comparison_data$temporal_mean), na.rm = TRUE)
total_valid <- sum(!is.na(comparison_data$limma_logfc) & !is.na(comparison_data$temporal_mean))

cat("Sign agreement (limma vs -temporal):", sign_agreement, "/", total_valid, 
    "(", round(100*sign_agreement/total_valid, 1), "%)\n")

cat("\n=== RECOMMENDATION ===\n")
cat("The temporal heatmap should either:\n")
cat("1. Flip the sign: show CONTROL - CELIAC differences\n")
cat("2. Update the legend: clarify it shows CELIAC - CONTROL differences\n")
cat("3. Add validation plot showing the relationship\n")