#!/usr/bin/env Rscript
# Italy Intersection Top 200 Heatmap
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Create top 200 most variable ORFs heatmap from intersection set (t0-24 AND t0-30)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== ITALY INTERSECTION TOP 200 HEATMAP ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD INTERSECTION MATRIX DATA
# ============================================================================

cat("1. Loading intersection matrix data...\n")

# Load the intersection heatmap matrix
intersection_matrix <- read.csv("italy_intersection_heatmap_matrix.csv", row.names = 1)

cat("Loaded intersection matrix:", nrow(intersection_matrix), "ORFs x", ncol(intersection_matrix), "timepoints\n")
cat("Data range:", round(range(intersection_matrix, na.rm = TRUE), 2), "\n")

# ============================================================================
# PART 2: IDENTIFY TOP 200 MOST VARIABLE ORFS
# ============================================================================

cat("2. Identifying top 200 most variable ORFs from intersection set...\n")

# Calculate variance for each ORF across all timepoints
orf_variances <- apply(intersection_matrix, 1, function(x) var(x, na.rm = TRUE))
orf_variances_clean <- orf_variances[!is.na(orf_variances)]

cat("Total ORFs with variance data:", length(orf_variances_clean), "\n")
cat("Variance range:", round(range(orf_variances_clean), 3), "\n")

# Get top 200 most variable ORFs
n_top <- min(200, length(orf_variances_clean))
top_variable_orfs <- names(sort(orf_variances_clean, decreasing = TRUE)[1:n_top])
top_intersection_matrix <- intersection_matrix[top_variable_orfs, ]

cat("Selected top", n_top, "most variable ORFs from intersection set\n")
cat("Top ORF variances (first 10):", round(sort(orf_variances_clean, decreasing = TRUE)[1:10], 3), "\n")

# ============================================================================
# PART 3: CREATE TOP 200 INTERSECTION HEATMAP
# ============================================================================

cat("3. Creating top 200 intersection heatmap...\n")

# Set up color scheme
data_range <- range(top_intersection_matrix, na.rm = TRUE)
max_abs <- max(abs(data_range), na.rm = TRUE)
color_limits <- c(-max_abs, max_abs)
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

cat("Color scale range:", round(color_limits, 2), "\n")

# Create the top 200 intersection heatmap
png("italy_intersection_heatmap_top200_t0-24_AND_t0-30.png", 
    width = 2400, height = max(800, n_top * 15), res = 150)

pheatmap(top_intersection_matrix,
         cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
         color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
         main = paste("Italy Cohort: Top", n_top, "Most Variable Intersection ORFs (t0-24 AND t0-30)\n",
                     "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 18, 
         fontsize_row = max(8, min(12, 1000/n_top)), 
         fontsize_col = 20,
         cex_main = 1.2,
         show_rownames = TRUE,  # Show ORF names for top 200
         show_colnames = TRUE, 
         border_color = NA)

dev.off()

cat("Top", n_top, "intersection ORFs heatmap created successfully!\n")

# ============================================================================
# PART 4: SAVE TOP 200 DATA AND STATISTICS
# ============================================================================

cat("4. Saving top 200 intersection data and statistics...\n")

# Save the top 200 intersection matrix
write.csv(top_intersection_matrix, "italy_intersection_heatmap_top200_matrix.csv")

# Save list of top 200 ORFs with their variance values
top200_info <- data.frame(
  ORF = top_variable_orfs,
  Variance = orf_variances_clean[top_variable_orfs],
  Rank = 1:n_top,
  stringsAsFactors = FALSE
)
write.csv(top200_info, "italy_intersection_top200_orfs_list.csv", row.names = FALSE)

# Create summary statistics for top 200
top200_summary <- data.frame(
  Timepoint = colnames(top_intersection_matrix),
  Valid_ORFs = colSums(!is.na(top_intersection_matrix)),
  Mean_Difference = round(colMeans(top_intersection_matrix, na.rm = TRUE), 3),
  SD_Difference = round(apply(top_intersection_matrix, 2, sd, na.rm = TRUE), 3),
  Range_Min = round(apply(top_intersection_matrix, 2, min, na.rm = TRUE), 3),
  Range_Max = round(apply(top_intersection_matrix, 2, max, na.rm = TRUE), 3),
  Analysis = "Top200_Intersection",
  stringsAsFactors = FALSE
)

write.csv(top200_summary, "italy_intersection_top200_summary.csv", row.names = FALSE)

# ============================================================================
# PART 5: COMPARISON WITH FULL INTERSECTION
# ============================================================================

cat("5. Comparing top 200 with full intersection set...\n")

# Load full intersection summary for comparison
full_intersection_summary <- read.csv("italy_intersection_heatmap_summary.csv")

# Create comparison table
comparison_table <- data.frame(
  Timepoint = top200_summary$Timepoint,
  Full_Intersection_Mean = full_intersection_summary$Mean_Difference,
  Top200_Intersection_Mean = top200_summary$Mean_Difference,
  Mean_Difference = round(top200_summary$Mean_Difference - full_intersection_summary$Mean_Difference, 3),
  Full_Intersection_SD = full_intersection_summary$SD_Difference,
  Top200_Intersection_SD = top200_summary$SD_Difference,
  stringsAsFactors = FALSE
)

write.csv(comparison_table, "italy_intersection_full_vs_top200_comparison.csv", row.names = FALSE)

cat("Comparison between full intersection and top 200:\n")
print(comparison_table)

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== ITALY INTERSECTION TOP 200 HEATMAP COMPLETED ===\n")
cat("Generated files:\n")
cat("- italy_intersection_heatmap_top200_t0-24_AND_t0-30.png: Top 200 intersection heatmap\n")
cat("- italy_intersection_heatmap_top200_matrix.csv: Top 200 intersection matrix data\n")
cat("- italy_intersection_top200_orfs_list.csv: List of top 200 ORFs with variance values\n")
cat("- italy_intersection_top200_summary.csv: Summary statistics for top 200\n")
cat("- italy_intersection_full_vs_top200_comparison.csv: Comparison with full intersection\n")

cat("\nKey statistics:\n")
cat("- ORFs selected:", n_top, "out of", nrow(intersection_matrix), "(top", round(n_top/nrow(intersection_matrix)*100, 1), "%)\n")
cat("- Variance range of selected ORFs:", round(range(orf_variances_clean[top_variable_orfs]), 3), "\n")
cat("- Data range in heatmap:", round(range(top_intersection_matrix, na.rm = TRUE), 2), "\n")

cat("\nTop 200 intersection summary by timepoint:\n")
print(top200_summary)

# Identify strongest patterns in top 200
strongest_positive <- top200_summary$Timepoint[which.max(top200_summary$Mean_Difference)]
strongest_negative <- top200_summary$Timepoint[which.min(top200_summary$Mean_Difference)]

cat("\n- Strongest CONTROL dominance in top 200:", strongest_positive, "(mean =", top200_summary$Mean_Difference[top200_summary$Timepoint == strongest_positive], ")\n")
cat("- Strongest CELIAC dominance in top 200:", strongest_negative, "(mean =", top200_summary$Mean_Difference[top200_summary$Timepoint == strongest_negative], ")\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")