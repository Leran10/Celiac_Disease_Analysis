#!/usr/bin/env Rscript
# Linear Trend Temporal Heatmap Analysis - US Cohort
# Author: Claude AI
# Date: 2025-08-11
# Purpose: Create heatmap showing temporal progression of ORFs with significant linear trends

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== LINEAR TREND TEMPORAL HEATMAP ANALYSIS - US COHORT ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND FIT LIMMA MODEL
# ============================================================================

cat("1. Loading US cohort data and extracting linear trend results...\n")

# Load abundance data
abundance_data_raw <- read.csv("US.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names (remove X prefix if present)
sample_names <- colnames(abundance_data)
sample_names_clean <- gsub("^X", "", sample_names)
colnames(abundance_data) <- sample_names_clean

# Load metadata
metadata_raw <- read.csv("US.metadata.clean.csv")
rownames(metadata_raw) <- metadata_raw$X
metadata <- metadata_raw

# Match samples
matching_samples <- intersect(colnames(abundance_data), rownames(metadata))
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")

# Set reference levels (CELIAC as reference)
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))

# Create design matrix with NUMERICAL timeline (for linear trends)
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_numeric + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

cat("Design matrix dimensions:", nrow(linear.model.design), "x", ncol(linear.model.design), "\n")

# Convert abundance data to matrix for limma
abundance_matrix <- as.matrix(abundance_filtered)

# ============================================================================
# PART 2: EXTRACT LINEAR TREND RESULTS
# ============================================================================

cat("\n2. Fitting limma model and extracting onset_timeline_numeric results...\n")

# Fit limma model with voomLmFit and patient blocking
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get results for onset_timeline_numeric coefficient (linear trends in CELIAC group)
linear_trend_results <- topTable(linear.model.fit, coef = "onset_timeline_numeric", 
                                 number = Inf, sort.by = "P")

# Get significant linear trend ORFs
significant_trend_orfs <- linear_trend_results[linear_trend_results$adj.P.Val < 0.05, ]

cat("Total ORFs analyzed:", nrow(linear_trend_results), "\n")
cat("ORFs with significant linear trends:", nrow(significant_trend_orfs), "\n")
cat("Percentage with significant trends:", round(nrow(significant_trend_orfs)/nrow(linear_trend_results)*100, 1), "%\n")

# Print trend direction summary
positive_trends <- sum(significant_trend_orfs$logFC > 0)
negative_trends <- sum(significant_trend_orfs$logFC < 0)
cat("Positive trends (increasing toward diagnosis):", positive_trends, "\n")
cat("Negative trends (decreasing toward diagnosis):", negative_trends, "\n")

# Save linear trend results
write.csv(significant_trend_orfs, "significant_linear_trend_orfs.csv", row.names = TRUE)

# ============================================================================
# PART 3: CREATE LINEAR TREND TEMPORAL MATRIX
# ============================================================================

cat("\n3. Creating linear trend temporal matrix...\n")

# Define timepoints from the data
timepoints <- sort(unique(metadata_filtered$onset_timeline_numeric))
cat("Timeline timepoints:", paste(timepoints, collapse = ", "), "\n")

# Get intercept values for significant trending ORFs
intercept_results <- topTable(linear.model.fit, coef = "(Intercept)", 
                              number = Inf, sort.by = "none")

# Create matrix to store predicted values
trend_matrix <- matrix(NA, 
                       nrow = nrow(significant_trend_orfs), 
                       ncol = length(timepoints),
                       dimnames = list(rownames(significant_trend_orfs), 
                                      ifelse(timepoints == 0, "T0", paste0("T", timepoints))))

cat("Creating linear trend predictions for", nrow(significant_trend_orfs), "ORFs...\n")

# Calculate predicted values for each ORF at each timepoint using linear equation
for(i in 1:nrow(significant_trend_orfs)) {
  orf_id <- rownames(significant_trend_orfs)[i]
  
  # Get slope (onset_timeline_numeric coefficient) and intercept
  slope <- significant_trend_orfs[orf_id, "logFC"]
  intercept <- intercept_results[orf_id, "logFC"]
  
  # Calculate predicted value at each timepoint: predicted = intercept + slope * timepoint
  for(j in 1:length(timepoints)) {
    timepoint <- timepoints[j]
    predicted_value <- intercept + slope * timepoint
    trend_matrix[i, j] <- predicted_value
  }
  
  if(i %% 500 == 0) {
    cat("Processed", i, "of", nrow(significant_trend_orfs), "ORFs\n")
  }
}

cat("Linear trend matrix created:", nrow(trend_matrix), "ORFs x", ncol(trend_matrix), "timepoints\n")
cat("Data range:", round(range(trend_matrix, na.rm = TRUE), 2), "\n")

# Save the trend matrix
write.csv(trend_matrix, "linear_trend_temporal_matrix.csv")

# ============================================================================
# PART 4: CREATE TREND HEATMAPS
# ============================================================================

cat("\n4. Creating linear trend temporal heatmaps...\n")

if(nrow(trend_matrix) > 0) {
  # Set up color scheme
  data_range <- range(trend_matrix, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # Heatmap 1: All significant linear trend ORFs
  if(nrow(trend_matrix) <= 1000) {  # Only if manageable number
    png("all_significant_linear_trend_orfs_heatmap.png", 
        width = 2200, height = max(1200, nrow(trend_matrix) * 8), res = 150)
    pheatmap(trend_matrix,
             cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
             color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
             main = paste("US Cohort: All", nrow(trend_matrix), "ORFs with Significant Linear Trends (CELIAC Group)\n",
                          "Red = Higher Abundance, Blue = Lower Abundance"),
             fontsize = 14, fontsize_row = 6, fontsize_col = 16,
             cex_main = 1.2,
             show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
    dev.off()
    cat("All significant trend ORFs heatmap created\n")
  }
  
  # Heatmap 2: Top variable linear trend ORFs
  trend_variances <- apply(trend_matrix, 1, function(x) var(x, na.rm = TRUE))
  trend_variances <- trend_variances[!is.na(trend_variances)]
  n_top <- min(200, length(trend_variances))
  top_variable_trend_orfs <- names(sort(trend_variances, decreasing = TRUE)[1:n_top])
  top_trend_matrix <- trend_matrix[top_variable_trend_orfs, ]
  
  png("top_variable_linear_trend_orfs_heatmap.png", 
      width = 2200, height = max(800, n_top * 15), res = 150)
  pheatmap(top_trend_matrix,
           cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
           color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("US Cohort: Top", n_top, "Most Variable Linear Trend ORFs (CELIAC Group)\n",
                        "Red = Higher Abundance, Blue = Lower Abundance"),
           fontsize = 18, fontsize_row = max(10, min(14, 1200/n_top)), fontsize_col = 20,
           cex_main = 1.5,
           show_rownames = TRUE, show_colnames = TRUE, border_color = NA)
  dev.off()
  cat("Top variable trend ORFs heatmap created\n")
  
  # Heatmap 3: Separate positive and negative trends
  positive_trend_orfs <- rownames(significant_trend_orfs)[significant_trend_orfs$logFC > 0]
  negative_trend_orfs <- rownames(significant_trend_orfs)[significant_trend_orfs$logFC < 0]
  
  if(length(positive_trend_orfs) > 0 && length(positive_trend_orfs) <= 500) {
    positive_matrix <- trend_matrix[positive_trend_orfs, ]
    png("positive_linear_trend_orfs_heatmap.png", 
        width = 2200, height = max(800, length(positive_trend_orfs) * 10), res = 150)
    pheatmap(positive_matrix,
             cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
             color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
             main = paste("US Cohort:", length(positive_trend_orfs), "ORFs with Positive Linear Trends (CELIAC Group)\n",
                          "Increasing toward diagnosis - Red = Higher, Blue = Lower"),
             fontsize = 16, fontsize_row = 8, fontsize_col = 18,
             cex_main = 1.3,
             show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
    dev.off()
    cat("Positive trend ORFs heatmap created\n")
  }
  
  if(length(negative_trend_orfs) > 0 && length(negative_trend_orfs) <= 500) {
    negative_matrix <- trend_matrix[negative_trend_orfs, ]
    png("negative_linear_trend_orfs_heatmap.png", 
        width = 2200, height = max(800, length(negative_trend_orfs) * 10), res = 150)
    pheatmap(negative_matrix,
             cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
             color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
             main = paste("US Cohort:", length(negative_trend_orfs), "ORFs with Negative Linear Trends (CELIAC Group)\n",
                          "Decreasing toward diagnosis - Red = Higher, Blue = Lower"),
             fontsize = 16, fontsize_row = 8, fontsize_col = 18,
             cex_main = 1.3,
             show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
    dev.off()
    cat("Negative trend ORFs heatmap created\n")
  }
}

# ============================================================================
# PART 5: COMPARISON WITH FITTED VALUES APPROACH
# ============================================================================

cat("\n5. Creating comparison summary...\n")

# Load previous interaction results for comparison
if(file.exists("US_limma_model_Sig_res_final.csv")) {
  interaction_results <- read.csv("US_limma_model_Sig_res_final.csv", row.names = 1)
  interaction_orfs <- rownames(interaction_results)
  
  # Compare overlaps
  trend_orfs <- rownames(significant_trend_orfs)
  overlap_orfs <- intersect(trend_orfs, interaction_orfs)
  
  comparison_summary <- data.frame(
    Analysis_Type = c("Linear Trends (onset_timeline_numeric)", 
                     "Interaction Effects (Dx.Status:onset_timeline_numeric)",
                     "Overlapping ORFs"),
    Number_of_ORFs = c(length(trend_orfs), length(interaction_orfs), length(overlap_orfs)),
    Percentage = c(round(length(trend_orfs)/nrow(abundance_matrix)*100, 1),
                   round(length(interaction_orfs)/nrow(abundance_matrix)*100, 1),
                   round(length(overlap_orfs)/min(length(trend_orfs), length(interaction_orfs))*100, 1)),
    stringsAsFactors = FALSE
  )
  
  cat("\n=== COMPARISON WITH INTERACTION ANALYSIS ===\n")
  print(comparison_summary)
  
  write.csv(comparison_summary, "linear_trend_vs_interaction_comparison.csv", row.names = FALSE)
}

# Create summary statistics
summary_stats <- data.frame(
  Metric = c("Total ORFs Analyzed", "Significant Linear Trend ORFs", "Positive Trends", "Negative Trends",
             "Percentage Significant", "Data Range Min", "Data Range Max", "Top Variable ORFs"),
  Value = c(nrow(linear_trend_results), nrow(significant_trend_orfs), positive_trends, negative_trends,
            paste0(round(nrow(significant_trend_orfs)/nrow(linear_trend_results)*100, 1), "%"),
            round(min(trend_matrix, na.rm = TRUE), 2), round(max(trend_matrix, na.rm = TRUE), 2), n_top),
  stringsAsFactors = FALSE
)

write.csv(summary_stats, "linear_trend_analysis_summary.csv", row.names = FALSE)

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== LINEAR TREND TEMPORAL HEATMAP ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
if(exists("n_top")) {
  cat("- top_variable_linear_trend_orfs_heatmap.png: Top", n_top, "most variable trend ORFs\n")
}
if(nrow(trend_matrix) <= 1000) {
  cat("- all_significant_linear_trend_orfs_heatmap.png: All significant trend ORFs\n")
}
if(length(positive_trend_orfs) > 0 && length(positive_trend_orfs) <= 500) {
  cat("- positive_linear_trend_orfs_heatmap.png: ORFs increasing toward diagnosis\n")
}
if(length(negative_trend_orfs) > 0 && length(negative_trend_orfs) <= 500) {
  cat("- negative_linear_trend_orfs_heatmap.png: ORFs decreasing toward diagnosis\n")
}
cat("- linear_trend_temporal_matrix.csv: Complete trend matrix data\n")
cat("- significant_linear_trend_orfs.csv: Significant trend ORFs with statistics\n")
cat("- linear_trend_analysis_summary.csv: Analysis summary statistics\n")
cat("- linear_trend_vs_interaction_comparison.csv: Comparison with interaction analysis\n")

cat("\nKey Results Summary:\n")
cat("- Total ORFs analyzed:", nrow(linear_trend_results), "\n")
cat("- Significant linear trend ORFs:", nrow(significant_trend_orfs), 
    "(", round(nrow(significant_trend_orfs)/nrow(linear_trend_results)*100, 1), "% of total)\n")
cat("- Positive trends (increasing toward diagnosis):", positive_trends, "\n")
cat("- Negative trends (decreasing toward diagnosis):", negative_trends, "\n")
cat("- Data range (predicted values):", round(min(trend_matrix, na.rm = TRUE), 2), 
    "to", round(max(trend_matrix, na.rm = TRUE), 2), "\n")

if(positive_trends > negative_trends) {
  cat("- Predominant pattern: ORFs INCREASE toward diagnosis in CELIAC group\n")
} else if(negative_trends > positive_trends) {
  cat("- Predominant pattern: ORFs DECREASE toward diagnosis in CELIAC group\n")
} else {
  cat("- Pattern: Balanced increase and decrease trends\n")
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")