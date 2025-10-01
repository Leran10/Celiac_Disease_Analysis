#!/usr/bin/env Rscript
# Regenerate Pattern Clustering Analysis with Updated 2,009 ORFs
# Author: Claude AI
# Date: 2025-08-05

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== REGENERATING PATTERN CLUSTERING ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD EXISTING LIMMA FITTED VALUES
# ============================================================================

cat("1. Loading existing limma fitted values...\n")

# Load the saved model results from our previous analysis
model_results <- readRDS("limma_model_results_fitted_values.rds")
cat("Loaded model results for", length(model_results), "ORFs\n")

# Load the heatmap matrix we created earlier
heatmap_matrix <- read.csv("temporal_heatmap_matrix_limma_fitted.csv", row.names = 1)
cat("Loaded heatmap matrix with", nrow(heatmap_matrix), "ORFs and", ncol(heatmap_matrix), "timepoints\n")

# ============================================================================
# PART 2: CALCULATE PATTERN FEATURES FOR CLUSTERING
# ============================================================================

cat("\n2. Calculating pattern features for clustering...\n")

# Function to calculate pattern features for each ORF
calculate_pattern_features <- function(orf_profile) {
  # Remove NAs for calculations
  valid_values <- orf_profile[!is.na(orf_profile)]
  
  if(length(valid_values) < 3) {
    return(list(slope = NA, variance = NA, early_late_diff = NA, 
                peak_timing = NA, monotonicity = NA, pattern_type = "insufficient_data"))
  }
  
  # 1. Slope (linear trend)
  timepoints <- 1:length(valid_values)
  slope <- coef(lm(valid_values ~ timepoints))[2]
  
  # 2. Variance (temporal variability)
  variance <- var(valid_values)
  
  # 3. Early-late difference
  n <- length(valid_values)
  early_mean <- mean(valid_values[1:min(3, n)])
  late_mean <- mean(valid_values[max(1, n-2):n])
  early_late_diff <- late_mean - early_mean
  
  # 4. Peak timing (normalized position of maximum absolute value)
  peak_timing <- which.max(abs(valid_values)) / length(valid_values)
  
  # 5. Monotonicity (how consistently increasing/decreasing)
  diffs <- diff(valid_values)
  monotonicity <- abs(sum(sign(diffs))) / length(diffs)
  
  # 6. Pattern type classification
  if(abs(slope) < 0.1 && variance < 1) {
    pattern_type <- "stable"
  } else if(slope > 0.1) {
    pattern_type <- "increasing"
  } else if(slope < -0.1) {
    pattern_type <- "decreasing"
  } else {
    pattern_type <- "complex"
  }
  
  return(list(slope = slope, variance = variance, early_late_diff = early_late_diff,
              peak_timing = peak_timing, monotonicity = monotonicity, pattern_type = pattern_type))
}

# Calculate features for all ORFs
pattern_features <- list()
for(i in 1:nrow(heatmap_matrix)) {
  orf_name <- rownames(heatmap_matrix)[i]
  orf_profile <- as.numeric(heatmap_matrix[i, ])
  pattern_features[[orf_name]] <- calculate_pattern_features(orf_profile)
  
  if(i %% 200 == 0) {
    cat("Calculated features for", i, "of", nrow(heatmap_matrix), "ORFs\n")
  }
}

# Convert to data frame
features_df <- data.frame(
  orf_name = names(pattern_features),
  slope = sapply(pattern_features, function(x) x$slope),
  variance = sapply(pattern_features, function(x) x$variance),
  early_late_diff = sapply(pattern_features, function(x) x$early_late_diff),
  peak_timing = sapply(pattern_features, function(x) x$peak_timing),
  monotonicity = sapply(pattern_features, function(x) x$monotonicity),
  pattern_type = sapply(pattern_features, function(x) x$pattern_type),
  stringsAsFactors = FALSE
)

# Remove ORFs with insufficient data
features_clean <- features_df[!is.na(features_df$slope), ]
cat("Pattern features calculated for", nrow(features_clean), "ORFs\n")

# ============================================================================
# PART 3: PERFORM CLUSTERING
# ============================================================================

cat("\n3. Performing pattern-based clustering...\n")

# Prepare feature matrix for clustering (standardize features)
feature_matrix <- features_clean[, c("slope", "variance", "early_late_diff", "peak_timing", "monotonicity")]
feature_matrix_scaled <- scale(feature_matrix)

# Remove any rows with NAs after scaling
complete_rows <- complete.cases(feature_matrix_scaled)
feature_matrix_clean <- feature_matrix_scaled[complete_rows, ]
features_for_clustering <- features_clean[complete_rows, ]

cat("Using", nrow(feature_matrix_clean), "ORFs for clustering\n")

# Determine optimal number of clusters using elbow method
wss <- sapply(1:10, function(k) {
  set.seed(123)
  kmeans(feature_matrix_clean, centers = k, nstart = 10)$tot.withinss
})

# Find elbow point (simple method)
optimal_k <- which.min(diff(diff(wss))) + 1
if(optimal_k < 3) optimal_k <- 4  # Ensure minimum clusters

cat("Optimal number of clusters:", optimal_k, "\n")

# Perform K-means clustering
set.seed(123)
kmeans_result <- kmeans(feature_matrix_clean, centers = optimal_k, nstart = 25)

# Add cluster assignments to features
features_for_clustering$cluster <- kmeans_result$cluster

# ============================================================================
# PART 4: CREATE CLUSTERED HEATMAP
# ============================================================================

cat("\n4. Creating pattern-based clustered heatmap...\n")

# Order ORFs by cluster and then by pattern characteristics within cluster
clustered_orfs <- features_for_clustering %>%
  arrange(cluster, pattern_type, desc(abs(slope))) %>%
  pull(orf_name)

# Filter heatmap matrix to clustered ORFs
heatmap_clustered <- heatmap_matrix[clustered_orfs, ]

# Create cluster annotation for rows
cluster_annotation <- data.frame(
  Cluster = paste0("Cluster_", features_for_clustering$cluster[match(clustered_orfs, features_for_clustering$orf_name)]),
  Pattern_Type = features_for_clustering$pattern_type[match(clustered_orfs, features_for_clustering$orf_name)]
)
rownames(cluster_annotation) <- clustered_orfs

# Create color palettes for annotations
cluster_colors <- rainbow(optimal_k, alpha = 0.8)
names(cluster_colors) <- paste0("Cluster_", 1:optimal_k)

pattern_colors <- c("increasing" = "#E31A1C", "decreasing" = "#1F78B4", 
                   "stable" = "#33A02C", "complex" = "#FF7F00")

annotation_colors <- list(
  Cluster = cluster_colors,
  Pattern_Type = pattern_colors
)

# Set color limits for heatmap
data_range <- range(heatmap_clustered, na.rm = TRUE)
max_abs <- max(abs(data_range), na.rm = TRUE)
color_limits <- c(-max_abs, max_abs)

# Create color palette
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Create the clustered heatmap
png("temporal_pattern_clusters_updated.png", width = 1400, height = 1200, res = 150)

pheatmap(heatmap_clustered,
         cluster_rows = FALSE,  # Don't cluster - use our pattern-based ordering
         cluster_cols = FALSE,
         scale = "none",
         color = colors,
         breaks = seq(color_limits[1], color_limits[2], length.out = 101),
         annotation_row = cluster_annotation,
         annotation_colors = annotation_colors,
         main = paste("Pattern-Based Clustering:", nrow(heatmap_clustered), "ORFs in", optimal_k, "Temporal Programs\n",
                      "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 8,
         fontsize_row = 4,
         fontsize_col = 10,
         show_rownames = FALSE,  # Too many ORFs to show names
         show_colnames = TRUE,
         border_color = NA,
         annotation_legend = TRUE)

dev.off()

cat("Pattern-based clustered heatmap saved\n")

# ============================================================================
# PART 5: CREATE SUMMARY STATISTICS
# ============================================================================

cat("\n5. Creating summary statistics...\n")

# Cluster summary
cluster_summary <- features_for_clustering %>%
  group_by(cluster, pattern_type) %>%
  summarise(
    n_orfs = n(),
    mean_slope = mean(slope, na.rm = TRUE),
    mean_variance = mean(variance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cluster, pattern_type)

write.csv(cluster_summary, "pattern_clustering_summary_updated.csv", row.names = FALSE)

# Save detailed features
write.csv(features_for_clustering, "pattern_features_updated.csv", row.names = FALSE)

cat("Summary statistics saved\n")

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== PATTERN CLUSTERING ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- temporal_pattern_clusters_updated.png: Pattern-based clustered heatmap\n")
cat("- pattern_clustering_summary_updated.csv: Cluster summary statistics\n")
cat("- pattern_features_updated.csv: Detailed pattern features for all ORFs\n")

cat("\nCluster Summary:\n")
print(cluster_summary)

cat("\nTotal ORFs analyzed:", nrow(heatmap_matrix), "\n")
cat("ORFs successfully clustered:", nrow(features_for_clustering), "\n")
cat("Number of clusters identified:", optimal_k, "\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")