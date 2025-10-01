#!/usr/bin/env Rscript
# Create Annotated Cluster Heatmap - US Cohort
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Create heatmap with cluster annotations to show temporal patterns clearly

library(dplyr)
library(pheatmap)
library(RColorBrewer)

cat("=== CREATING ANNOTATED CLUSTER HEATMAP ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# Load the existing heatmap data
cat("Loading existing heatmap matrix data...\n")
heatmap_matrix_clean <- read.csv("temporal_heatmap_matrix_US_limma_fitted.csv", row.names = 1)
heatmap_matrix_clean <- as.matrix(heatmap_matrix_clean)

cat("Loaded heatmap matrix:", nrow(heatmap_matrix_clean), "ORFs x", ncol(heatmap_matrix_clean), "timepoints\n")

# Calculate variances and get top 200 variable ORFs (same as original)
orf_variances <- apply(heatmap_matrix_clean, 1, function(x) var(x, na.rm = TRUE))
orf_variances <- orf_variances[!is.na(orf_variances)]
n_top <- min(200, length(orf_variances))
top_variable_orfs <- names(sort(orf_variances, decreasing = TRUE)[1:n_top])
top_heatmap_matrix <- heatmap_matrix_clean[top_variable_orfs, ]

cat("Selected top", n_top, "most variable ORFs for cluster analysis\n")

# Set up colors (same as original)
data_range <- range(heatmap_matrix_clean, na.rm = TRUE)
max_abs <- max(abs(data_range), na.rm = TRUE)
color_limits <- c(-max_abs, max_abs)
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Perform hierarchical clustering to get the row order
cat("Performing hierarchical clustering...\n")
row_distances <- dist(top_heatmap_matrix)
row_clustering <- hclust(row_distances)
clustered_order <- row_clustering$order

# Get the clustered matrix (same order as the original heatmap)
clustered_matrix <- top_heatmap_matrix[clustered_order, ]

cat("Matrix clustered, analyzing temporal patterns...\n")

# Analyze temporal patterns to define clusters based on when transitions occur
# For each ORF, find when the major transition from red to blue occurs
transition_timepoints <- numeric(nrow(clustered_matrix))
early_values <- numeric(nrow(clustered_matrix))
late_values <- numeric(nrow(clustered_matrix))

for(i in 1:nrow(clustered_matrix)) {
  orf_profile <- clustered_matrix[i, ]
  valid_values <- orf_profile[!is.na(orf_profile)]
  
  if(length(valid_values) >= 3) {
    # Find the timepoint where the biggest change occurs
    differences <- abs(diff(valid_values))
    max_change_idx <- which.max(differences) + 1
    
    # Convert back to original timepoint index
    transition_timepoints[i] <- max_change_idx
    
    # Calculate early vs late average values
    n_vals <- length(valid_values)
    early_mean <- mean(valid_values[1:min(4, n_vals)])  # First ~4 timepoints
    late_mean <- mean(valid_values[max(1, n_vals-3):n_vals])  # Last ~4 timepoints
    
    early_values[i] <- early_mean
    late_values[i] <- late_mean
  }
}

# Define clusters based on transition timing and pattern
cluster_labels <- character(nrow(clustered_matrix))

for(i in 1:nrow(clustered_matrix)) {
  transition_time <- transition_timepoints[i]
  early_val <- early_values[i]
  late_val <- late_values[i]
  transition_magnitude <- late_val - early_val
  
  # Early Switchers: Transition occurs in first 1/3 of timeline (early timepoints)
  # and shows clear red-to-blue pattern (positive transition magnitude)
  if(transition_time <= 4 && transition_magnitude > 1) {
    cluster_labels[i] <- "Early Switchers"
  }
  # Late Responders: Transition occurs in last 1/3 of timeline (late timepoints)
  # or shows delayed/weak patterns
  else if(transition_time >= 10 || abs(transition_magnitude) < 1) {
    cluster_labels[i] <- "Late Responders"
  }
  # Gradual Progressors: Everything else - middle transitions with moderate changes
  else {
    cluster_labels[i] <- "Gradual Progressors"
  }
}

# Create annotation data frame
cluster_annotation <- data.frame(
  Temporal_Pattern = factor(cluster_labels, 
                           levels = c("Early Switchers", "Gradual Progressors", "Late Responders")),
  row.names = rownames(clustered_matrix)
)

# Define cluster colors
cluster_colors <- c(
  "Early Switchers" = "#E31A1C",      # Red
  "Gradual Progressors" = "#FF7F00",   # Orange  
  "Late Responders" = "#1F78B4"       # Blue
)

annotation_colors <- list(Temporal_Pattern = cluster_colors)

cat("Cluster assignment summary:\n")
cluster_summary <- table(cluster_labels)
print(cluster_summary)

# Save cluster assignments for reference
cluster_info <- data.frame(
  ORF = rownames(clustered_matrix),
  Cluster = cluster_labels,
  Transition_Timepoint = transition_timepoints,
  Early_Mean = round(early_values, 2),
  Late_Mean = round(late_values, 2),
  Transition_Magnitude = round(late_values - early_values, 2),
  stringsAsFactors = FALSE
)

write.csv(cluster_info, "temporal_cluster_assignments.csv", row.names = FALSE)

# Create the annotated heatmap
cat("Creating annotated cluster heatmap...\n")

png("top_variable_orfs_annotated_clusters_heatmap_US.png", 
    width = 2400, height = max(800, n_top * 15), res = 150)

pheatmap(clustered_matrix,
         cluster_rows = FALSE,  # Don't re-cluster, use existing order
         cluster_cols = FALSE, 
         scale = "none",
         color = colors, 
         breaks = seq(color_limits[1], color_limits[2], length.out = 101),
         annotation_row = cluster_annotation,
         annotation_colors = annotation_colors,
         main = paste("US Cohort: Top", n_top, "Variable ORFs with Temporal Pattern Clusters\n",
                      "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 18, 
         fontsize_row = max(8, min(12, 1000/n_top)), 
         fontsize_col = 20,
         cex_main = 1.5,
         show_rownames = FALSE,  # Too many to show individual names
         show_colnames = TRUE, 
         border_color = NA,
         annotation_names_row = TRUE,
         annotation_legend = TRUE)

dev.off()

cat("Annotated cluster heatmap created successfully!\n")

# Create a summary plot showing cluster characteristics
library(ggplot2)

# Prepare data for summary plot
cluster_summary_df <- data.frame(
  Cluster = names(cluster_summary),
  Count = as.numeric(cluster_summary),
  Color = cluster_colors[names(cluster_summary)]
)

# Summary bar plot
p_summary <- ggplot(cluster_summary_df, aes(x = reorder(Cluster, -Count), y = Count, fill = Cluster)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = cluster_colors) +
  geom_text(aes(label = Count), vjust = -0.3, size = 5) +
  labs(title = "Temporal Pattern Clusters - ORF Distribution",
       subtitle = paste("Total:", sum(cluster_summary), "ORFs classified into temporal programs"),
       x = "Temporal Cluster", y = "Number of ORFs") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))

ggsave("temporal_clusters_summary.png", p_summary, width = 10, height = 6, dpi = 300)

# Print final summary
cat("\n=== ANNOTATED CLUSTER HEATMAP COMPLETED ===\n")
cat("Generated files:\n")
cat("- top_variable_orfs_annotated_clusters_heatmap_US.png: Main annotated heatmap with cluster bars\n")
cat("- temporal_cluster_assignments.csv: Detailed cluster assignments and statistics\n")
cat("- temporal_clusters_summary.png: Summary plot of cluster distribution\n")

cat("\nCluster Characteristics:\n")
for(cluster in names(cluster_summary)) {
  cluster_orfs <- cluster_info[cluster_info$Cluster == cluster, ]
  avg_transition <- round(mean(cluster_orfs$Transition_Magnitude, na.rm = TRUE), 2)
  avg_timepoint <- round(mean(cluster_orfs$Transition_Timepoint, na.rm = TRUE), 1)
  
  cat("-", cluster, ":", cluster_summary[cluster], "ORFs\n")
  cat("  Average transition magnitude:", avg_transition, "\n")
  cat("  Average transition timepoint:", avg_timepoint, "\n")
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")