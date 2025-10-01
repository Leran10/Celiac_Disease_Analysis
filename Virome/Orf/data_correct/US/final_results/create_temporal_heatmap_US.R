#!/usr/bin/env Rscript
# Temporal Heatmap Analysis - US Cohort Only
# Author: Claude AI
# Date: 2025-07-29

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/US/final_results")

cat("=== TEMPORAL HEATMAP ANALYSIS - US COHORT ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading US cohort data...\n")

# Load limma results
limma_results <- read.csv("../US.limma_model_res.csv")
cat("Loaded limma results with", nrow(limma_results), "ORFs\n")

# Load ORF abundance data
orf_data_raw <- read.csv("../US.orf.abundance.clean.csv")
rownames(orf_data_raw) <- orf_data_raw$X
orf_data <- orf_data_raw[, -1]

# Fix sample names (remove X prefix if present)
orf_samples <- colnames(orf_data)
orf_samples_clean <- gsub("^X", "", orf_samples)
colnames(orf_data) <- orf_samples_clean

cat("Loaded ORF abundance data with", nrow(orf_data), "ORFs and", ncol(orf_data), "samples\n")

# Load metadata
metadata <- read.csv("../US.metadata.clean.csv")
cat("Loaded metadata with", nrow(metadata), "samples\n")

# Check sample matching
metadata_samples <- metadata$X
common_samples <- intersect(orf_samples_clean, metadata_samples)
cat("Sample matching: Found", length(common_samples), "common samples\n")

if(length(common_samples) == 0) {
  stop("No common samples found between abundance data and metadata!")
}

# Filter data to common samples
orf_data_filtered <- orf_data[, common_samples]
metadata_filtered <- metadata[metadata$X %in% common_samples, ]

cat("Filtered to", ncol(orf_data_filtered), "samples and", nrow(metadata_filtered), "metadata rows\n")

# ============================================================================
# PART 2: FILTER SIGNIFICANT ORFS
# ============================================================================

cat("\n2. Filtering significant ORFs (adj.p < 0.01)...\n")

# Filter significant ORFs
significant_orfs <- limma_results[limma_results$adj.P.Val < 0.01, ]
cat("Found", nrow(significant_orfs), "significant ORFs (", 
    round(nrow(significant_orfs)/nrow(limma_results)*100, 1), "% of total)\n")

if(nrow(significant_orfs) == 0) {
  stop("No significant ORFs found at adj.p < 0.01 threshold!")
}

# Get ORF names that match between limma results and abundance data
significant_orf_names <- significant_orfs$X
available_orfs <- rownames(orf_data_filtered)
matching_orfs <- intersect(significant_orf_names, available_orfs)

cat("ORF matching:", length(matching_orfs), "significant ORFs found in abundance data\n")

if(length(matching_orfs) == 0) {
  stop("No significant ORFs found in abundance data!")
}

# Filter ORF data to significant ORFs
significant_orf_data <- orf_data_filtered[matching_orfs, ]

# ============================================================================
# PART 3: CREATE TEMPORAL HEATMAP DATA
# ============================================================================

cat("\n3. Creating temporal heatmap data...\n")

# Create timepoint mapping
metadata_filtered$timepoint_numeric <- as.numeric(metadata_filtered$onset_timeline_numeric)

# Get unique timepoints and sort them
timepoints <- sort(unique(metadata_filtered$timepoint_numeric), decreasing = FALSE)
cat("Found", length(timepoints), "unique timepoints:", paste(timepoints, collapse = ", "), "\n")

# Initialize heatmap matrix
heatmap_matrix <- matrix(NA, nrow = length(matching_orfs), ncol = length(timepoints))
rownames(heatmap_matrix) <- matching_orfs
colnames(heatmap_matrix) <- paste0("T", timepoints)

# Calculate log2 differences (CONTROL - CELIAC) for each ORF at each timepoint
for(i in 1:length(matching_orfs)) {
  orf_name <- matching_orfs[i]
  
  for(j in 1:length(timepoints)) {
    tp <- timepoints[j]
    
    # Get samples for this timepoint
    tp_samples <- metadata_filtered$X[metadata_filtered$timepoint_numeric == tp]
    
    if(length(tp_samples) > 0) {
      # Get abundance data for this ORF at this timepoint
      orf_abundance <- significant_orf_data[orf_name, tp_samples]
      tp_metadata <- metadata_filtered[metadata_filtered$X %in% tp_samples, ]
      
      # Calculate group means
      celiac_samples <- tp_metadata$X[tp_metadata$Dx.Status == "CELIAC"]
      control_samples <- tp_metadata$X[tp_metadata$Dx.Status == "CONTROL"]
      
      if(length(celiac_samples) > 0 && length(control_samples) > 0) {
        celiac_mean <- mean(as.numeric(orf_abundance[celiac_samples]), na.rm = TRUE)
        control_mean <- mean(as.numeric(orf_abundance[control_samples]), na.rm = TRUE)
        
        # Store difference: CONTROL - CELIAC (positive = higher in CONTROL)
        # Handle sparse data by allowing small positive values
        if(!is.na(celiac_mean) && !is.na(control_mean)) {
          # Add small pseudocount to handle zeros
          celiac_adj <- celiac_mean + 0.001
          control_adj <- control_mean + 0.001
          heatmap_matrix[i, j] <- log2(control_adj) - log2(celiac_adj)
        }
      }
    }
  }
}

# Remove ORFs with too many missing values (more permissive for sparse data)
valid_orfs <- rowSums(!is.na(heatmap_matrix)) >= (ncol(heatmap_matrix) * 0.3)
heatmap_matrix_clean <- heatmap_matrix[valid_orfs, ]

cat("Created heatmap matrix with", nrow(heatmap_matrix_clean), "ORFs and", ncol(heatmap_matrix_clean), "timepoints\n")

# ============================================================================
# PART 4: CREATE ALL SIGNIFICANT ORFS HEATMAP
# ============================================================================

cat("\n4. Creating all significant ORFs heatmap...\n")

if(nrow(heatmap_matrix_clean) > 0) {
  
  # Set color limits to capture the data range
  data_range <- range(heatmap_matrix_clean, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  
  # Create color palette
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # Create heatmap
  png("all_significant_orfs_temporal_heatmap_US.png", width = 1400, height = 1000, res = 150)
  
  pheatmap(heatmap_matrix_clean,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "none",
           color = colors,
           breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Temporal Heatmap:", nrow(heatmap_matrix_clean), "Significant ORFs (US Cohort Only)\n",
                        "adj.p < 0.01, Mixed-effects models with patient random effects\n",
                        "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 8,
           fontsize_row = 6,
           fontsize_col = 10,
           show_rownames = FALSE,
           show_colnames = TRUE,
           border_color = NA)
  
  dev.off()
  
  cat("All significant ORFs heatmap saved\n")
}

# ============================================================================
# PART 5: CREATE TOP VARIABLE ORFS HEATMAP
# ============================================================================

cat("\n5. Creating top variable ORFs heatmap...\n")

if(nrow(heatmap_matrix_clean) > 0) {
  
  # Calculate variance for each ORF
  orf_variances <- apply(heatmap_matrix_clean, 1, function(x) var(x, na.rm = TRUE))
  orf_variances <- orf_variances[!is.na(orf_variances)]
  
  # Select top variable ORFs (top 200 or all if fewer)
  n_top <- min(200, length(orf_variances))
  top_variable_orfs <- names(sort(orf_variances, decreasing = TRUE)[1:n_top])
  
  # Create subset for top variable ORFs
  top_heatmap_matrix <- heatmap_matrix_clean[top_variable_orfs, ]
  
  cat("Selected top", n_top, "most variable ORFs\n")
  
  # Create heatmap with ORF names visible
  png("top_variable_significant_orfs_heatmap_US.png", width = 1400, height = max(800, n_top * 15), res = 150)
  
  pheatmap(top_heatmap_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "none",
           color = colors,
           breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Top", n_top, "Most Variable Significant ORFs (US Cohort Only)\n",
                        "adj.p < 0.01, Ordered by temporal variance\n",
                        "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 8,
           fontsize_row = max(4, min(8, 800/n_top)),
           fontsize_col = 10,
           show_rownames = TRUE,
           show_colnames = TRUE,
           border_color = NA)
  
  dev.off()
  
  cat("Top variable ORFs heatmap saved\n")
}

# ============================================================================
# PART 6: SAVE SUMMARY DATA
# ============================================================================

cat("\n6. Saving summary data...\n")

# Save heatmap matrix
write.csv(heatmap_matrix_clean, "temporal_heatmap_matrix_US.csv")

# Save significant ORFs info
significant_orfs_summary <- data.frame(
  orf_name = matching_orfs,
  adj_p_val = significant_orfs$adj.P.Val[match(matching_orfs, significant_orfs$X)],
  log_fc = significant_orfs$logFC[match(matching_orfs, significant_orfs$X)],
  in_heatmap = matching_orfs %in% rownames(heatmap_matrix_clean),
  stringsAsFactors = FALSE
)

write.csv(significant_orfs_summary, "significant_orfs_summary_US.csv", row.names = FALSE)

cat("Summary data saved\n")

# ============================================================================
# PART 7: FINAL SUMMARY
# ============================================================================

cat("\n=== TEMPORAL HEATMAP ANALYSIS COMPLETED (US COHORT) ===\n")
cat("Generated files:\n")
cat("- all_significant_orfs_temporal_heatmap_US.png: All significant ORFs heatmap\n")
cat("- top_variable_significant_orfs_heatmap_US.png: Top variable ORFs heatmap\n")
cat("- temporal_heatmap_matrix_US.csv: Heatmap data matrix\n")
cat("- significant_orfs_summary_US.csv: Summary of significant ORFs\n")

cat("\nKey statistics (US Cohort):\n")
cat("- Total ORFs analyzed:", nrow(limma_results), "\n")
cat("- Significant ORFs (adj.p < 0.01):", nrow(significant_orfs), 
    "(", round(nrow(significant_orfs)/nrow(limma_results)*100, 1), "%)\n")
cat("- ORFs in heatmap:", nrow(heatmap_matrix_clean), "\n")
cat("- Samples analyzed:", ncol(orf_data_filtered), "\n")
cat("- Timepoints:", length(timepoints), "\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")