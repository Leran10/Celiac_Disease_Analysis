#!/usr/bin/env Rscript
# Create temporal heatmap using limma mixed-effects models (same as _final plots)
# Author: Claude AI
# Date: 2025-08-05

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== CREATING TEMPORAL HEATMAP WITH LIMMA MODELS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading data files...\n")

# Load abundance data
abundance_data_raw <- read.csv("total_orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]  # Remove first column

# Fix sample names (remove X prefix if present)
sample_names <- colnames(abundance_data)
sample_names_clean <- gsub("^X", "", sample_names)
colnames(abundance_data) <- sample_names_clean

cat("Loaded abundance data with", nrow(abundance_data), "ORFs and", ncol(abundance_data), "samples\n")

# Load metadata
metadata_raw <- read.csv("total_metadata.clean.csv")
rownames(metadata_raw) <- metadata_raw$X
metadata <- metadata_raw  # Keep all columns for modeling

cat("Loaded metadata with", nrow(metadata), "samples\n")

# Load significant ORFs
topTable_results <- read.csv("total_limma_model_Sig_res_final.csv", row.names = 1)
significant_orfs <- rownames(topTable_results)
cat("Loaded", length(significant_orfs), "significant ORFs from topTable\n")

# ============================================================================
# PART 2: PREPARE DATA FOR MODELING
# ============================================================================

cat("\n2. Preparing data for limma modeling...\n")

# Match samples between abundance and metadata
abundance_samples <- colnames(abundance_data)
metadata_samples <- rownames(metadata)
matching_samples <- intersect(abundance_samples, metadata_samples)

cat("Sample matching: Found", length(matching_samples), "common samples\n")

if(length(matching_samples) == 0) {
  stop("No matching samples found!")
}

# Filter to matching samples
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

# Ensure Dx.Status is a factor with CELIAC as reference
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))
cat("Dx.Status levels:", paste(levels(metadata_filtered$Dx.Status), collapse = ", "), "\n")

# Check required columns
required_cols <- c("Dx.Status", "onset_timeline_numeric", "Country", "Sex", 
                   "Age.at.Gluten.Introduction..months.", "HLA.Category", 
                   "feeding_first_year", "Delivery.Mode", "patientID")

missing_cols <- setdiff(required_cols, colnames(metadata_filtered))
if(length(missing_cols) > 0) {
  cat("Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  stop("Required columns missing from metadata!")
}

cat("All required columns present\n")

# ============================================================================
# PART 3: CREATE DESIGN MATRIX
# ============================================================================

cat("\n3. Creating design matrix...\n")

# Create design matrix exactly as specified
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_numeric + Country + Sex + 
                                      Age.at.Gluten.Introduction..months. + HLA.Category + 
                                      feeding_first_year + Delivery.Mode, 
                                    data = metadata_filtered)

cat("Design matrix dimensions:", nrow(linear.model.design), "x", ncol(linear.model.design), "\n")
cat("Design matrix columns:", paste(colnames(linear.model.design), collapse = ", "), "\n")

# ============================================================================
# PART 4: FIT LIMMA MODELS
# ============================================================================

cat("\n4. Fitting limma mixed-effects models...\n")

# Convert abundance data to matrix for limma
abundance_matrix <- as.matrix(abundance_filtered)

cat("Abundance matrix dimensions:", nrow(abundance_matrix), "x", ncol(abundance_matrix), "\n")

# Fit limma model with voomLmFit and patient blocking
cat("Running voomLmFit with patient blocking...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)

cat("Running eBayes...\n")
linear.model.fit <- eBayes(linear.model.fit)

cat("Model fitting completed\n")

# ============================================================================
# PART 6: CALCULATE TEMPORAL DIFFERENCES (SAME AS _FINAL PLOTS)
# ============================================================================

cat("\n6. Calculating temporal differences using _final methodology...\n")

# Define the calculate_temporal_differences function (same as used for _final plots)
calculate_temporal_differences <- function(model_results, time_bins) {
  
  # Pre-allocate matrix
  heatmap_matrix <- matrix(NA, 
                          nrow = length(model_results), 
                          ncol = length(time_bins),
                          dimnames = list(names(model_results), 
                                         ifelse(time_bins == 0, "T0", paste0("T", time_bins))))
  
  for(i in 1:length(model_results)) {
    orf_id <- names(model_results)[i]
    orf_data <- model_results[[orf_id]]
    
    # Assign samples to exact timepoints (with tolerance for rounding)
    orf_data$time_bin <- NA
    for(timepoint in time_bins) {
      matches <- abs(orf_data$onset_timeline_numeric - timepoint) < 1
      orf_data$time_bin[matches] <- timepoint
    }
    
    for(j in 1:length(time_bins)) {
      timepoint <- time_bins[j]
      
      # Get fitted abundances for this timepoint
      time_data <- orf_data[!is.na(orf_data$time_bin) & orf_data$time_bin == timepoint, ]
      
      if(nrow(time_data) >= 2) {  # Require at least 2 samples total
        # Calculate mean fitted abundances by group
        group_means <- time_data %>%
          group_by(Dx.Status) %>%
          summarise(
            mean_fitted_abundance = mean(fitted_abundance, na.rm = TRUE),
            n_samples = n(),
            .groups = "drop"
          ) %>%
          filter(n_samples >= 1)
        
        # Calculate difference: CELIAC - CONTROL 
        # (Since CELIAC is reference group, we want positive = higher in CONTROL, negative = higher in CELIAC)
        # So we use celiac_fitted - control_fitted to get the right sign
        if(nrow(group_means) == 2) {
          celiac_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CELIAC"]
          control_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CONTROL"]
          
          if(length(celiac_fitted) > 0 && length(control_fitted) > 0) {
            heatmap_matrix[i, j] <- celiac_fitted - control_fitted
          }
        }
      }
    }
    
    # Progress update
    if(i %% 100 == 0) {
      cat("Processed", i, "of", length(model_results), "ORFs for temporal differences\n")
    }
  }
  
  return(heatmap_matrix)
}

# Define time bins (same as original analysis)
time_bins <- sort(unique(metadata_filtered$onset_timeline_numeric))
cat("Time bins:", paste(time_bins, collapse = ", "), "\n")

# Calculate temporal differences
heatmap_matrix <- calculate_temporal_differences(model_results, time_bins)

# Filter out ORFs with too many missing values
valid_orfs <- rowSums(!is.na(heatmap_matrix)) >= (ncol(heatmap_matrix) * 0.3)
heatmap_matrix_clean <- heatmap_matrix[valid_orfs, ]

cat("Final heatmap matrix:", nrow(heatmap_matrix_clean), "ORFs x", ncol(heatmap_matrix_clean), "timepoints\n")

# ============================================================================
# PART 7: CREATE HEATMAP VISUALIZATIONS (SAME AS _FINAL PLOTS)
# ============================================================================

cat("\n7. Creating heatmap visualizations...\n")

if(nrow(heatmap_matrix_clean) > 0) {
  
  # Set color limits
  data_range <- range(heatmap_matrix_clean, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  
  # Create color palette
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # ===== ALL SIGNIFICANT ORFS HEATMAP =====
  cat("Creating all significant ORFs heatmap...\n")
  
  png("all_significant_orfs_temporal_heatmap_limma_fitted.png", width = 1400, height = 1000, res = 150)
  
  pheatmap(heatmap_matrix_clean,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "none",
           color = colors,
           breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Temporal Heatmap:", nrow(heatmap_matrix_clean), "Significant ORFs (Limma Fitted Values)\n",
                        "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 8,
           fontsize_row = 6,
           fontsize_col = 10,
           show_rownames = FALSE,
           show_colnames = TRUE,
           border_color = NA)
  
  dev.off()
  cat("All ORFs heatmap saved\n")
  
  # ===== TOP VARIABLE ORFS HEATMAP =====
  cat("Creating top variable ORFs heatmap...\n")
  
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
  png("top_variable_significant_orfs_heatmap_limma_fitted.png", 
      width = 1400, height = max(800, n_top * 15), res = 150)
  
  pheatmap(top_heatmap_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "none",
           color = colors,
           breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Top", n_top, "Most Variable Significant ORFs (Limma Fitted Values)\n",
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
# PART 8: SAVE DATA AND SUMMARY
# ============================================================================

cat("\n8. Saving data and summary...\n")

# Save heatmap matrix
write.csv(heatmap_matrix_clean, "temporal_heatmap_matrix_limma_fitted.csv")

# Save model fitted values for future use
saveRDS(model_results, "limma_model_results_fitted_values.rds")

# Save summary
summary_stats <- data.frame(
  metric = c("Total ORFs in abundance data", "Significant ORFs from topTable", 
             "ORFs in final heatmap", "Total samples", "Unique timepoints"),
  count = c(nrow(abundance_matrix), length(significant_orfs), 
            nrow(heatmap_matrix_clean), ncol(abundance_matrix), length(time_bins))
)

write.csv(summary_stats, "temporal_heatmap_summary_limma_fitted.csv", row.names = FALSE)

cat("Data saved\n")

# ============================================================================
# PART 9: FINAL SUMMARY
# ============================================================================

cat("\n=== TEMPORAL HEATMAP ANALYSIS COMPLETED (LIMMA FITTED VALUES) ===\n")
cat("Generated files:\n")
cat("- all_significant_orfs_temporal_heatmap_limma_fitted.png: All significant ORFs heatmap\n")
cat("- top_variable_significant_orfs_heatmap_limma_fitted.png: Top variable ORFs heatmap\n")
cat("- temporal_heatmap_matrix_limma_fitted.csv: Heatmap data matrix\n")
cat("- limma_model_results_fitted_values.rds: Model fitted values for future use\n")
cat("- temporal_heatmap_summary_limma_fitted.csv: Summary statistics\n")

cat("\nKey statistics:\n")
print(summary_stats)

if(nrow(heatmap_matrix_clean) > 0) {
  cat("\nData range:", round(min(heatmap_matrix_clean, na.rm = TRUE), 3), 
      "to", round(max(heatmap_matrix_clean, na.rm = TRUE), 3), "\n")
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")