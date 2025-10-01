#!/usr/bin/env Rscript
# Italy Cohort Temporal Heatmap Analysis - Top 200 Most Variable ORFs
# Author: Claude AI
# Date: 2025-08-11
# Purpose: Create temporal heatmap for Italy cohort to compare with US cohort patterns

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== ITALY COHORT TEMPORAL HEATMAP ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND FIT LIMMA MODELS
# ============================================================================

cat("1. Loading Italy cohort data and fitting limma models...\n")

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

# Ensure Dx.Status is a factor with CELIAC as reference (consistent with US cohort)
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
cat("Running voomLmFit with patient blocking...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get all ORF results (not just significant ones)
topTable_results <- topTable(linear.model.fit, coef = "Dx.StatusCONTROL:onset_timeline_numeric", 
                             number = Inf, sort.by = "P")
significant_orfs <- topTable_results[topTable_results$adj.P.Val < 0.05, ]

cat("Identified", nrow(significant_orfs), "significant ORFs at adj.p < 0.05 (interaction term)\n")
cat("Total ORFs available for analysis:", nrow(topTable_results), "\n")
cat("Will create heatmap using top 200 most variable ORFs regardless of significance\n")

# ============================================================================
# PART 2: CREATE TEMPORAL HEATMAPS WITH LIMMA FITTED VALUES
# ============================================================================

cat("\n2. Creating temporal heatmaps with limma fitted values...\n")

# Get fitted values from the model
fitted_values <- fitted(linear.model.fit)

# Create model_results list for ALL ORFs (not just significant ones)
all_orf_names <- rownames(fitted_values)
available_orfs <- rownames(fitted_values)

cat("Creating model_results list for all", length(all_orf_names), "ORFs...\n")

model_results <- list()
for(orf_id in all_orf_names) {
  orf_fitted <- fitted_values[orf_id, ]
  model_results[[orf_id]] <- data.frame(
    onset_timeline_numeric = metadata_filtered$onset_timeline_numeric,
    Dx.Status = metadata_filtered$Dx.Status,
    fitted_abundance = as.numeric(orf_fitted),
    stringsAsFactors = FALSE
  )
}

cat("Created model_results list with", length(model_results), "ORFs\n")

# Define the calculate_temporal_differences function (same as US cohort)
calculate_temporal_differences <- function(model_results, time_bins) {
  heatmap_matrix <- matrix(NA, 
                          nrow = length(model_results), 
                          ncol = length(time_bins),
                          dimnames = list(names(model_results), 
                                         ifelse(time_bins == 0, "T0", paste0("T", time_bins))))
  
  for(i in 1:length(model_results)) {
    orf_id <- names(model_results)[i]
    orf_data <- model_results[[orf_id]]
    
    orf_data$time_bin <- NA
    for(timepoint in time_bins) {
      matches <- abs(orf_data$onset_timeline_numeric - timepoint) < 1
      orf_data$time_bin[matches] <- timepoint
    }
    
    for(j in 1:length(time_bins)) {
      timepoint <- time_bins[j]
      time_data <- orf_data[!is.na(orf_data$time_bin) & orf_data$time_bin == timepoint, ]
      
      if(nrow(time_data) >= 2) {
        group_means <- time_data %>%
          group_by(Dx.Status) %>%
          summarise(
            mean_fitted_abundance = mean(fitted_abundance, na.rm = TRUE),
            n_samples = n(),
            .groups = "drop"
          ) %>%
          filter(n_samples >= 1)
        
        if(nrow(group_means) == 2) {
          celiac_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CELIAC"]
          control_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CONTROL"]
          
          if(length(celiac_fitted) > 0 && length(control_fitted) > 0) {
            # Use same calculation as US cohort: control_fitted - celiac_fitted
            # This accounts for CELIAC reference level
            heatmap_matrix[i, j] <- control_fitted - celiac_fitted
          }
        }
      }
    }
    
    if(i %% 500 == 0) {
      cat("Processed", i, "of", length(model_results), "ORFs for temporal differences\n")
    }
  }
  
  return(heatmap_matrix)
}

# Calculate temporal differences for all ORFs
time_bins <- sort(unique(metadata_filtered$onset_timeline_numeric))
cat("Time bins:", paste(time_bins, collapse = ", "), "\n")

heatmap_matrix <- calculate_temporal_differences(model_results, time_bins)

# Filter out ORFs with too many missing values (same threshold as US cohort)
valid_orfs <- rowSums(!is.na(heatmap_matrix)) >= (ncol(heatmap_matrix) * 0.3)
heatmap_matrix_clean <- heatmap_matrix[valid_orfs, ]

cat("Final heatmap matrix:", nrow(heatmap_matrix_clean), "ORFs x", ncol(heatmap_matrix_clean), "timepoints\n")

# Save the heatmap matrix
write.csv(heatmap_matrix_clean, "temporal_heatmap_matrix_Italy_limma_fitted.csv")

# ============================================================================
# PART 3: CREATE TOP 200 VARIABLE ORFS HEATMAP
# ============================================================================

cat("\n3. Creating top 200 most variable ORFs heatmap...\n")

if(nrow(heatmap_matrix_clean) > 0) {
  # Calculate variances for all ORFs
  orf_variances <- apply(heatmap_matrix_clean, 1, function(x) var(x, na.rm = TRUE))
  orf_variances <- orf_variances[!is.na(orf_variances)]
  
  # Select top 200 most variable ORFs (or all if less than 200)
  n_top <- min(200, length(orf_variances))
  top_variable_orfs <- names(sort(orf_variances, decreasing = TRUE)[1:n_top])
  top_heatmap_matrix <- heatmap_matrix_clean[top_variable_orfs, ]
  
  cat("Selected top", n_top, "most variable ORFs for heatmap\n")
  
  # Set up color scheme (same as US cohort)
  data_range <- range(heatmap_matrix_clean, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  cat("Data range for heatmap:", round(data_range[1], 2), "to", round(data_range[2], 2), "\n")
  
  # Generate the temporal heatmap (same dimensions as US cohort)
  png("top_variable_significant_orfs_heatmap_Italy_limma_fitted.png", 
      width = 2200, height = max(800, n_top * 15), res = 150)
  pheatmap(top_heatmap_matrix,
           cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
           color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Italy Cohort: Top", n_top, "Most Variable ORFs (Limma Fitted Values)\n",
                        "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 18, fontsize_row = max(10, min(14, 1200/n_top)), fontsize_col = 20,
           cex_main = 1.5,
           show_rownames = TRUE, show_colnames = TRUE, border_color = NA)
  dev.off()
  
  cat("Temporal heatmap created successfully\n")
  
  # Also create a version with all valid ORFs for comparison
  if(nrow(heatmap_matrix_clean) < 1000) {  # Only if manageable number
    png("all_orfs_temporal_heatmap_Italy_limma_fitted.png", 
        width = 2200, height = max(1200, nrow(heatmap_matrix_clean) * 8), res = 150)
    pheatmap(heatmap_matrix_clean,
             cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
             color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
             main = paste("Italy Cohort: All", nrow(heatmap_matrix_clean), "ORFs - Temporal Patterns (Limma Fitted Values)\n",
                          "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
             fontsize = 14, fontsize_row = 6, fontsize_col = 16,
             cex_main = 1.2,
             show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
    dev.off()
    
    cat("All ORFs heatmap also created for reference\n")
  }
  
} else {
  cat("ERROR: No valid ORFs found for heatmap creation\n")
}

# ============================================================================
# PART 4: SUMMARY AND COMPARISON DATA
# ============================================================================

cat("\n4. Creating summary statistics for comparison with US cohort...\n")

# Calculate summary statistics
summary_stats <- data.frame(
  Metric = c("Total Samples", "Total Patients", "Total ORFs", "Valid ORFs for Heatmap", 
             "Top Variable ORFs", "Significant ORFs (Interaction)", "Data Range Min", "Data Range Max"),
  Italy_Cohort = c(nrow(metadata_filtered), 
                   length(unique(metadata_filtered$patientID)),
                   nrow(abundance_matrix),
                   nrow(heatmap_matrix_clean),
                   n_top,
                   nrow(significant_orfs),
                   round(min(heatmap_matrix_clean, na.rm = TRUE), 2),
                   round(max(heatmap_matrix_clean, na.rm = TRUE), 2)),
  stringsAsFactors = FALSE
)

write.csv(summary_stats, "italy_temporal_heatmap_summary.csv", row.names = FALSE)

# Save the model results for future use
saveRDS(model_results, "limma_model_results_fitted_values_Italy.rds")

# Print summary
cat("\n=== ITALY COHORT TEMPORAL HEATMAP ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- top_variable_significant_orfs_heatmap_Italy_limma_fitted.png: Top", n_top, "variable ORFs heatmap\n")
if(nrow(heatmap_matrix_clean) < 1000) {
  cat("- all_orfs_temporal_heatmap_Italy_limma_fitted.png: All ORFs heatmap for reference\n")
}
cat("- temporal_heatmap_matrix_Italy_limma_fitted.csv: Heatmap matrix data\n")
cat("- italy_temporal_heatmap_summary.csv: Summary statistics\n")
cat("- limma_model_results_fitted_values_Italy.rds: Model results for future use\n")

cat("\nKey Results Summary:\n")
cat("- Total samples:", nrow(metadata_filtered), "\n")
cat("- Total patients:", length(unique(metadata_filtered$patientID)), "\n")
cat("- Total ORFs analyzed:", nrow(abundance_matrix), "\n")
cat("- Valid ORFs for heatmap:", nrow(heatmap_matrix_clean), "\n")
cat("- Top variable ORFs displayed:", n_top, "\n")
cat("- Significant ORFs (interaction):", nrow(significant_orfs), "\n")
cat("- Data range (heatmap):", round(min(heatmap_matrix_clean, na.rm = TRUE), 2), "to", 
    round(max(heatmap_matrix_clean, na.rm = TRUE), 2), "\n")

cat("\nReady for comparison with US cohort temporal patterns!\n")
cat("Analysis completed at:", format(Sys.time()), "\n")