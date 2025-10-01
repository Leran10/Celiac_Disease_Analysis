#!/usr/bin/env Rscript
# Italy Complete ORF Consistency Analysis - All ORFs
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Analyze ORF consistency patterns in Italy cohort (comparable to US)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== ITALY COMPLETE ORF CONSISTENCY ANALYSIS - ALL ORFS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND USE EXISTING FACTOR MODEL
# ============================================================================

cat("1. Loading Italy data and fitted factor model...\n")

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
rownames(metadata_raw) <- metadata_raw[, 1]  # Use first column (sample names)
metadata <- metadata_raw

# Match samples (abundance column names should match metadata row names)
cat("Abundance columns (first 5):", head(colnames(abundance_data), 5), "\n")
cat("Metadata rows (first 5):", head(rownames(metadata), 5), "\n")

matching_samples <- intersect(colnames(abundance_data), rownames(metadata))
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")

# Set factor levels (same as previous analysis)
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))
metadata_filtered$onset_timeline_combined <- factor(metadata_filtered$onset_timeline_combined, 
                                                   levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0"))

cat("Timeline levels (reversed for consistency):", levels(metadata_filtered$onset_timeline_combined), "\n")

# Recreate the factor model to get fitted values
cat("Recreating factor model to get fitted values...\n")

# Set factor levels exactly as in the previous analysis
metadata_filtered$feeding_first_year <- factor(metadata_filtered$feeding_first_year,
                                               levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
metadata_filtered$HLA.Category <- factor(metadata_filtered$HLA.Category,
                                         levels = c("Standard Risk","High Risk","Low/No Risk"))
metadata_filtered$Sex <- factor(metadata_filtered$Sex,
                                levels = c("Female","Male"))
metadata_filtered$Delivery.Mode <- factor(metadata_filtered$Delivery.Mode,
                                          levels = c("Vaginal","C-Section"))
metadata_filtered$Age.at.Gluten.Introduction..months. <- as.numeric(metadata_filtered$Age.at.Gluten.Introduction..months.)

# Create design matrix
italy_design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                            Age.at.Gluten.Introduction..months. + HLA.Category + 
                            feeding_first_year + Delivery.Mode,
                            metadata_filtered)

# Convert abundance data to matrix
abundance_matrix <- as.matrix(abundance_filtered)

# Fit model
italy_model_fit <- voomLmFit(abundance_matrix, italy_design, 
                            block = metadata_filtered$patientID, plot = FALSE)
italy_model_fit <- eBayes(italy_model_fit)

# Get fitted values
fitted_values <- fitted(italy_model_fit)

cat("Factor model recreated and fitted values obtained\n")

# ============================================================================
# PART 2: LOAD SIGNIFICANT AND NON-SIGNIFICANT ORFS
# ============================================================================

cat("2. Loading ALL significant and non-significant ORFs...\n")

# Load ALL union significant ORFs
union_orfs <- read.csv("italy_factor_interaction_union_orfs.csv")
all_significant_orfs <- union_orfs$ORF

cat("ALL interaction-significant ORFs:", length(all_significant_orfs), "\n")

# Get all ORFs in the dataset
all_orfs_in_dataset <- rownames(abundance_filtered)
cat("Total ORFs in dataset:", length(all_orfs_in_dataset), "\n")

# Identify ALL non-significant ORFs
all_non_significant_orfs <- setdiff(all_orfs_in_dataset, all_significant_orfs)
cat("ALL non-significant ORFs:", length(all_non_significant_orfs), "\n")

# ============================================================================
# PART 3: CREATE HEATMAP MATRIX FUNCTION
# ============================================================================

create_complete_heatmap_matrix <- function(orf_list, name_suffix) {
  cat("Creating heatmap matrix for", length(orf_list), name_suffix, "ORFs...\n")
  
  all_timepoints <- levels(metadata_filtered$onset_timeline_combined)
  
  heatmap_matrix <- matrix(NA, 
                          nrow = length(orf_list), 
                          ncol = length(all_timepoints),
                          dimnames = list(orf_list, all_timepoints))
  
  for(i in 1:length(orf_list)) {
    orf_id <- orf_list[i]
    
    if(orf_id %in% rownames(fitted_values)) {
      orf_fitted <- fitted_values[orf_id, ]
      
      for(j in 1:length(all_timepoints)) {
        timepoint <- all_timepoints[j]
        
        timepoint_samples <- metadata_filtered$onset_timeline_combined == timepoint
        
        if(sum(timepoint_samples) > 0) {
          timepoint_data <- data.frame(
            fitted_abundance = as.numeric(orf_fitted[timepoint_samples]),
            Dx.Status = metadata_filtered$Dx.Status[timepoint_samples]
          )
          
          group_means <- timepoint_data %>%
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
              heatmap_matrix[i, j] <- control_fitted - celiac_fitted
            }
          }
        }
      }
    }
    
    if(i %% 500 == 0) {
      cat("Processed", i, "of", length(orf_list), name_suffix, "ORFs\n")
    }
  }
  
  # Remove ORFs with too many missing values
  valid_orfs <- rowSums(!is.na(heatmap_matrix)) >= (ncol(heatmap_matrix) * 0.5)
  heatmap_matrix_clean <- heatmap_matrix[valid_orfs, ]
  
  cat("Final", name_suffix, "matrix:", nrow(heatmap_matrix_clean), "x", ncol(heatmap_matrix_clean), "\n")
  
  return(heatmap_matrix_clean)
}

# ============================================================================
# PART 4: CREATE MATRICES FOR ALL ORFS
# ============================================================================

cat("3. Creating matrices for ALL Italy ORFs...\n")

# Create matrices for ALL ORFs
all_significant_matrix <- create_complete_heatmap_matrix(all_significant_orfs, "significant")
all_non_significant_matrix <- create_complete_heatmap_matrix(all_non_significant_orfs, "non-significant")

# ============================================================================
# PART 5: CONSISTENCY ANALYSIS FUNCTIONS
# ============================================================================

calculate_complete_consistency <- function(matrix_data, name) {
  cat("=== Complete Consistency Analysis for", name, "ORFs ===\n")
  
  # 1. Calculate correlation between timepoints
  timepoint_correlations <- cor(matrix_data, use = "complete.obs")
  mean_correlation <- mean(timepoint_correlations[upper.tri(timepoint_correlations)], na.rm = TRUE)
  
  # 2. Calculate proportion of ORFs following same direction at each timepoint
  timepoint_directions <- apply(matrix_data, 2, function(x) {
    pos_prop <- sum(x > 0, na.rm = TRUE) / sum(!is.na(x))
    neg_prop <- sum(x < 0, na.rm = TRUE) / sum(!is.na(x))
    return(c(positive = pos_prop, negative = neg_prop))
  })
  
  # 3. Overall consistency score
  consistency_scores <- apply(timepoint_directions, 2, max)
  overall_consistency <- mean(consistency_scores)
  
  # 4. Calculate mean pattern and individual ORF correlations
  mean_pattern <- colMeans(matrix_data, na.rm = TRUE)
  
  individual_pattern_correlations <- c()
  for(i in 1:nrow(matrix_data)) {
    orf_values <- matrix_data[i, ]
    orf_values_clean <- orf_values[!is.na(orf_values)]
    mean_pattern_clean <- mean_pattern[names(mean_pattern) %in% names(orf_values_clean)]
    
    if(length(orf_values_clean) >= 3 && length(mean_pattern_clean) >= 3) {
      pattern_cor <- cor(orf_values_clean, mean_pattern_clean, use = "complete.obs")
      if(!is.na(pattern_cor)) {
        individual_pattern_correlations <- c(individual_pattern_correlations, pattern_cor)
      }
    }
  }
  
  mean_pattern_correlation <- mean(individual_pattern_correlations, na.rm = TRUE)
  
  # 5. Variance analysis
  orf_variances <- apply(matrix_data, 1, function(x) var(x, na.rm = TRUE))
  orf_variances_clean <- orf_variances[!is.na(orf_variances)]
  
  # Outlier identification (5% thresholds)
  n_total_orfs <- length(individual_pattern_correlations)
  n_high_variance <- sum(orf_variances_clean >= quantile(orf_variances_clean, 0.95, na.rm = TRUE), na.rm = TRUE)
  n_low_correlation <- sum(individual_pattern_correlations <= quantile(individual_pattern_correlations, 0.05, na.rm = TRUE), na.rm = TRUE)
  
  cat("Total ORFs analyzed:", n_total_orfs, "\n")
  cat("Mean inter-timepoint correlation:", round(mean_correlation, 3), "\n")
  cat("Overall consistency score:", round(overall_consistency, 3), "\n")
  cat("Mean pattern correlation:", round(mean_pattern_correlation, 3), "\n")
  cat("High variance outliers (top 5%):", n_high_variance, "\n")
  cat("Low correlation outliers (bottom 5%):", n_low_correlation, "\n")
  cat("Data range:", round(range(matrix_data, na.rm = TRUE), 2), "\n")
  cat("Variance range:", round(range(orf_variances_clean, na.rm = TRUE), 3), "\n")
  cat("Pattern correlation range:", round(range(individual_pattern_correlations, na.rm = TRUE), 3), "\n\n")
  
  return(list(
    total_orfs = n_total_orfs,
    mean_correlation = mean_correlation,
    consistency_score = overall_consistency,
    mean_pattern_correlation = mean_pattern_correlation,
    high_variance_outliers = n_high_variance,
    low_correlation_outliers = n_low_correlation,
    data_range = range(matrix_data, na.rm = TRUE),
    variance_range = range(orf_variances_clean, na.rm = TRUE),
    correlation_range = range(individual_pattern_correlations, na.rm = TRUE),
    individual_correlations = individual_pattern_correlations,
    individual_variances = orf_variances_clean
  ))
}

# ============================================================================
# PART 6: RUN COMPLETE CONSISTENCY ANALYSIS
# ============================================================================

cat("4. Running complete consistency analysis for Italy...\n")

italy_significant_consistency <- calculate_complete_consistency(all_significant_matrix, "ALL Significant (Italy)")
italy_non_significant_consistency <- calculate_complete_consistency(all_non_significant_matrix, "ALL Non-significant (Italy)")

# ============================================================================
# PART 7: COMPARISON AND SUMMARY
# ============================================================================

cat("5. Creating comprehensive Italy consistency comparison...\n")

# Create detailed comparison
italy_comparison <- data.frame(
  Metric = c(
    "Total ORFs in Italy dataset", 
    "ALL significant ORFs (Italy)", 
    "ALL non-significant ORFs (Italy)",
    "Significant ORFs analyzed",
    "Non-significant ORFs analyzed",
    "Significant: Mean inter-timepoint correlation",
    "Non-significant: Mean inter-timepoint correlation", 
    "Significant: Consistency score",
    "Non-significant: Consistency score",
    "Significant: Mean pattern correlation",
    "Non-significant: Mean pattern correlation",
    "Significant: High variance outliers",
    "Non-significant: High variance outliers",
    "Significant: Low correlation outliers", 
    "Non-significant: Low correlation outliers"
  ),
  Value = c(
    length(all_orfs_in_dataset),
    length(all_significant_orfs),
    length(all_non_significant_orfs),
    italy_significant_consistency$total_orfs,
    italy_non_significant_consistency$total_orfs,
    round(italy_significant_consistency$mean_correlation, 3),
    round(italy_non_significant_consistency$mean_correlation, 3),
    round(italy_significant_consistency$consistency_score, 3),
    round(italy_non_significant_consistency$consistency_score, 3),
    round(italy_significant_consistency$mean_pattern_correlation, 3),
    round(italy_non_significant_consistency$mean_pattern_correlation, 3),
    italy_significant_consistency$high_variance_outliers,
    italy_non_significant_consistency$high_variance_outliers,
    italy_significant_consistency$low_correlation_outliers,
    italy_non_significant_consistency$low_correlation_outliers
  ),
  stringsAsFactors = FALSE
)

# Save results
write.csv(italy_comparison, "italy_complete_orf_consistency_comparison.csv", row.names = FALSE)

# Save detailed correlation data
italy_sig_correlations_df <- data.frame(
  ORF = 1:length(italy_significant_consistency$individual_correlations),
  Pattern_Correlation = italy_significant_consistency$individual_correlations,
  Variance = italy_significant_consistency$individual_variances[1:length(italy_significant_consistency$individual_correlations)],
  Group = "Significant"
)

italy_nonsig_correlations_df <- data.frame(
  ORF = 1:length(italy_non_significant_consistency$individual_correlations),
  Pattern_Correlation = italy_non_significant_consistency$individual_correlations,
  Variance = italy_non_significant_consistency$individual_variances[1:length(italy_non_significant_consistency$individual_correlations)],
  Group = "Non-significant"
)

italy_all_correlations_df <- rbind(italy_sig_correlations_df, italy_nonsig_correlations_df)
write.csv(italy_all_correlations_df, "italy_all_orfs_pattern_correlations.csv", row.names = FALSE)

# Create comparison plots
png("italy_complete_consistency_comparison.png", width = 2400, height = 1600, res = 150)
par(mfrow = c(2, 2))

# Pattern correlation histograms
hist(italy_significant_consistency$individual_correlations, breaks = 50, 
     main = paste("Italy Significant ORFs Pattern Correlations (n =", italy_significant_consistency$total_orfs, ")"),
     xlab = "Correlation with Mean Pattern", col = "lightblue", border = "darkblue")
abline(v = mean(italy_significant_consistency$individual_correlations), col = "red", lwd = 2)

hist(italy_non_significant_consistency$individual_correlations, breaks = 50,
     main = paste("Italy Non-significant ORFs Pattern Correlations (n =", italy_non_significant_consistency$total_orfs, ")"),
     xlab = "Correlation with Mean Pattern", col = "lightcoral", border = "darkred")
abline(v = mean(italy_non_significant_consistency$individual_correlations), col = "red", lwd = 2)

# Variance histograms
hist(italy_significant_consistency$individual_variances, breaks = 50,
     main = "Italy Significant ORFs Variance Distribution",
     xlab = "Variance", col = "lightgreen", border = "darkgreen")

hist(italy_non_significant_consistency$individual_variances, breaks = 50,
     main = "Italy Non-significant ORFs Variance Distribution", 
     xlab = "Variance", col = "lightyellow", border = "orange")

dev.off()

# Create heatmaps for visualization
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Significant ORFs heatmap (sample for visualization)
sample_size_sig <- min(500, nrow(all_significant_matrix))
sample_indices_sig <- sample(1:nrow(all_significant_matrix), sample_size_sig)
sample_significant_matrix <- all_significant_matrix[sample_indices_sig, ]

png("italy_consistency_significant_sample_heatmap.png", width = 2000, height = 1500, res = 150)
data_range_sig <- range(sample_significant_matrix, na.rm = TRUE)
max_abs_sig <- max(abs(data_range_sig))
pheatmap(sample_significant_matrix,
         cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
         color = colors, breaks = seq(-max_abs_sig, max_abs_sig, length.out = 101),
         main = paste("Italy: Sample of Interaction-Significant ORFs (n =", nrow(sample_significant_matrix), ")\nRed = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 14, show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
dev.off()

# Non-significant ORFs heatmap (sample)
sample_size_nonsig <- min(500, nrow(all_non_significant_matrix))
sample_indices_nonsig <- sample(1:nrow(all_non_significant_matrix), sample_size_nonsig)
sample_non_significant_matrix <- all_non_significant_matrix[sample_indices_nonsig, ]

png("italy_consistency_non_significant_sample_heatmap.png", width = 2000, height = 1500, res = 150)
data_range_nonsig <- range(sample_non_significant_matrix, na.rm = TRUE)
max_abs_nonsig <- max(abs(data_range_nonsig))
pheatmap(sample_non_significant_matrix,
         cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
         color = colors, breaks = seq(-max_abs_nonsig, max_abs_nonsig, length.out = 101),
         main = paste("Italy: Sample of Non-Interaction-Significant ORFs (n =", nrow(sample_non_significant_matrix), ")\nRed = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 14, show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
dev.off()

# ============================================================================
# PART 8: FINAL SUMMARY
# ============================================================================

cat("\n=== ITALY COMPLETE ORF CONSISTENCY ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- italy_complete_orf_consistency_comparison.csv: Detailed comparison metrics\n")
cat("- italy_all_orfs_pattern_correlations.csv: Individual ORF correlation data\n")
cat("- italy_complete_consistency_comparison.png: Comparison plots\n")
cat("- italy_consistency_significant_sample_heatmap.png: Sample of significant ORFs\n")
cat("- italy_consistency_non_significant_sample_heatmap.png: Sample of non-significant ORFs\n")

cat("\n=== ITALY KEY FINDINGS (NO SAMPLING BIAS) ===\n")
print(italy_comparison)

cat("\nAnalysis completed at:", format(Sys.time()), "\n")