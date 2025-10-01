#!/usr/bin/env Rscript
# Complete ORF Consistency Analysis - All ORFs
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Recalculate consistency metrics using ALL ORFs (no sampling bias)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== COMPLETE ORF CONSISTENCY ANALYSIS - ALL ORFS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND FIT MODEL
# ============================================================================

cat("1. Loading data and fitting factor model...\n")

# Load abundance data
abundance_data_raw <- read.csv("US.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names
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

# Set reference levels with reversed timeline order
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))
metadata_filtered$onset_timeline_combined <- factor(metadata_filtered$onset_timeline_combined, 
                                                   levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0"))

# Create design matrix
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

# Convert abundance data to matrix
abundance_matrix <- as.matrix(abundance_filtered)

# Fit factor model
cat("Fitting factor model with patient blocking...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get fitted values
fitted_values <- fitted(linear.model.fit)

cat("Factor model fitted successfully\n")

# ============================================================================
# PART 2: LOAD ALL SIGNIFICANT AND NON-SIGNIFICANT ORFS
# ============================================================================

cat("2. Loading ALL significant and non-significant ORFs...\n")

# Load ALL union significant ORFs (no sampling)
union_orfs <- read.csv("factor_interaction_union_orfs.csv")
all_significant_orfs <- union_orfs$ORF

cat("ALL interaction-significant ORFs:", length(all_significant_orfs), "\n")

# Get all ORFs in the dataset
all_orfs_in_dataset <- rownames(abundance_matrix)
cat("Total ORFs in dataset:", length(all_orfs_in_dataset), "\n")

# Identify ALL non-significant ORFs (no sampling)
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
    
    if(i %% 1000 == 0) {
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

cat("3. Creating matrices for ALL ORFs...\n")

# Create matrices for ALL ORFs (no sampling)
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

cat("4. Running complete consistency analysis...\n")

all_significant_consistency <- calculate_complete_consistency(all_significant_matrix, "ALL Significant")
all_non_significant_consistency <- calculate_complete_consistency(all_non_significant_matrix, "ALL Non-significant")

# ============================================================================
# PART 7: COMPARISON AND SUMMARY
# ============================================================================

cat("5. Creating comprehensive comparison...\n")

# Create detailed comparison
complete_comparison <- data.frame(
  Metric = c(
    "Total ORFs in dataset", 
    "ALL significant ORFs", 
    "ALL non-significant ORFs",
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
    all_significant_consistency$total_orfs,
    all_non_significant_consistency$total_orfs,
    round(all_significant_consistency$mean_correlation, 3),
    round(all_non_significant_consistency$mean_correlation, 3),
    round(all_significant_consistency$consistency_score, 3),
    round(all_non_significant_consistency$consistency_score, 3),
    round(all_significant_consistency$mean_pattern_correlation, 3),
    round(all_non_significant_consistency$mean_pattern_correlation, 3),
    all_significant_consistency$high_variance_outliers,
    all_non_significant_consistency$high_variance_outliers,
    all_significant_consistency$low_correlation_outliers,
    all_non_significant_consistency$low_correlation_outliers
  ),
  stringsAsFactors = FALSE
)

# Save results
write.csv(complete_comparison, "complete_orf_consistency_comparison.csv", row.names = FALSE)

# Save detailed correlation data
sig_correlations_df <- data.frame(
  ORF = 1:length(all_significant_consistency$individual_correlations),
  Pattern_Correlation = all_significant_consistency$individual_correlations,
  Variance = all_significant_consistency$individual_variances[1:length(all_significant_consistency$individual_correlations)],
  Group = "Significant"
)

nonsig_correlations_df <- data.frame(
  ORF = 1:length(all_non_significant_consistency$individual_correlations),
  Pattern_Correlation = all_non_significant_consistency$individual_correlations,
  Variance = all_non_significant_consistency$individual_variances[1:length(all_non_significant_consistency$individual_correlations)],
  Group = "Non-significant"
)

all_correlations_df <- rbind(sig_correlations_df, nonsig_correlations_df)
write.csv(all_correlations_df, "all_orfs_pattern_correlations.csv", row.names = FALSE)

# Create comparison plots
png("complete_consistency_comparison.png", width = 2400, height = 1600, res = 150)
par(mfrow = c(2, 2))

# Pattern correlation histograms
hist(all_significant_consistency$individual_correlations, breaks = 50, 
     main = paste("Significant ORFs Pattern Correlations (n =", all_significant_consistency$total_orfs, ")"),
     xlab = "Correlation with Mean Pattern", col = "lightblue", border = "darkblue")
abline(v = mean(all_significant_consistency$individual_correlations), col = "red", lwd = 2)

hist(all_non_significant_consistency$individual_correlations, breaks = 50,
     main = paste("Non-significant ORFs Pattern Correlations (n =", all_non_significant_consistency$total_orfs, ")"),
     xlab = "Correlation with Mean Pattern", col = "lightcoral", border = "darkred")
abline(v = mean(all_non_significant_consistency$individual_correlations), col = "red", lwd = 2)

# Variance histograms
hist(all_significant_consistency$individual_variances, breaks = 50,
     main = "Significant ORFs Variance Distribution",
     xlab = "Variance", col = "lightgreen", border = "darkgreen")

hist(all_non_significant_consistency$individual_variances, breaks = 50,
     main = "Non-significant ORFs Variance Distribution", 
     xlab = "Variance", col = "lightyellow", border = "orange")

dev.off()

# ============================================================================
# PART 8: FINAL SUMMARY
# ============================================================================

cat("=== COMPLETE ORF CONSISTENCY ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- complete_orf_consistency_comparison.csv: Detailed comparison metrics\n")
cat("- all_orfs_pattern_correlations.csv: Individual ORF correlation data\n")
cat("- complete_consistency_comparison.png: Comparison plots\n\n")

cat("=== KEY FINDINGS (NO SAMPLING BIAS) ===\n")
print(complete_comparison)

cat("\nAnalysis completed at:", format(Sys.time()), "\n")