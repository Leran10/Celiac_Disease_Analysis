#!/usr/bin/env Rscript
# Italy Factor Format Interaction Analysis
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Create factor format interaction analysis for Italy cohort (comparable to US)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== ITALY FACTOR FORMAT INTERACTION ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD AND PREPARE DATA
# ============================================================================

cat("1. Loading Italy cohort data...\n")

# Load abundance data
Italy.orf.abundance.raw <- read.csv("Italy.orf.abundance.clean.csv")
rownames(Italy.orf.abundance.raw) <- Italy.orf.abundance.raw$X
Italy.orf.abundance.clean <- Italy.orf.abundance.raw[, -1]

# Load metadata
Italy.metadata.raw <- read.csv("Italy.metadata.clean.csv")
rownames(Italy.metadata.raw) <- Italy.metadata.raw$X
Italy.metadata.clean <- Italy.metadata.raw

# Match samples between abundance and metadata
matching_samples <- intersect(colnames(Italy.orf.abundance.clean), rownames(Italy.metadata.clean))
Italy.orf.abundance.filtered <- Italy.orf.abundance.clean[, matching_samples]
Italy.metadata.filtered <- Italy.metadata.clean[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")
cat("Total ORFs:", nrow(Italy.orf.abundance.filtered), "\n")

# ============================================================================
# PART 2: SET FACTOR LEVELS (AS SPECIFIED)
# ============================================================================

cat("2. Setting factor levels as specified...\n")

# Set factor levels exactly as specified
Italy.metadata.clean$feeding_first_year <- factor(Italy.metadata.clean$feeding_first_year,
                                                  levels = c("Breast_fed","Formula","Breastmilk_and_formula"))

Italy.metadata.clean$HLA.Category <- factor(Italy.metadata.clean$HLA.Category,
                                           levels = c("Standard Risk","High Risk","Low/No Risk"))

Italy.metadata.clean$Sex <- factor(Italy.metadata.clean$Sex,
                                  levels = c("Female","Male"))

Italy.metadata.clean$Delivery.Mode <- factor(Italy.metadata.clean$Delivery.Mode,
                                            levels = c("Vaginal","C-Section"))

Italy.metadata.clean$Age.at.Gluten.Introduction..months. <- as.numeric(Italy.metadata.clean$Age.at.Gluten.Introduction..months.)

Italy.metadata.clean$Dx.Status <- factor(Italy.metadata.clean$Dx.Status,
                                        levels = c("CELIAC","CONTROL"))

Italy.metadata.clean$onset_timeline_combined <- factor(Italy.metadata.clean$onset_timeline_combined,
                                                      levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))

cat("Factor levels set:\n")
cat("- Dx.Status levels:", levels(Italy.metadata.clean$Dx.Status), "\n")
cat("- Timeline levels:", levels(Italy.metadata.clean$onset_timeline_combined), "\n")
cat("- HLA levels:", levels(Italy.metadata.clean$HLA.Category), "\n")
cat("- Feeding levels:", levels(Italy.metadata.clean$feeding_first_year), "\n")
cat("- Sex levels:", levels(Italy.metadata.clean$Sex), "\n")
cat("- Delivery levels:", levels(Italy.metadata.clean$Delivery.Mode), "\n")

# ============================================================================
# PART 3: CREATE DESIGN MATRIX AND FIT MODEL (AS SPECIFIED)
# ============================================================================

cat("3. Creating design matrix and fitting model...\n")

# Create design matrix exactly as specified
Italy.categorical.model.design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                                              Age.at.Gluten.Introduction..months. + HLA.Category + 
                                              feeding_first_year + Delivery.Mode,
                                              Italy.metadata.clean)

cat("Design matrix dimensions:", nrow(Italy.categorical.model.design), "x", ncol(Italy.categorical.model.design), "\n")
cat("Design matrix column names:\n")
print(colnames(Italy.categorical.model.design))

# Convert abundance data to matrix
Italy.orf.abundance.matrix <- as.matrix(Italy.orf.abundance.filtered)

# Fit model exactly as specified
cat("Fitting factor model with patient blocking...\n")
Italy.categorical.model.fit <- voomLmFit(Italy.orf.abundance.clean,
                                       Italy.categorical.model.design,
                                       block = Italy.metadata.clean$patientID,
                                       plot = FALSE)
Italy.categorical.model.fit <- eBayes(Italy.categorical.model.fit)

cat("Factor model fitted successfully\n")

# ============================================================================
# PART 4: EXTRACT INTERACTION COEFFICIENTS
# ============================================================================

cat("4. Extracting interaction coefficients...\n")

# Get coefficient names
coef_names <- colnames(Italy.categorical.model.fit$coefficients)
cat("All coefficients:\n")
print(coef_names)

# Find interaction coefficients (Dx.Status:onset_timeline_combined terms)
interaction_coefs <- coef_names[grep("Dx.StatusCONTROL:onset_timeline_combined", coef_names)]
cat("Interaction coefficients found:\n")
print(interaction_coefs)

if(length(interaction_coefs) == 0) {
  cat("WARNING: No interaction coefficients found with expected pattern\n")
  cat("Searching for alternative patterns...\n")
  # Try alternative patterns
  interaction_coefs_alt1 <- coef_names[grep("Dx.Status.*:.*onset_timeline", coef_names)]
  interaction_coefs_alt2 <- coef_names[grep("onset_timeline.*:.*Dx.Status", coef_names)]
  cat("Alternative pattern 1:", interaction_coefs_alt1, "\n")
  cat("Alternative pattern 2:", interaction_coefs_alt2, "\n")
  
  if(length(interaction_coefs_alt1) > 0) {
    interaction_coefs <- interaction_coefs_alt1
  } else if(length(interaction_coefs_alt2) > 0) {
    interaction_coefs <- interaction_coefs_alt2
  }
}

# ============================================================================
# PART 5: EXTRACT RESULTS FOR EACH INTERACTION TIMEPOINT
# ============================================================================

cat("5. Extracting results for each interaction timepoint...\n")

if(length(interaction_coefs) > 0) {
  # Create results for each interaction coefficient
  interaction_results_list <- list()
  
  for(coef_name in interaction_coefs) {
    cat("Processing coefficient:", coef_name, "\n")
    
    # Get topTable results for this coefficient
    coef_index <- which(colnames(Italy.categorical.model.fit$coefficients) == coef_name)
    
    if(length(coef_index) > 0) {
      results <- topTable(Italy.categorical.model.fit, coef = coef_index, 
                         number = Inf, adjust.method = "BH", sort.by = "P")
      
      # Add ORF ID column
      results$ORF <- rownames(results)
      
      # Save individual coefficient results
      clean_coef_name <- gsub("[^A-Za-z0-9_-]", "_", coef_name)
      filename <- paste0("italy_interaction_", clean_coef_name, "_results.csv")
      write.csv(results, filename, row.names = FALSE)
      
      # Extract significant ORFs (padj < 0.05)
      significant_orfs <- results[results$adj.P.Val < 0.05, ]
      sig_filename <- paste0("italy_significant_orfs_", clean_coef_name, ".csv")
      write.csv(significant_orfs, sig_filename, row.names = FALSE)
      
      cat("- Coefficient:", coef_name, "\n")
      cat("  - Total ORFs:", nrow(results), "\n")
      cat("  - Significant ORFs (padj < 0.05):", nrow(significant_orfs), "\n")
      cat("  - Results saved to:", filename, "\n")
      cat("  - Significant ORFs saved to:", sig_filename, "\n\n")
      
      # Store for union analysis
      interaction_results_list[[coef_name]] <- list(
        all_results = results,
        significant_orfs = significant_orfs$ORF
      )
    }
  }
  
  # ============================================================================
  # PART 6: CREATE UNION OF SIGNIFICANT ORFS
  # ============================================================================
  
  cat("6. Creating union of significant ORFs across all interaction timepoints...\n")
  
  # Collect all significant ORFs
  all_significant_orfs <- c()
  timepoint_specific_orfs <- list()
  
  for(coef_name in names(interaction_results_list)) {
    sig_orfs <- interaction_results_list[[coef_name]]$significant_orfs
    all_significant_orfs <- c(all_significant_orfs, sig_orfs)
    timepoint_specific_orfs[[coef_name]] <- sig_orfs
    
    cat("- ", coef_name, ":", length(sig_orfs), "significant ORFs\n")
  }
  
  # Get unique ORFs (union)
  union_significant_orfs <- unique(all_significant_orfs)
  
  cat("Total unique significant ORFs (union):", length(union_significant_orfs), "\n")
  
  # Save union results
  union_df <- data.frame(
    ORF = union_significant_orfs,
    stringsAsFactors = FALSE
  )
  write.csv(union_df, "italy_factor_interaction_union_orfs.csv", row.names = FALSE)
  
  # Create breakdown analysis
  breakdown_results <- data.frame(
    Timepoint = names(timepoint_specific_orfs),
    Number_of_ORFs = sapply(timepoint_specific_orfs, length),
    stringsAsFactors = FALSE
  )
  
  # Add union summary
  breakdown_results <- rbind(breakdown_results, 
                           data.frame(Timepoint = "Union_All_Timepoints", 
                                    Number_of_ORFs = length(union_significant_orfs)))
  
  write.csv(breakdown_results, "italy_factor_interaction_breakdown.csv", row.names = FALSE)
  
  # ============================================================================
  # PART 7: SAVE COMPLETE MODEL RESULTS
  # ============================================================================
  
  cat("7. Saving complete model results...\n")
  
  # Extract all coefficients results
  all_results <- topTable(Italy.categorical.model.fit, number = Inf, 
                         adjust.method = "BH", sort.by = "none")
  all_results$ORF <- rownames(all_results)
  
  write.csv(all_results, "italy_factor_complete_model_results.csv", row.names = FALSE)
  
  # Save fitted values
  fitted_values <- fitted(Italy.categorical.model.fit)
  saveRDS(fitted_values, "italy_factor_model_fitted_values.rds")
  
  cat("Complete model results saved\n")
  
} else {
  cat("ERROR: No interaction coefficients found. Cannot proceed with analysis.\n")
  cat("Please check the model design and factor levels.\n")
}

# ============================================================================
# PART 8: SUMMARY
# ============================================================================

cat("\n=== ITALY FACTOR INTERACTION ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")

if(length(interaction_coefs) > 0) {
  cat("- italy_factor_interaction_union_orfs.csv: Union of all significant interaction ORFs\n")
  cat("- italy_factor_interaction_breakdown.csv: Breakdown by timepoint\n")
  cat("- italy_factor_complete_model_results.csv: Complete model results\n")
  cat("- italy_factor_model_fitted_values.rds: Fitted values for heatmap generation\n")
  
  for(coef_name in interaction_coefs) {
    clean_coef_name <- gsub("[^A-Za-z0-9_-]", "_", coef_name)
    cat("- italy_interaction_", clean_coef_name, "_results.csv: Individual coefficient results\n")
    cat("- italy_significant_orfs_", clean_coef_name, ".csv: Significant ORFs for this coefficient\n")
  }
  
  cat("\nSummary:\n")
  print(breakdown_results)
  
  cat("\nNow ready for consistency analysis comparable to US cohort!\n")
} else {
  cat("Analysis failed - no interaction coefficients found\n")
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")