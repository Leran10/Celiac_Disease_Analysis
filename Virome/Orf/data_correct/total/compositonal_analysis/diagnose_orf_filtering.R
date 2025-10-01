#!/usr/bin/env Rscript
# Diagnose why ORFs were filtered out of temporal analysis
# Author: Claude AI

library(dplyr)

setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/compositonal_analysis")

cat("=== DIAGNOSING ORF FILTERING ===\n")

# Load data
limma_results <- read.csv("../limma_model_res.csv", row.names = 1)
metadata <- read.csv("../total.metadata.cleaned.csv", row.names = 1)
orf_data <- read.csv("../total.orf.abundance.table_0.75_prevFiltered_temporal_cleaned.csv", row.names = 1)

# Clean column names
orf_samples_clean <- gsub("^X", "", colnames(orf_data))
colnames(orf_data) <- orf_samples_clean

# Step 1: Significant ORFs
significant_orfs <- limma_results[limma_results$adj.P.Val < 0.05, ]
cat("Step 1 - Significant ORFs:", nrow(significant_orfs), "\n")

# Step 2: ORF-data intersection
common_samples <- intersect(orf_samples_clean, rownames(metadata))
orf_data_matched <- orf_data[, common_samples]
orf_ids_available <- intersect(rownames(significant_orfs), rownames(orf_data_matched))
cat("Step 2 - ORFs in abundance data:", length(orf_ids_available), "\n")
cat("         Lost at intersection:", nrow(significant_orfs) - length(orf_ids_available), "\n")

# Step 3: Check sample/patient requirements for each ORF
metadata_clean <- metadata[common_samples, ] %>%
  mutate(sample_id = rownames(.)) %>%
  filter(complete.cases(.))

failed_sample_req <- 0
failed_patient_req <- 0
both_failed <- 0

for(orf_id in orf_ids_available) {
  # Get abundance for this ORF
  orf_abundances <- orf_data_matched[orf_id, ]
  
  # Create data frame
  model_data <- data.frame(
    sample_id = names(orf_abundances),
    abundance = as.numeric(orf_abundances),
    stringsAsFactors = FALSE
  )
  
  # Merge with metadata
  model_data <- merge(model_data, metadata_clean, by = "sample_id")
  model_data <- model_data[complete.cases(model_data), ]
  
  # Check requirements
  sample_req <- nrow(model_data) >= 20
  patient_req <- length(unique(model_data$patientID)) >= 10
  
  if(!sample_req && !patient_req) {
    both_failed <- both_failed + 1
  } else if(!sample_req) {
    failed_sample_req <- failed_sample_req + 1
  } else if(!patient_req) {
    failed_patient_req <- failed_patient_req + 1
  }
}

cat("Step 3 - Sample/Patient requirements:\n")
cat("         Failed sample requirement (< 20 samples):", failed_sample_req, "\n")
cat("         Failed patient requirement (< 10 patients):", failed_patient_req, "\n") 
cat("         Failed both requirements:", both_failed, "\n")
cat("         Total failed:", failed_sample_req + failed_patient_req + both_failed, "\n")
cat("         Should pass to modeling:", length(orf_ids_available) - (failed_sample_req + failed_patient_req + both_failed), "\n")

# Final check
cat("\nActual modeled ORFs: 1880\n")
cat("Expected to pass: ", length(orf_ids_available) - (failed_sample_req + failed_patient_req + both_failed), "\n")

# Summary
cat("\n=== SUMMARY ===\n")
cat("Starting significant ORFs: 2707\n")
cat("ORFs not in abundance data: ", nrow(significant_orfs) - length(orf_ids_available), "\n")
cat("ORFs failing sample/patient req: ", failed_sample_req + failed_patient_req + both_failed, "\n")
cat("ORFs that should be modelable: ", length(orf_ids_available) - (failed_sample_req + failed_patient_req + both_failed), "\n")
cat("Actually modeled: 1880\n")
cat("Difference (likely convergence issues): ", 
    (length(orf_ids_available) - (failed_sample_req + failed_patient_req + both_failed)) - 1880, "\n")