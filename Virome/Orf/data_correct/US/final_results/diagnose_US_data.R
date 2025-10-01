#!/usr/bin/env Rscript
# Diagnostic script for US data
# Author: Claude AI

library(dplyr)

setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/US/final_results")

# Load data
limma_results <- read.csv("../US.limma_model_res.csv")
orf_data_raw <- read.csv("../US.orf.abundance.clean.csv")
rownames(orf_data_raw) <- orf_data_raw$X
orf_data <- orf_data_raw[, -1]
orf_samples_clean <- gsub("^X", "", colnames(orf_data))
colnames(orf_data) <- orf_samples_clean

metadata <- read.csv("../US.metadata.clean.csv")

cat("=== DIAGNOSTIC ANALYSIS ===\n")
cat("Limma results:", nrow(limma_results), "ORFs\n")
cat("ORF abundance:", nrow(orf_data), "ORFs,", ncol(orf_data), "samples\n")
cat("Metadata:", nrow(metadata), "samples\n")

# Check significant ORFs
significant_orfs <- limma_results[limma_results$adj.P.Val < 0.01, ]
cat("Significant ORFs:", nrow(significant_orfs), "\n")

# Check sample distribution by group and timepoint
cat("\nSample distribution:\n")
sample_dist <- metadata %>%
  group_by(Dx.Status, onset_timeline_numeric) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  arrange(onset_timeline_numeric)

print(sample_dist)

# Check if samples have both groups at each timepoint
cat("\nTimepoints with both CELIAC and CONTROL:\n")
both_groups <- sample_dist %>%
  group_by(onset_timeline_numeric) %>%
  summarise(
    n_groups = n(),
    has_celiac = any(Dx.Status == "CELIAC"),
    has_control = any(Dx.Status == "CONTROL"),
    both_groups = has_celiac & has_control,
    .groups = "drop"
  )

print(both_groups)

# Check abundance data sparsity
cat("\nAbundance data sparsity check:\n")
first_10_orfs <- rownames(orf_data)[1:10]
for(orf in first_10_orfs) {
  zero_count <- sum(orf_data[orf, ] == 0, na.rm = TRUE)
  total_count <- ncol(orf_data)
  cat("ORF", orf, ": ", zero_count, "/", total_count, " zeros (", 
      round(zero_count/total_count*100, 1), "%)\n")
}