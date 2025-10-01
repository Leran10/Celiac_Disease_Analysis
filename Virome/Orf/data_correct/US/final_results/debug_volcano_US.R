#!/usr/bin/env Rscript
# Debug volcano plots for US cohort
library(dplyr)

setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/US/final_results")

# Load data
limma_results <- read.csv("../US.limma_model_res.csv")
phold_data <- read.csv("../../total/final_results/phold_per_cds_predictions.tsv", sep = "\t")

cat("Limma columns:", paste(colnames(limma_results), collapse = ", "), "\n")
cat("Phold columns:", paste(colnames(phold_data), collapse = ", "), "\n")

# Filter significant ORFs
sig_data <- limma_results[limma_results$adj.P.Val < 0.01, ]
cat("Significant ORFs:", nrow(sig_data), "\n")

# Merge
merged_data <- merge(sig_data, phold_data, by.x = "X", by.y = "cds_id", all.x = TRUE)
cat("Merged rows:", nrow(merged_data), "\n")
cat("Merged columns:", paste(colnames(merged_data), collapse = ", "), "\n")

# Check function column
if("function" %in% colnames(merged_data)) {
  cat("Function column exists\n")
  func_vals <- merged_data[["function"]]
  cat("Function column class:", class(func_vals), "\n")
  cat("Function column length:", length(func_vals), "\n")
  cat("First 5 function values:", paste(head(func_vals, 5), collapse = ", "), "\n")
  cat("Number of non-NA function values:", sum(!is.na(func_vals)), "\n")
} else {
  cat("Function column missing!\n")
}