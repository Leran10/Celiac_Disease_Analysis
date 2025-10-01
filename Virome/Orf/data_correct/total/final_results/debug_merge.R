#!/usr/bin/env Rscript
# Debug the merge issue

setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

# Load data
limma_results <- read.csv("total_limma_model_res.csv", row.names = 1)
limma_results$orf_id <- rownames(limma_results)

phold_data <- read.delim("phold_per_cds_predictions.tsv", stringsAsFactors = FALSE)

# Check data
cat("Limma results:\n")
cat("  Rows:", nrow(limma_results), "\n")
cat("  First 5 IDs:", head(limma_results$orf_id, 5), "\n")

cat("\nPhold data:\n")
cat("  Rows:", nrow(phold_data), "\n")
cat("  First 5 IDs:", head(phold_data$cds_id, 5), "\n")
cat("  Column names:", colnames(phold_data), "\n")

# Test merge
merged_test <- merge(limma_results, phold_data, 
                     by.x = "orf_id", by.y = "cds_id", 
                     all.x = TRUE)

cat("\nMerged data:\n")
cat("  Rows:", nrow(merged_test), "\n")
cat("  Columns:", ncol(merged_test), "\n")
cat("  Has function column:", "function" %in% colnames(merged_test), "\n")
cat("  Function column summary:\n")
if("function" %in% colnames(merged_test)) {
  print(table(merged_test$`function`, useNA = "always"))
}