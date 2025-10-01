#!/usr/bin/env Rscript
# Debug the ordering issue for function plot

setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

# Load data (simplified version)
limma_results <- read.csv("total_limma_model_res.csv", row.names = 1)
limma_results$orf_id <- rownames(limma_results)

phold_data <- read.delim("phold_per_cds_predictions.tsv", stringsAsFactors = FALSE)
names(phold_data)[names(phold_data) == "function."] <- "func_category"

merged_data <- merge(limma_results, phold_data, 
                     by.x = "orf_id", by.y = "cds_id", 
                     all.x = TRUE)

merged_data$func_category <- ifelse(is.na(merged_data$func_category) | merged_data$func_category == "", 
                                   "no annotation", merged_data$func_category)

sig_data <- merged_data[merged_data$adj.P.Val < 0.01, ]

# Get top functions for coloring
function_counts <- table(sig_data$func_category)
top_functions <- names(sort(function_counts, decreasing = TRUE))[1:8]

sig_data$function_category <- ifelse(sig_data$func_category %in% top_functions, 
                                    sig_data$func_category, "Other/Unknown")

# Define colors
function_colors <- c(
  "head and packaging" = "#E31A1C",
  "DNA, RNA and nucleotide metabolism" = "#1F78B4", 
  "tail" = "#33A02C",
  "transcription regulation" = "#FF7F00",
  "lysis" = "#6A3D9A",
  "connector" = "#A6CEE3",
  "integration and excision" = "#FB9A99",
  "other" = "#FDBF6F",
  "Other/Unknown" = "grey70"
)

# Debug the ordering
gray_function_categories <- names(function_colors)[grepl("grey", function_colors)]
cat("Gray categories detected:", gray_function_categories, "\n")

# Check category distribution
cat("\nFunction category counts:\n")
print(table(sig_data$function_category))

# Check which points are gray vs colored
is_gray <- sig_data$function_category %in% gray_function_categories
cat("\nGray points:", sum(is_gray), "\n")
cat("Colored points:", sum(!is_gray), "\n")

# Test ordering
sig_data_ordered <- sig_data[order(sig_data$function_category %in% gray_function_categories, decreasing = TRUE), ]

cat("\nFirst 10 categories in ordered data:\n")
print(head(sig_data_ordered$function_category, 10))

cat("\nLast 10 categories in ordered data:\n")
print(tail(sig_data_ordered$function_category, 10))

cat("\nFirst 10 should be gray, last 10 should be colored\n")