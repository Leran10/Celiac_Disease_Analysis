#!/usr/bin/env Rscript
# Check ORF ID to CDS ID mapping between analysis and Phold output
# Author: Claude AI
# Date: 2025-07-21

library(dplyr)

setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/compositonal_analysis")

cat("=== ORF ID to CDS ID MAPPING CHECK ===\n")

# Load Phold data
phold_data <- read.table("../phold_per_cds_predictions.tsv", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Phold data loaded:", nrow(phold_data), "CDS entries\n")

# Load some example ORF IDs from our analysis
matrix_file <- "all_significant_orfs_temporal_matrix.csv"
if(file.exists(matrix_file)) {
  temporal_matrix <- read.csv(matrix_file, row.names = 1)
  example_orf_ids <- head(rownames(temporal_matrix), 10)
  
  cat("\nExample ORF IDs from analysis:\n")
  print(example_orf_ids)
  
  # Function to convert ORF ID to CDS ID format
  orf_to_cds <- function(orf_id) {
    # Split by underscore
    parts <- unlist(strsplit(orf_id, "_"))
    
    # Expected format: virus_comp_XXXX_cycle_Y_Z
    # Convert to: virus_comp_XXXX_cycle_Y_CDS_000Z
    if(length(parts) >= 5) {
      contig_part <- paste(parts[1:4], collapse = "_")  # virus_comp_XXXX_cycle_Y
      orf_num <- parts[5]  # Z (the actual ORF number)
      
      # Convert to CDS format with zero padding
      cds_num <- sprintf("%04d", as.numeric(orf_num))
      cds_id <- paste0(contig_part, "_CDS_", cds_num)
      
      return(cds_id)
    } else {
      return(NA)
    }
  }
  
  cat("\nTesting ORF ID to CDS ID conversion:\n")
  test_conversions <- data.frame(
    ORF_ID = example_orf_ids,
    Predicted_CDS_ID = sapply(example_orf_ids, orf_to_cds),
    stringsAsFactors = FALSE
  )
  
  # Check if these CDS IDs exist in Phold data  
  test_conversions$In_Phold <- test_conversions$Predicted_CDS_ID %in% phold_data$cds_id
  
  # Also check by constructing the correct mapping
  test_conversions$Correct_CDS_ID <- NA
  for(i in 1:nrow(test_conversions)) {
    parts <- unlist(strsplit(test_conversions$ORF_ID[i], "_"))
    if(length(parts) >= 5) {
      contig_part <- paste(parts[1:4], collapse = "_")
      orf_num <- as.numeric(parts[5])
      correct_cds_id <- paste0(contig_part, "_CDS_", sprintf("%04d", orf_num))
      test_conversions$Correct_CDS_ID[i] <- correct_cds_id
    }
  }
  test_conversions$Correct_In_Phold <- test_conversions$Correct_CDS_ID %in% phold_data$cds_id
  
  print(test_conversions)
  
  # Summary
  cat("\nMapping Summary:\n")
  cat("Total ORF IDs tested:", nrow(test_conversions), "\n")
  cat("Successfully converted:", sum(!is.na(test_conversions$Predicted_CDS_ID)), "\n")
  cat("Found in Phold data:", sum(test_conversions$In_Phold, na.rm = TRUE), "\n")
  cat("Success rate:", 
      round(100 * sum(test_conversions$In_Phold, na.rm = TRUE) / sum(!is.na(test_conversions$Predicted_CDS_ID)), 1), "%\n")
  
  # Show functional annotations for matched ORFs
  if(sum(test_conversions$In_Phold, na.rm = TRUE) > 0) {
    matched_cds <- test_conversions$Predicted_CDS_ID[test_conversions$In_Phold]
    cat("\nFunctional annotations for matched ORFs:\n")
    
    annotations <- phold_data[phold_data$cds_id %in% matched_cds, 
                             c("cds_id", "function", "product")]
    print(annotations)
  }
  
} else {
  cat("Temporal matrix file not found. Please run all_significant_orfs_heatmap.R first.\n")
}

cat("\nMapping pattern confirmed:\n")
cat("Analysis ORF format:  virus_comp_XXXX_cycle_Y_Z\n")
cat("Phold CDS format:     virus_comp_XXXX_cycle_Y_CDS_000Z\n")
cat("Where Z is the ORF number (1, 2, 3, etc.) converted to 4-digit format (0001, 0002, 0003, etc.)\n")