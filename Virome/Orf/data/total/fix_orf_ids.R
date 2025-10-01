#!/usr/bin/env Rscript
# Fix the ORF ID column to ensure uniqueness
# Author: Claude AI
# Date: 2025-07-21

library(dplyr)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total")

cat("=== FIXING ORF ID COLUMN ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# Load the table with problematic ORF IDs
cat("1. Loading table with ORF IDs...\n")
phold_data <- read.table("phold_per_cds_predictions_with_orf_ids.tsv", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                        fill = TRUE, comment.char = "")

cat("Original data:", nrow(phold_data), "rows\n")

# Check uniqueness before fix
cat("CDS IDs unique:", length(unique(phold_data$cds_id)) == nrow(phold_data), "\n")
cat("ORF IDs unique:", length(unique(phold_data$orf_id)) == nrow(phold_data), "\n")
cat("Unique ORF IDs:", length(unique(phold_data$orf_id)), "out of", nrow(phold_data), "\n")

# Fixed conversion function
cds_to_orf_fixed <- function(cds_id) {
  # Input format: virus_comp_XXXX_cycle_Y_CDS_ZZZZ
  # Output format: virus_comp_XXXX_cycle_Y_Z
  
  # Split by underscore
  parts <- unlist(strsplit(cds_id, "_"))
  
  if(length(parts) >= 7 && parts[6] == "CDS") {
    # Extract components correctly
    contig_part <- paste(parts[1:5], collapse = "_")  # virus_comp_XXXX_cycle_Y
    cds_number <- parts[7]  # ZZZZ (like 0001, 0002, etc.)
    
    # Convert CDS number to integer (removes leading zeros)
    orf_number <- as.numeric(cds_number)
    
    # Build complete ORF ID with contig part
    orf_id <- paste0(contig_part, "_", orf_number)
    
    return(orf_id)
  } else {
    # If format doesn't match, return the original
    return(cds_id)
  }
}

# Apply the fixed conversion
cat("2. Applying fixed conversion...\n")
phold_data$orf_id <- sapply(phold_data$cds_id, cds_to_orf_fixed)

# Check results
successful_conversions <- sum(!is.na(phold_data$orf_id))
cat("Successful conversions:", successful_conversions, "out of", nrow(phold_data), "\n")

# Check uniqueness after fix
cat("ORF IDs unique after fix:", length(unique(phold_data$orf_id)) == nrow(phold_data), "\n")
cat("Unique ORF IDs after fix:", length(unique(phold_data$orf_id)), "out of", nrow(phold_data), "\n")

# Show some examples
cat("\nExample conversions:\n")
examples <- head(phold_data[, c("cds_id", "orf_id")], 10)
print(examples)

# Save the corrected table
output_file <- "phold_per_cds_predictions_with_orf_ids_fixed.tsv"
cat("3. Saving corrected table to", output_file, "...\n")

write.table(phold_data, output_file, 
           sep = "\t", row.names = FALSE, quote = FALSE, na = "")

cat("Fixed table saved successfully!\n")
cat("Analysis completed at:", format(Sys.time()), "\n")