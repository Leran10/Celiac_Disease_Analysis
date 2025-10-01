#!/usr/bin/env Rscript
# Add ORF ID column to Phold predictions table
# Converts CDS format to ORF format used in analysis
# Author: Claude AI
# Date: 2025-07-21

library(dplyr)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total")

cat("=== ADDING ORF ID COLUMN TO PHOLD TABLE ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# Load the Phold predictions
cat("1. Loading Phold data...\n")
phold_data <- read.table("phold_per_cds_predictions.tsv", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Original Phold data:", nrow(phold_data), "rows,", ncol(phold_data), "columns\n")

# Function to convert CDS ID to ORF ID format
cds_to_orf <- function(cds_id) {
  # Input format: virus_comp_XXXX_cycle_Y_CDS_ZZZZ
  # Output format: virus_comp_XXXX_cycle_Y_Z
  
  # Split by underscore
  parts <- unlist(strsplit(cds_id, "_"))
  
  if(length(parts) >= 7 && parts[6] == "CDS") {
    # Extract components
    contig_part <- paste(parts[1:5], collapse = "_")  # virus_comp_XXXX_cycle_Y
    cds_number <- parts[7]  # ZZZZ (like 0001, 0002, etc.)
    
    # Convert CDS number to integer (removes leading zeros)
    orf_number <- as.numeric(cds_number)
    
    # Construct ORF ID
    orf_id <- paste0(contig_part, "_", orf_number)
    
    return(orf_id)
  } else {
    # If format doesn't match expected pattern, return NA
    return(NA)
  }
}

# Apply conversion to create ORF ID column
cat("2. Converting CDS IDs to ORF ID format...\n")
phold_data$orf_id <- sapply(phold_data$cds_id, cds_to_orf)

# Check conversion results
successful_conversions <- sum(!is.na(phold_data$orf_id))
cat("Successful conversions:", successful_conversions, "out of", nrow(phold_data), 
    "(", round(100 * successful_conversions / nrow(phold_data), 1), "%)\n")

# Show some examples
cat("\nExample conversions:\n")
examples <- head(phold_data[!is.na(phold_data$orf_id), c("cds_id", "orf_id")], 10)
print(examples)

# Move orf_id column to be the third column (after contig_id and cds_id)
cat("3. Reordering columns...\n")
phold_with_orf <- phold_data %>%
  select(contig_id, cds_id, orf_id, everything())

cat("Final table:", nrow(phold_with_orf), "rows,", ncol(phold_with_orf), "columns\n")

# Save the enhanced table
output_file <- "phold_per_cds_predictions_with_orf_ids.tsv"
cat("4. Saving enhanced table to", output_file, "...\n")

write.table(phold_with_orf, output_file, 
           sep = "\t", row.names = FALSE, quote = FALSE, na = "")

cat("Enhanced Phold table saved successfully!\n")

# Create summary statistics
cat("\nSummary Statistics:\n")
cat("- Total CDS entries:", nrow(phold_with_orf), "\n")
cat("- Successfully converted to ORF format:", sum(!is.na(phold_with_orf$orf_id)), "\n")
cat("- Unique contigs:", length(unique(phold_with_orf$contig_id)), "\n")
functional_annotations <- sum(!is.na(phold_with_orf$function) & phold_with_orf$function != "unknown function")
cat("- Functional annotations available:", functional_annotations, "\n")

# Show column names for verification
cat("\nColumn names in enhanced table:\n")
print(colnames(phold_with_orf))

cat("\nAnalysis completed at:", format(Sys.time()), "\n")