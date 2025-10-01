#!/usr/bin/env Rscript
# Test the CDS to ORF conversion function
# Author: Claude AI

# Test function
cds_to_orf <- function(cds_id) {
  # Input format: virus_comp_XXXX_cycle_Y_CDS_ZZZZ
  # Output format: virus_comp_XXXX_cycle_Y_Z
  
  # Split by underscore
  parts <- unlist(strsplit(cds_id, "_"))
  
  cat("Parsing:", cds_id, "\n")
  cat("Parts:", paste(parts, collapse=" | "), "\n")
  cat("Length:", length(parts), "\n")
  cat("Part 6:", if(length(parts) >= 6) parts[6] else "NA", "\n")
  
  if(length(parts) >= 7 && parts[6] == "CDS") {
    # Extract components
    contig_part <- paste(parts[1:5], collapse = "_")  # virus_comp_XXXX_cycle_Y  
    cds_number <- parts[7]  # ZZZZ (like 0001, 0002, etc.)
    
    # Convert CDS number to integer (removes leading zeros)
    orf_number <- as.numeric(cds_number)
    
    # Build ORF ID
    orf_id <- paste0(contig_part, "_", orf_number)
    
    cat("Result:", orf_id, "\n\n")
    return(orf_id)
  } else {
    cat("No match\n\n")
    return(NA)
  }
}

# Test with real examples
test_cases <- c(
  "virus_comp_2179_cycle_1_CDS_0001",
  "virus_comp_315_cycle_1_CDS_0001", 
  "virus_comp_2572_cycle_1_CDS_0004"
)

for(test_case in test_cases) {
  result <- cds_to_orf(test_case)
}