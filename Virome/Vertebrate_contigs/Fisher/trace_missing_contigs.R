library(dplyr)
library(tibble)

## Missing contigs
missing_contigs <- c("contig_3669", "contig_43336", "contig_8153")

cat("=== Tracing where the three contigs were lost ===\n\n")

## Step 1: Check original RDS file
cat("STEP 1: Checking original 75% coverage filtered file\n")
original_file <- "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/wide_rpkm_genes_75Cov.rds"

if (file.exists(original_file)) {
  US_vertebrate_contig_processed <- readRDS(original_file) %>%
    column_to_rownames("ORFID")

  cat("Original dimensions:", nrow(US_vertebrate_contig_processed), "contigs x",
      ncol(US_vertebrate_contig_processed), "samples\n\n")

  for (contig in missing_contigs) {
    if (contig %in% rownames(US_vertebrate_contig_processed)) {
      contig_values <- US_vertebrate_contig_processed[contig, ]
      non_zero_count <- sum(contig_values > 0, na.rm = TRUE)
      total_count <- sum(contig_values, na.rm = TRUE)

      cat(contig, ":\n")
      cat("  - PRESENT in original file\n")
      cat("  - Non-zero in", non_zero_count, "out of", ncol(US_vertebrate_contig_processed), "samples\n")
      cat("  - Total abundance:", total_count, "\n")

      if (non_zero_count > 0) {
        cat("  - Samples where present:\n")
        present_samples <- names(contig_values[contig_values > 0])
        for (sample in present_samples) {
          cat("    *", sample, "=", contig_values[sample], "\n")
        }
      }
      cat("\n")
    } else {
      cat(contig, ": NOT in original 75% coverage file\n\n")
    }
  }

  ## Step 2: Check which samples were removed
  cat("\n=== STEP 2: Checking which samples were removed ===\n")

  # Find samples with all zeros
  sample_all_0 <- names(colSums(US_vertebrate_contig_processed)[colSums(US_vertebrate_contig_processed) == 0])
  cat("Samples with all zeros (removed in step 1):", length(sample_all_0), "\n")

  # Remove zero samples
  US_vertebrate_contig_processed <- US_vertebrate_contig_processed %>%
    select(-all_of(sample_all_0))

  # Clean column names (same as original script)
  US_vertebrate_contig_count <- US_vertebrate_contig_processed
  colnames(US_vertebrate_contig_count) <- gsub("NovaSeq_N966_|NovaSeq_N978_|NovaSeq_N979_|NovaSeq_N980_|NovaSeq_N981_|NovaSeq_N982_|NovaSeq_N983_|_stats","",colnames(US_vertebrate_contig_count))
  colnames(US_vertebrate_contig_count) <- stringr::str_replace_all(colnames(US_vertebrate_contig_count),"Leonard_Human_stool_Celiac","Celiac_Leonard_Stool")
  colnames(US_vertebrate_contig_count) <- sapply(strsplit(colnames(US_vertebrate_contig_count),"_Stool_|_Stool_Celiac_"),"[",2)

  cat("\nAfter removing all-zero samples:", ncol(US_vertebrate_contig_count), "samples\n")

  cat("\nChecking missing contigs after sample removal:\n")
  for (contig in missing_contigs) {
    if (contig %in% rownames(US_vertebrate_contig_count)) {
      non_zero_count <- sum(US_vertebrate_contig_count[contig, ] > 0, na.rm = TRUE)
      cat(contig, ": Non-zero in", non_zero_count, "samples\n")
    }
  }

} else {
  cat("ERROR: Could not find original file at:", original_file, "\n")
  cat("Please check the file path.\n")
}

cat("\n=== CONCLUSION ===\n")
cat("The three missing contigs were lost because:\n")
cat("1. They were either not in the original 75% coverage filtered file, OR\n")
cat("2. They became all zeros after removing samples with 'Unknown' metadata\n")
cat("   and were removed in the 'all_0_contigs' filtering step.\n")
