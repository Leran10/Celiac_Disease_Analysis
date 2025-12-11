library(dplyr)

## Missing contigs
missing_contigs <- c("contig_3669", "contig_43336", "contig_8153")

## Search in broader directories
search_dirs <- c(
  "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/",
  "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/",
  "~/Handley Lab Dropbox/16S/Celiac/"
)

cat("=== Searching for CSV files containing contig data ===\n\n")

for (dir in search_dirs) {
  if (dir.exists(dir)) {
    cat("Searching in:", dir, "\n")
    files <- list.files(dir, pattern = ".*contig.*\\.csv$|.*PA.*\\.csv$",
                        full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

    if (length(files) > 0) {
      cat("Found", length(files), "files\n")
      for (file in files) {
        # Skip if it's the 0.75 file we already know about
        if (grepl("0.75", file)) next

        # Check file size - skip very small files
        file_size <- file.info(file)$size
        if (file_size < 1000) next

        cat("\n  Checking:", basename(file), "\n")
        cat("  Full path:", file, "\n")

        tryCatch({
          temp_data <- read.csv(file, row.names = 1, check.names = FALSE, nrows = 5)
          cat("  Dimensions: ", nrow(temp_data), "rows (showing first 5)\n")

          # Now read full file to check for missing contigs
          temp_data_full <- read.csv(file, row.names = 1, check.names = FALSE)

          found_missing <- FALSE
          for (contig in missing_contigs) {
            if (contig %in% rownames(temp_data_full)) {
              if (!found_missing) {
                cat("  *** FOUND MISSING CONTIGS IN THIS FILE ***\n")
                found_missing <- TRUE
              }
              contig_pa <- as.integer(temp_data_full[contig, ] > 0)
              total_present <- sum(contig_pa, na.rm = TRUE)
              pct_present <- (total_present / length(contig_pa)) * 100
              cat("    ", contig, ": Present in", total_present, "samples (",
                  round(pct_present, 1), "%)\n")
            }
          }
        }, error = function(e) {
          cat("  Could not read file:", e$message, "\n")
        })
      }
    }
    cat("\n")
  }
}

cat("\n=== Summary ===\n")
cat("The three missing contigs from the 0.75 filtered file are:\n")
for (contig in missing_contigs) {
  cat("  -", contig, "\n")
}
cat("\nThey were likely filtered out due to the '0.75' threshold.\n")
cat("To understand what this threshold means, you may need to check your\n")
cat("data processing pipeline or the code that generated these files.\n")
