library(dplyr)

## Load the filtered data
count_mat <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_contig_PA.table_0.75.csv",
                      row.names = 1, check.names = FALSE)

## Contigs of interest
Astroviridae_contigs <- c("contig_3357","contig_34440","contig_3567","contig_3669",
                          "contig_43336","contig_4647","contig_49047","contig_49065","contig_8153")

## Check which contigs are present in the filtered table
cat("=== Checking which contigs are in the filtered table (0.75) ===\n")
for (contig in Astroviridae_contigs) {
  in_table <- contig %in% rownames(count_mat)
  cat(contig, ":", ifelse(in_table, "PRESENT", "MISSING"), "\n")
}

## Try to find unfiltered or less filtered versions
cat("\n=== Looking for other data files in the Fisher/US directory ===\n")
files <- list.files("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/",
                    pattern = "*.csv", full.names = TRUE)
cat("Found files:\n")
print(files)

## Check if there's an unfiltered version
possible_files <- c(
  "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_contig_PA.table.csv",
  "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_contig_PA.csv",
  "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/vertebrate_contig_PA.csv"
)

for (file in possible_files) {
  if (file.exists(file)) {
    cat("\n=== Checking file:", file, "===\n")
    temp_data <- read.csv(file, row.names = 1, check.names = FALSE)
    cat("Number of contigs:", nrow(temp_data), "\n")
    cat("Number of samples:", ncol(temp_data), "\n")

    cat("\nChecking Astroviridae contigs in this file:\n")
    for (contig in Astroviridae_contigs) {
      if (contig %in% rownames(temp_data)) {
        contig_pa <- as.integer(temp_data[contig, ] > 0)
        total_present <- sum(contig_pa, na.rm = TRUE)
        pct_present <- (total_present / length(contig_pa)) * 100
        cat(contig, ": Present in", total_present, "samples (", round(pct_present, 1), "%)\n")
      } else {
        cat(contig, ": NOT IN THIS FILE\n")
      }
    }
  }
}

## Check parent directory for original files
cat("\n=== Checking parent directory ===\n")
parent_files <- list.files("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/",
                           pattern = "*.csv", full.names = TRUE, recursive = FALSE)
cat("Files in parent directory:\n")
print(parent_files)
