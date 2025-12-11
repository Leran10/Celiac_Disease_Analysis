library(dplyr)

## Load data
count_mat <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_contig_PA.table_0.75.csv",
                      row.names = 1, check.names = FALSE)
meta <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_metadata_0.75.csv",
                 row.names = 1, check.names = FALSE)

## Contigs of interest
Astroviridae_contigs <- c("contig_3357","contig_34440","contig_3567","contig_3669",
                          "contig_43336","contig_4647","contig_49047","contig_49065","contig_8153")

## Ensure metadata matches count table samples
meta <- meta[colnames(count_mat), , drop = FALSE]
meta$SampleID <- rownames(meta)

## Check all contigs
for (contig in Astroviridae_contigs) {
  # Extract presence/absence
  contig_pa <- as.integer(count_mat[contig, ] > 0)

  # Create data frame
  plot_data <- data.frame(
    SampleID = colnames(count_mat),
    Presence = contig_pa
  )

  # Merge with metadata
  plot_data <- plot_data %>%
    left_join(meta %>% select(SampleID, onset_timeline_combined, Dx.Status),
              by = "SampleID")

  # Aggregate by timepoint and Dx.Status
  plot_data_agg <- plot_data %>%
    group_by(onset_timeline_combined, Dx.Status) %>%
    summarise(
      total_PA = sum(Presence, na.rm = TRUE),
      .groups = "drop"
    )

  # Find max count
  max_count <- max(plot_data_agg$total_PA)

  # Show results
  cat("\n", contig, "- Maximum PA count:", max_count, "\n")

  if (max_count > 1) {
    cat("Top timepoints with presence:\n")
    print(plot_data_agg %>% filter(total_PA > 0) %>% arrange(desc(total_PA)) %>% head())
  }
}

# Overall statistics
cat("\n\n=== Overall Contig Presence Statistics ===\n")
for (contig in Astroviridae_contigs) {
  contig_pa <- as.integer(count_mat[contig, ] > 0)
  total_present <- sum(contig_pa)
  pct_present <- (total_present / length(contig_pa)) * 100
  cat(contig, "- Present in", total_present, "out of", length(contig_pa),
      "samples (", round(pct_present, 1), "%)\n")
}
