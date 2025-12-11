library(dplyr)

## Load data
count_mat <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_contig_PA.table_0.75.csv",
                      row.names = 1, check.names = FALSE)
meta <- read.csv("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Fisher/US/US_vertebrate_metadata_0.75.csv",
                 row.names = 1, check.names = FALSE)

## Ensure metadata matches count table samples
meta <- meta[colnames(count_mat), , drop = FALSE]
meta$SampleID <- rownames(meta)

## Check one contig
test_contig <- "contig_3357"

# Extract presence/absence
contig_pa <- as.integer(count_mat[test_contig, ] > 0)

# Create data frame
plot_data <- data.frame(
  SampleID = colnames(count_mat),
  Presence = contig_pa
)

# Merge with metadata
plot_data <- plot_data %>%
  left_join(meta %>% select(SampleID, onset_timeline_combined, Dx.Status),
            by = "SampleID")

# Show first few rows
print("First few rows of data:")
print(head(plot_data, 20))

# Aggregate by timepoint and Dx.Status
plot_data_agg <- plot_data %>%
  group_by(onset_timeline_combined, Dx.Status) %>%
  summarise(
    total_PA = sum(Presence, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

print("\nAggregated data by timepoint and disease status:")
print(plot_data_agg)

# Check samples per timepoint
samples_per_timepoint <- meta %>%
  group_by(onset_timeline_combined, Dx.Status) %>%
  summarise(n_samples = n(), .groups = "drop")

print("\nNumber of samples per timepoint:")
print(samples_per_timepoint)
