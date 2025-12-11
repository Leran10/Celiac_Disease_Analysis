# Debug Italy Abundance Analysis
library(dplyr)
library(tidyr)
library(tibble)
library(edgeR)
library(glmmTMB)
library(purrr)
library(broom.mixed)

# Load the data
Italy_vertebrate_contig_abundance_0.75_0.03 <- read.csv(
  "/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_abundance_0.75_0.03.csv",
  row.names = 1
)

Italy_vertebrate_metadata_0.75_0.03 <- read.csv(
  "/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_metadata_0.75_0.03",
  row.names = 1
)

print("=== Data dimensions ===")
print(paste("Abundance matrix:", nrow(Italy_vertebrate_contig_abundance_0.75_0.03), "contigs x",
            ncol(Italy_vertebrate_contig_abundance_0.75_0.03), "samples"))
print(paste("Metadata:", nrow(Italy_vertebrate_metadata_0.75_0.03), "samples x",
            ncol(Italy_vertebrate_metadata_0.75_0.03), "variables"))

print("\n=== Metadata column names ===")
print(colnames(Italy_vertebrate_metadata_0.75_0.03))

print("\n=== Check factor levels ===")
print("Dx.Status levels:")
print(table(Italy_vertebrate_metadata_0.75_0.03$Dx.Status))
print("onset_timeline_combined levels:")
print(table(Italy_vertebrate_metadata_0.75_0.03$onset_timeline_combined))

# Prepare DGEList
Italy.vertebrate_y <- DGEList(counts = Italy_vertebrate_contig_abundance_0.75_0.03,
                              samples = Italy_vertebrate_metadata_0.75_0.03)
Italy.vertebrate_y <- calcNormFactors(Italy.vertebrate_y, method = "TMMwsp")

# Create long-format data
vertebrate_contig_Italy_abundance_long <- Italy_vertebrate_contig_abundance_0.75_0.03 %>%
  as.data.frame() %>%
  rownames_to_column("vertebrate_contig_id") %>%
  pivot_longer(cols = -vertebrate_contig_id, names_to = "sample_id", values_to = "count") %>%
  left_join(Italy.vertebrate_y$samples %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  dplyr::mutate(offset = log(lib.size * norm.factors))

print("\n=== Long format data structure ===")
print(dim(vertebrate_contig_Italy_abundance_long))
print(head(vertebrate_contig_Italy_abundance_long, 10))

print("\n=== Trying to fit ONE model for the first contig ===")
test_data <- vertebrate_contig_Italy_abundance_long %>%
  filter(vertebrate_contig_id == vertebrate_contig_Italy_abundance_long$vertebrate_contig_id[1])

print(paste("Testing with contig:", test_data$vertebrate_contig_id[1]))
print(paste("Number of observations:", nrow(test_data)))
print("Data summary:")
print(summary(test_data))

# Try fitting the model
test_model <- tryCatch({
  glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. +
          HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
          family = nbinom2, data = test_data)
}, error = function(e) {
  print(paste("ERROR:", e$message))
  return(NULL)
})

if (!is.null(test_model)) {
  print("\n=== Model fitted successfully! ===")
  print(summary(test_model))

  # Try tidying
  tidy_result <- tryCatch({
    tidy(test_model, conf.int = TRUE)
  }, error = function(e) {
    print(paste("Tidy ERROR:", e$message))
    return(NULL)
  })

  if (!is.null(tidy_result)) {
    print("\n=== Tidy results ===")
    print(tidy_result)
  }
} else {
  print("\n=== Model fitting FAILED ===")
}
