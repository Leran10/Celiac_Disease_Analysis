# Section: Load required libraries... ----
library(tidyverse)
library(readr)


# Section: Set working directory and paths... ----
p_wd <- Sys.getenv("PROJ_CELIAC_PATH", "/Users/dhohaabid/Documents/proj_celiac")
p_data <- file.path(p_wd, "data")
p_output <- file.path(p_wd, "output")


# Section: Set paths for input and output files... ----
# Read the metadata file
p_f_metadata <- file.path(
  p_data, "metadata", "orig",
  "Updated_Metadata_with_Onset_Timeline.csv"
)
p_f_output <- file.path(p_data, "metadata", "metadata_20250908.RDS")



# Section: Read the metadata file... ----
df_metadata <- read_csv(p_f_metadata, show_col_types = FALSE, col_names = TRUE)
cat("Original metadata dimensions:", dim(df_metadata), "\n")
cat("Original column names:\n")
print(colnames(df_metadata))

# Section: Rename columns... ----
column_mapping <- c(
  "sample_id" = "...1",
  "run_number" = "RunNumber",
  "sample_name" = "Sample_full_name",
  "patient_id" = "patientID",
  "sample_external_id" = "Sample.Name.External.ID",
  "classification_original" = "Classification",
  "classification_updated" = "Updated.Classification",
  "country" = "Country",
  "sex" = "Sex",
  "delivery_mode" = "Delivery.Mode",
  "hla_genotype" = "HLA",
  "hla_risk_category" = "HLA.Category",
  "age_at_onset" = "Onset.Dx",
  "disease_status" = "Dx.Status",
  "feeding_type_0_6" = "Feeding.Type..0.6.months.",
  "feeding_type_7_12" = "Feeding.Type..7.12.months.",
  "age_at_solid_introduction" = "Age.at.Solid.Introduction..months.",
  "age_at_gluten_introduction" = "Age.at.Gluten.Introduction..months.",
  "pair_id" = "Pair",
  "celiac_group" = "celiacs_group",
  "age" = "month",
  "disease_onset_details" = "disease_onset_details",
  "status_timepoint" = "Status_timepoint",
  "feeding_type_first_year" = "feeding_first_year",
  "onset_timeline" = "onset_timeline",
  "onset_timeline_combined" = "onset_timeline_combined",
  "time_to_onset" = "onset_timeline_numeric"
)
# Rename columns
df_metadata_cleaned <- df_metadata %>%
  rename(!!!column_mapping)
cat("Renamed column names:\n")
print(colnames(df_metadata_cleaned))

# Section: Check for unique sample IDs... ----
cat("Total samples:", nrow(df_metadata_cleaned), "\n")
cat("Unique sample IDs:", length(unique(df_metadata_cleaned$sample_id)), "\n")
cat(
  "Are all sample IDs unique?",
  length(unique(df_metadata_cleaned$sample_id)) == nrow(df_metadata_cleaned),
  "\n"
)
# If there are duplicates, show them
duplicated_ids <-
  df_metadata_cleaned$sample_id[duplicated(df_metadata_cleaned$sample_id)]
if (length(duplicated_ids) > 0) {
  cat("DUPLICATED SAMPLE IDs FOUND:\n")
  print(duplicated_ids)
  cat("Rows with duplicated sample IDs:\n")
  print(df_metadata_cleaned[
    df_metadata_cleaned$sample_id %in% duplicated_ids,
    c("sample_id", "patient_id", "age")
  ])
} else {
  cat("✓ All sample IDs are unique!\n")
}

# Section: Clean special characters in hla_risk_category column... ----
cat("Cleaning special characters in hla_risk_category column...\n")
df_metadata_cleaned$hla_genotype <-
  gsub("[/., ]", "_", df_metadata_cleaned$hla_genotype)

# Section: Remove rows with NA, unknown, or empty values... ----
# There are predictors that we are interested in, so we need to remove the rows
# with NA, unknown, or empty values for these predictors
vars_to_check <- c(
  "disease_status", "time_to_onset", "sex",
  "delivery_mode", "age_at_gluten_introduction",
  "age_at_solid_introduction",
  "hla_risk_category", "feeding_type_first_year",
  "patient_id"
)
# Create logical vector for complete cases (no NAs)
complete_cases <- complete.cases(df_metadata_cleaned[, vars_to_check])
# Create logical vector for no empty/unknown values
no_problematic_values <-
  apply(df_metadata_cleaned[, vars_to_check], 1, function(row) {
    all(!is.na(row) &
      row != "" &
      !grepl("^\\s*$", row) & # not just whitespace
      !grepl("^unknown$", row, ignore.case = TRUE) &
      !grepl("^n/a$", row, ignore.case = TRUE))
  })
# Combine both conditions
df_metadata_cleaned <-
  df_metadata_cleaned[complete_cases & no_problematic_values, ]
cat("Number of rows before removing missing values:", nrow(df_metadata), "\n")
cat(
  "Number of rows after removing missing values:",
  nrow(df_metadata_cleaned), "\n"
)

# Section: Display first few rows to verify... ----
cat("First 3 rows of renamed metadata:\n")
print(head(df_metadata_cleaned, 3))


# Section: Summary statistics... ----
cat("Total samples:", nrow(df_metadata_cleaned), "\n")
cat("Total variables:", ncol(df_metadata_cleaned), "\n")

# Section: Also save as TSV file... ----
p_f_output_tsv <- file.path(p_data, "metadata", "metadata_20250908.tsv")
cat("Save cleaned metadata as TSV in: ", p_f_output_tsv, "\n")
write.table(df_metadata_cleaned, p_f_output_tsv,
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)
# Set sample_id as row names
df_metadata_cleaned <- df_metadata_cleaned %>%
  column_to_rownames("sample_id")
# Save the cleaned metadata
cat("Save cleaned metadata file in: ", p_f_output, "\n")
saveRDS(df_metadata_cleaned, p_f_output)