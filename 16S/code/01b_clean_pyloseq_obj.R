#!/usr/bin/env Rscript


# tell what samples we have been excluding..
# Section: load libraries... ----
library("phyloseq")
library("tidyverse")

# Section: set up paths & source utility functions... ----
p_dir_wd <- Sys.getenv("PROJ_CELIAC_PATH",
                       "/Users/dhohaabid/Documents/proj_celiac")
p_dir_out <- file.path(p_dir_wd, "output")
# input files
p_f_meta <- file.path(p_dir_wd, "data", "metadata", "metadata_20250908.RDS")
p_f_ps <- file.path(p_dir_out, "ps_gg2.RDS")
# output files
p_f_ps_clean_merged <- file.path(p_dir_out, "ps_gg2_cleaned_merged.RDS")
p_f_ps_clean_italy <- file.path(p_dir_out, "ps_gg2_cleaned_italy.RDS")
p_f_ps_clean_usa <- file.path(p_dir_out, "ps_gg2_cleaned_usa.RDS")

source(file.path(p_dir_wd, "code", "utils",
                 "identify_paired_unpaired_samples.R"))


# Section: read input files... ----
# Read phyloseq object from input directory
ps <- readRDS(p_f_ps)
cat("Reading metadata file:", p_f_ps, "\n")
cat("Samples:", nsamples(ps), "\n")
cat("Taxa:", ntaxa(ps), "\n")
# Read metadata from RDS file
df_meta <- readRDS(p_f_meta)
cat("Reading metadata file:", p_f_meta, "\n")
cat("Samples:", nrow(df_meta), "\n")
cat("Variables:", ncol(df_meta), "\n\n")


# Section: determine common sample ids... ----
# and check whether all samples belongs to patients 
# that are matched with a control id
common_sample_ids <- intersect(sample_names(ps), rownames(df_meta))
cat("Common sample ids:", length(common_sample_ids), "\n")
# check whether all samples belongs to a patient with a control id
df_meta <- df_meta[common_sample_ids, , drop = FALSE]
l_check_paired <- identify_paired_unpaired_samples(df_meta = df_meta)
if (isTRUE(l_check_paired$flag_paired)) {
  cat("All common samples are properly paired!\n")
} else {
  cat("Some common samples are not paired. Details:\n")
  if (length(l_check_paired$v_patient_ids_unpaired_controls) > 0) {
    cat("Controls missing pair_id:",
        paste(l_check_paired$v_patient_ids_unpaired_controls, collapse = ", "), "\n")}
  if (length(l_check_paired$v_patient_ids_unpaired_cases) > 0) {
    cat("Cases missing controls:",
        paste(l_check_paired$v_patient_ids_unpaired_cases, collapse = ", "), "\n")}
  # Remove unpaired samples
    cat("    Removing unpaired samples...\n")
    v_unpaired_patient_ids <- c(
      l_check_paired$v_patient_ids_unpaired_controls,
      l_check_paired$v_patient_ids_unpaired_cases
    )

    if (length(v_unpaired_patient_ids) > 0) {
      cat(
        "    Removing", length(v_unpaired_patient_ids), "unpaired samples: ",
        paste(v_unpaired_patient_ids, collapse = ", "), "\n"
      )
      df_meta <- df_meta %>%
        filter(!patient_id %in% v_unpaired_patient_ids)
      cat("    Number of rows after removing unpaired samples: ", nrow(df_meta), "\n")
      common_sample_ids <- rownames(df_meta)
    }
}


# Section: attach metadata to phyloseq object... ----
cat("\n=== Section: attach metadata to phyloseq object... ===\n")
cat("Pruning phyloseq to only include common sample ids...\n")
ps <- prune_samples(common_sample_ids, ps)
cat("Attaching metadata to phyloseq object...\n")
sample_data(ps) <- sample_data(df_meta)


# Section: Remove samples with low read counts from the phyloseq object... ----
cat("\n=== Section: Remove samples with low read counts from the phyloseq object... ===\n")
ps_clean_samples <- prune_samples(sample_sums(ps) > 1000, ps)
cat("# of samples after removing samples with low read counts:", nsamples(ps_clean_samples), "\n")
cat("# of taxa after removing samples with low read counts:", ntaxa(ps_clean_samples), "\n")
# Get filtered metedata from phyloseq object and make sure that it is
# of class dataframe by removing any sample_data attributes
# that might cause issues with dplyr
df_metadata_fltr <- as.data.frame(sample_data(ps_clean_samples))
attr(df_metadata_fltr, "class") <- "data.frame"
# check whether all samples belongs to a patient with a control id
l_check_paired <- identify_paired_unpaired_samples(df_meta = df_metadata_fltr)
if (isTRUE(l_check_paired$flag_paired)) {
  cat("All common samples are properly paired!\n")
} else {
  cat("Some common samples are not paired. Details:\n")
  if (length(l_check_paired$v_patient_ids_unpaired_controls) > 0) {
    cat("Controls missing pair_id:",
      paste(l_check_paired$v_patient_ids_unpaired_controls, collapse = ", "), "\n")}
  if (length(l_check_paired$v_patient_ids_unpaired_cases) > 0) {
    cat("Cases missing controls:",
      paste(l_check_paired$v_patient_ids_unpaired_cases, collapse = ", "), "\n")}
}


# Section: remove taxa with no read counts...----
cat("\n=== Section: remove taxa with no read counts... ===\n")
ps_clean_taxa <- prune_taxa(taxa_sums(ps_clean_samples) > 0, ps_clean_samples)
cat("# of samples after removing taxa with no reads:", nsamples(ps_clean_taxa), "\n")
cat("# of taxa after removing taxa with no reads:", ntaxa(ps_clean_taxa), "\n")


# Section: save cleaned phyloseq objects...----
cat("\n=== Section: save cleaned phyloseq objects to output directory ===\n")
cat("Save merged cohort: to ", p_f_ps_clean_merged, "\n")
saveRDS(ps_clean_taxa, file = p_f_ps_clean_merged)
cat("Save Italy cohort: to ", p_f_ps_clean_italy, "\n")
ps_clean_taxa_italy <- subset_samples(ps_clean_taxa, country == "ITALY")
saveRDS(ps_clean_taxa_italy, file = p_f_ps_clean_italy)
cat("Save USA cohort: to ", p_f_ps_clean_usa, "\n")
ps_clean_taxa_usa <- subset_samples(ps_clean_taxa, country == "USA")
saveRDS(ps_clean_taxa_usa, file = p_f_ps_clean_usa)