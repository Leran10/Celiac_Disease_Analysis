#!/usr/bin/env Rscript

# This script contains utility functions for the celiac project
# Author: Dhoha Abid
# Date: 09/18/2025


# Section: load libraries... ----
library(tidyverse)
library(dplyr)
library(readr)
library(tibble)


# Section: Function to identify paired/unpaired patients across all data
identify_paired_unpaired_samples <- function(df_meta) {
  # Separate df_controls and df_cases
  df_controls <- df_meta %>% filter(disease_status == "CONTROL")
  df_cases <- df_meta %>% filter(disease_status == "CELIAC")
  # Find control samples with missing/invalid pair_id
  df_controls_missing_pair <- df_controls %>%
    filter(is.na(pair_id) | pair_id == "" | pair_id == "NA")
  v_patient_ids_unpaired_controls <-
    df_controls_missing_pair$patient_id %>% unique()
  # Extract pair numbers from pair_id and create control-pair mapping
  df_controls_with_pair <- df_controls %>%
    filter(!is.na(pair_id) & pair_id != "" & pair_id != "NA")
  # Create a data frame with control samples and their pair IDs
  df_control_pair_mapping <- unique(data.frame(
    control_patient_id = df_controls_with_pair$patient_id,
    pair_id = df_controls_with_pair$pair_id,
    stringsAsFactors = FALSE
  ))
  # Find corresponding case samples and loop over all pair_ids to find
  # whether they are in the case samples.
  # Initialize results
  df_patient_ids_paired_details <- tibble(
    patient_id_case = character(0),
    patient_id_control = character(0)
  )
  v_patient_ids_paired <- c()
  for (i in seq_len(nrow(df_control_pair_mapping))) {
    patient_id_control <- df_control_pair_mapping$control_patient_id[i]
    pair_num <- df_control_pair_mapping$pair_id[i]
    # Extract patient id from pair_id
    l_part <- strsplit(pair_num, "_")[[1]]
    patient_id_case <- paste0("0", l_part[1], "_GEMM_", l_part[2])
    # Check whether case_patient_id is a case sample
    df_samples_paired <- df_cases %>% filter(patient_id == patient_id_case)
    if (nrow(df_samples_paired) > 0) {
      df_patient_ids_paired_details <- add_row(df_patient_ids_paired_details,
        patient_id_case = patient_id_case,
        patient_id_control = patient_id_control
      )
      v_patient_ids_paired <- c(v_patient_ids_paired, patient_id_case)
    }
  }
  # Identify unpaired patients
  v_patient_ids_all <- df_cases$patient_id %>% unique()
  v_patient_ids_unpaired_cases <-
    setdiff(v_patient_ids_all, v_patient_ids_paired)

  # Check if there are unpaired samples
  if (length(v_patient_ids_unpaired_cases) > 0 ||
    length(v_patient_ids_unpaired_controls) > 0) {
    return(list(
      flag_paired = FALSE,
      v_patient_ids_unpaired_controls = v_patient_ids_unpaired_controls,
      v_patient_ids_unpaired_cases = v_patient_ids_unpaired_cases,
      df_patient_ids_paired_details = df_patient_ids_paired_details
    ))
  } else {
    return(list(
      flag_paired = TRUE,
      df_patient_ids_paired_details = df_patient_ids_paired_details
    ))
  }
}


# Section: function to add significance to p-values... ----
add_significance <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  }
  if (p_value <= 0.001) {
    return("***")
  }
  if (p_value <= 0.01) {
    return("**")
  }
  if (p_value <= 0.05) {
    return("*")
  }
  return("")
}


# Section: save results in excel file... ----
save_results_in_excel_file <- function(df_results, v_metrics, p_f_excel) {
  # Create workbook
  wb <- openxlsx::createWorkbook()
  # Add a sheet per metric with formatted headers and auto-width columns
  for (metric_name in v_metrics) {
    df_results_metric <- df_results[df_results$metric == metric_name, ]
    openxlsx::addWorksheet(wb, sheetName = metric_name)
    openxlsx::writeData(wb, sheet = metric_name, x = df_results_metric)
    openxlsx::addStyle(
      wb,
      sheet = metric_name,
      style = openxlsx::createStyle(textDecoration = "bold"),
      rows = 1,
      cols = seq_len(ncol(df_results_metric))
    )
    openxlsx::setColWidths(
      wb,
      sheet = metric_name,
      cols = seq_len(ncol(df_results_metric)),
      widths = "auto"
    )
  }
  openxlsx::saveWorkbook(wb, p_f_excel, overwrite = TRUE)
  cat("Excel file saved to:", p_f_excel, "\n")
}
