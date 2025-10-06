#!/usr/bin/env Rscript


# Section: Main function - Differential abundance analysis using ANCOM-BC... ----
main <- function() {
  cat("Loading phyloseq objects...\n")
  ps_usa <- readRDS(p_f_ps_usa)
  ps_italy <- readRDS(p_f_ps_italy)
  cat("Filtering phyloseq objects by 3% prevalence...\n")
  ps_usa_filtered <- prepare_data(ps_usa, "USA", prev_threshold = 0.03)
  ps_italy_filtered <- prepare_data(ps_italy, "Italy", prev_threshold = 0.03)
  cat("Running ANCOM-BC models...\n")
  l_results_usa <- fit_model(ps_usa_filtered)
  l_results_italy <- fit_model(ps_italy_filtered)
  # Section: save and write results to excel files
  cat("Saving results to Excel files...\n")
  save_l_results(l_results_usa, "USA", p_f_excel_usa)
  save_l_results(l_results_italy, "Italy", p_f_excel_italy)
  cat("Plotting differential abundance figures...\n")
  # Create volcano plots for celiac covariate for both cohorts
  p_usa_celiac <- create_plot_volcano(l_results_usa, "USA", "celiac")
  p_italy_celiac <- create_plot_volcano(l_results_italy, "Italy", "celiac")

  # Save plots to PDF files
  if (!is.null(p_usa_celiac)) {
    ggsave(p_f_figs_usa, p_usa_celiac, width = 10, height = 8, dpi = 300)
    cat("USA celiac volcano plot saved to:", p_f_figs_usa, "\n")
  }

  if (!is.null(p_italy_celiac)) {
    ggsave(p_f_figs_italy, p_italy_celiac, width = 10, height = 8, dpi = 300)
    cat("Italy celiac volcano plot saved to:", p_f_figs_italy, "\n")
  }
}


# Section: Setup function... ----
setup <- function() {
  # Load required libraries
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(openxlsx)
  library(ANCOMBC)
  
  # Set random seed for reproducibility
  set.seed(42)
  
  # Set up paths
  p_dir_wd <<- Sys.getenv(
    "PROJ_CELIAC_PATH",
    "/Users/dhohaabid/Documents/proj_celiac"
  )
  p_dir_out <<- file.path(p_dir_wd, "output")
  p_dir_fig <<- file.path(p_dir_wd, "figures")
  p_dir_res <<- file.path(p_dir_wd, "results")
  
  # Define input files
  p_f_ps_usa <<- file.path(p_dir_out, "ps_gg2_cleaned_usa.RDS")
  p_f_ps_italy <<- file.path(p_dir_out, "ps_gg2_cleaned_italy.RDS")
  
  # Define output files
  p_f_figs_usa <<- file.path(p_dir_fig, "04_differential_abundance_usa.pdf")
  p_f_figs_italy <<- file.path(p_dir_fig, "04_differential_abundance_italy.pdf")
  p_f_excel_usa <<- file.path(p_dir_res, "04_differential_abundance_usa.xlsx")
  p_f_excel_italy <<- file.path(p_dir_res, "04_differential_abundance_italy.xlsx")
  
  # Source utility functions
  # to be able to call the function: add_significance to add sig. indicators
  source(file.path(p_dir_wd, "code", "utils", "utils.R"))
  # to able to extract parameters for plotting like size, font, etc.
  source(file.path(p_dir_wd, "code", "utils", "config.R"))
}


# Section: Define helper functions for this file... ----
# prepare_data: Function to prepare phyloseq object for ANCOM-BC
prepare_data <- function(ps, cohort_name, prev_threshold = 0.03) {
  cat("  Preparing data for", cohort_name, "...\n")
  
  # Step 1: Filter taxa based on prevalence threshold
  cat("    Step 1: Filtering taxa by prevalence...\n")
  v_prevalence <- apply(
    otu_table(ps), ifelse(taxa_are_rows(ps), 1, 2),
    function(x) sum(x > 0) / length(x)
  )
  v_is_taxa_kept <- v_prevalence >= prev_threshold
  ps_filtered <- prune_taxa(v_is_taxa_kept, ps)
  cat("      Total taxa before filtering:", ntaxa(ps), "\n")
  cat("      Taxa after", prev_threshold * 100, "% prevalence filter:", ntaxa(ps_filtered), "\n")

  # Step 2: Set variable types and reference levels
  cat("    Step 2: Setting variable types and reference levels...\n")
  
  # Get metadata
  df_metadata <- data.frame(sample_data(ps_filtered))
  # Set numeric variables as numeric
  df_metadata$time_to_onset <-
    as.numeric(as.character(df_metadata$time_to_onset))
  df_metadata$age_at_gluten_introduction <-
    as.numeric(as.character(df_metadata$age_at_gluten_introduction))
  df_metadata$age_at_solid_introduction <-
    as.numeric(as.character(df_metadata$age_at_solid_introduction))
  # Set categorical variables as factors with proper reference levels
  df_metadata$disease_status <-
    relevel(factor(df_metadata$disease_status), ref = "CONTROL")
  df_metadata$sex <-
    relevel(factor(df_metadata$sex), ref = "Male")
  df_metadata$hla_risk_category <-
    relevel(factor(df_metadata$hla_risk_category), ref = "Standard Risk")
  df_metadata$delivery_mode <-
    relevel(factor(df_metadata$delivery_mode), ref = "Vaginal")
  df_metadata$feeding_type_first_year <-
    relevel(factor(df_metadata$feeding_type_first_year), ref = "Breast_fed")
  # Update phyloseq object with modified sample data
  sample_data(ps_filtered) <- sample_data(df_metadata)

  ps_filtered
}

# fit_model: use ANCOM-BC and extract comprehensive results
fit_model <- function(ps) {
  # Use the same formula structure as alpha diversity analysis
  str_formula <- paste0(
    "disease_status * time_to_onset + sex +",
    "age_at_gluten_introduction + age_at_solid_introduction + ",
    "hla_risk_category + delivery_mode + feeding_type_first_year"
  )
  # Run ANCOM-BC
  tryCatch(
    {
      model <- ancombc2(
        data = ps,
        assay_name = "counts",
        tax_level = NULL, # Use current level (ASV)
        fix_formula = str_formula,
        rand_formula = "(1|patient_id)", # Random effect for patient
        p_adj_method = "fdr",
        pseudo_sens = TRUE,
        s0_perc = 0.05, # Sensitivity parameter
        group = "disease_status",
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,  # Significance level 
        n_cl = 1, # Number of cores
        verbose = TRUE
      )
      
      if (is.null(model)) {
        return(NULL)
      }
      
      # Extract results - ANCOMBC2 has different structure
      res <- model$res
      
      cat("    Creating comprehensive results dataframe with all covariates and taxonomic ranks...\n")

      # Get taxonomic information for all ranks 
      tax_table <- as.data.frame(tax_table(ps))
      tax_table$taxon <- rownames(tax_table)
      
      # Create a list were results of each effect are stored
      l_results <- list()

      # disease_status
      l_results[["celiac"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_disease_statusCELIAC"]],
        std_error = res[["se_disease_statusCELIAC"]],
        p_value = res[["p_disease_statusCELIAC"]],
        p_value_adj = res[["q_disease_statusCELIAC"]],
        W_test = res[["W_disease_statusCELIAC"]],
        diff = res[["diff_disease_statusCELIAC"]],
        diff_robust = res[["diff_robust_disease_statusCELIAC"]],
        passed_ss = res[["passed_ss_disease_statusCELIAC"]],
        sig_adj = sapply(res[["q_disease_statusCELIAC"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_disease_statusCELIAC"]]), ]
      
      # time_to_onset
      l_results[["time_to_onset"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_time_to_onset"]],
        std_error = res[["se_time_to_onset"]],
        p_value = res[["p_time_to_onset"]],
        p_value_adj = res[["q_time_to_onset"]],
        W_test = res[["W_time_to_onset"]],
        diff = res[["diff_time_to_onset"]],
        diff_robust = res[["diff_robust_time_to_onset"]],
        passed_ss = res[["passed_ss_time_to_onset"]],
        sig_adj = sapply(res[["q_time_to_onset"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_time_to_onset"]]), ]
      
      # disease_status:time_to_onset interaction
      l_results[["celiac:time_to_onset"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_disease_statusCELIAC:time_to_onset"]],
        std_error = res[["se_disease_statusCELIAC:time_to_onset"]],
        p_value = res[["p_disease_statusCELIAC:time_to_onset"]],
        p_value_adj = res[["q_disease_statusCELIAC:time_to_onset"]],
        W_test = res[["W_disease_statusCELIAC:time_to_onset"]],
        diff = res[["diff_disease_statusCELIAC:time_to_onset"]],
        diff_robust = res[["diff_robust_disease_statusCELIAC:time_to_onset"]],
        passed_ss = res[["passed_ss_disease_statusCELIAC:time_to_onset"]],
        sig_adj = sapply(res[["q_disease_statusCELIAC:time_to_onset"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_disease_statusCELIAC:time_to_onset"]]), ]
      
      # sex
      l_results[["female"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_sexFemale"]],
        std_error = res[["se_sexFemale"]],
        p_value = res[["p_sexFemale"]],
        p_value_adj = res[["q_sexFemale"]],
        W_test = res[["W_sexFemale"]],
        diff = res[["diff_sexFemale"]],
        diff_robust = res[["diff_robust_sexFemale"]],
        passed_ss = res[["passed_ss_sexFemale"]],
        sig_adj = sapply(res[["q_sexFemale"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_sexFemale"]]), ]
      
      # age_at_gluten_introduction
      l_results[["age_at_gluten_introduction"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_age_at_gluten_introduction"]],
        std_error = res[["se_age_at_gluten_introduction"]],
        p_value = res[["p_age_at_gluten_introduction"]],
        p_value_adj = res[["q_age_at_gluten_introduction"]],
        W_test = res[["W_age_at_gluten_introduction"]],
        diff = res[["diff_age_at_gluten_introduction"]],
        diff_robust = res[["diff_robust_age_at_gluten_introduction"]],
        passed_ss = res[["passed_ss_age_at_gluten_introduction"]],
        sig_adj = sapply(res[["q_age_at_gluten_introduction"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_age_at_gluten_introduction"]]), ]
      
      # age_at_solid_introduction
      l_results[["age_at_solid_introduction"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_age_at_solid_introduction"]],
        std_error = res[["se_age_at_solid_introduction"]],
        p_value = res[["p_age_at_solid_introduction"]],
        p_value_adj = res[["q_age_at_solid_introduction"]],
        W_test = res[["W_age_at_solid_introduction"]],
        diff = res[["diff_age_at_solid_introduction"]],
        diff_robust = res[["diff_robust_age_at_solid_introduction"]],
        passed_ss = res[["passed_ss_age_at_solid_introduction"]],
        sig_adj = sapply(res[["q_age_at_solid_introduction"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_age_at_solid_introduction"]]), ]
      
      # hla_risk_category (high risk)
      l_results[["hla_high_risk"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_hla_risk_categoryHigh Risk"]],
        std_error = res[["se_hla_risk_categoryHigh Risk"]],
        p_value = res[["p_hla_risk_categoryHigh Risk"]],
        p_value_adj = res[["q_hla_risk_categoryHigh Risk"]],
        W_test = res[["W_hla_risk_categoryHigh Risk"]],
        diff = res[["diff_hla_risk_categoryHigh Risk"]],
        diff_robust = res[["diff_robust_hla_risk_categoryHigh Risk"]],
        passed_ss = res[["passed_ss_hla_risk_categoryHigh Risk"]],
        sig_adj = sapply(res[["q_hla_risk_categoryHigh Risk"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_hla_risk_categoryHigh Risk"]]), ]
      
      # hla_risk_category (low risk)
      l_results[["hla_low_risk"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_hla_risk_categoryLow/No Risk"]],
        std_error = res[["se_hla_risk_categoryLow/No Risk"]],
        p_value = res[["p_hla_risk_categoryLow/No Risk"]],
        p_value_adj = res[["q_hla_risk_categoryLow/No Risk"]],
        W_test = res[["W_hla_risk_categoryLow/No Risk"]],
        diff = res[["diff_hla_risk_categoryLow/No Risk"]],
        diff_robust = res[["diff_robust_hla_risk_categoryLow/No Risk"]],
        passed_ss = res[["passed_ss_hla_risk_categoryLow/No Risk"]],
        sig_adj = sapply(res[["q_hla_risk_categoryLow/No Risk"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_hla_risk_categoryLow/No Risk"]]), ]

      # delivery_mode
      l_results[["c_section"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_delivery_modeC-Section"]],
        std_error = res[["se_delivery_modeC-Section"]],
        p_value = res[["p_delivery_modeC-Section"]],
        p_value_adj = res[["q_delivery_modeC-Section"]],
        W_test = res[["W_delivery_modeC-Section"]],
        diff = res[["diff_delivery_modeC-Section"]],
        diff_robust = res[["diff_robust_delivery_modeC-Section"]],
        passed_ss = res[["passed_ss_delivery_modeC-Section"]],
        sig_adj = sapply(res[["q_delivery_modeC-Section"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_delivery_modeC-Section"]]), ]
      
      # feeding_type_first_year (formula)
      l_results[["formula"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_feeding_type_first_yearFormula"]],
        std_error = res[["se_feeding_type_first_yearFormula"]],
        p_value = res[["p_feeding_type_first_yearFormula"]],
        p_value_adj = res[["q_feeding_type_first_yearFormula"]],
        W_test = res[["W_feeding_type_first_yearFormula"]],
        diff = res[["diff_feeding_type_first_yearFormula"]],
        diff_robust = res[["diff_robust_feeding_type_first_yearFormula"]],
        passed_ss = res[["passed_ss_feeding_type_first_yearFormula"]],
        sig_adj = sapply(res[["q_feeding_type_first_yearFormula"]], add_significance),
        stringsAsFactors = FALSE
      )[order(res[["q_feeding_type_first_yearFormula"]]), ]

      # feeding_type_first_year (breastmilk and formula)
      l_results[["breastmilk_and_formula"]] <- data.frame(
        taxon = res$taxon,
        kingdon = tax_table$Kingdom[match(res$taxon, tax_table$taxon)],
        phylum = tax_table$Phylum[match(res$taxon, tax_table$taxon)],
        class = tax_table$Class[match(res$taxon, tax_table$taxon)],
        order = tax_table$Order[match(res$taxon, tax_table$taxon)],
        family = tax_table$Family[match(res$taxon, tax_table$taxon)],
        genus = tax_table$Genus[match(res$taxon, tax_table$taxon)],
        species = tax_table$Species[match(res$taxon, tax_table$taxon)],
        lfc = res[["lfc_feeding_type_first_yearBreastmilk_and_formula"]],
        std_error = res[["se_feeding_type_first_yearBreastmilk_and_formula"]],
        p_value = res[["p_feeding_type_first_yearBreastmilk_and_formula"]],
        p_value_adj = res[["q_feeding_type_first_yearBreastmilk_and_formula"]],
        W_test = res[["W_feeding_type_first_yearBreastmilk_and_formula"]],
        diff = res[["diff_feeding_type_first_yearBreastmilk_and_formula"]],
        diff_robust = res[["diff_robust_feeding_type_first_yearBreastmilk_and_formula"]],
        sig_adj = sapply(res[["q_feeding_type_first_yearBreastmilk_and_formula"]], add_significance),
        passed_ss = res[["passed_ss_feeding_type_first_yearBreastmilk_and_formula"]],
        stringsAsFactors = FALSE
      )[order(res[["q_feeding_type_first_yearBreastmilk_and_formula"]]), ]
      
    },
    error = function(e) {
      cat("    Error in ANCOM-BC analysis:", e$message, "\n")
      return(NULL)
    }
  )
  l_results
}

# Function to save results for a single cohort
save_l_results <- function(l_results, cohort_name, output_file) {
  if (is.null(l_results)) {
    cat("No results available for", cohort_name, "cohort\n")
    return()
  }
  
  cat("Processing", cohort_name, "cohort...\n")
  wb <- createWorkbook()
  
  # Loop through each covariate in the results list
  for (covariate_name in names(l_results)) {
    covariate_df <- l_results[[covariate_name]]
    
    if (!is.null(covariate_df) && nrow(covariate_df) > 0) {
      # Use covariate name as sheet name, replacing colons with underscores
      sheet_name <- gsub(":", "x", covariate_name)
      
      # Add worksheet
      addWorksheet(wb, sheetName = sheet_name)
      
      # Write data to worksheet
      writeData(wb, sheet = sheet_name, x = covariate_df)
      
      # Add some basic formatting for significant results
      if ("p_value_adj" %in% colnames(covariate_df)) {
        sig_rows <- which(covariate_df$p_value_adj < 0.05) + 1 # +1 for header
        if (length(sig_rows) > 0) {
          for (row in sig_rows) {
            addStyle(wb,
              sheet = sheet_name,
              style = createStyle(bgFill = "#FFE6E6"), # Light red for significant
              rows = row, cols = 1:ncol(covariate_df)
            )
          }
        }
      }
      
      cat("  Added sheet:", sheet_name, "with", nrow(covariate_df), "rows\n")
    }
  }
  
  # Save the workbook
  cat("Excel file saved to:", output_file, "\n")
  saveWorkbook(wb, output_file, overwrite = TRUE)
}

# Function to create volcano plot for a single covariate
create_plot_volcano <- function(l_results, cohort_name, covariate_name,
                                output_file = NULL) {
  if (is.null(l_results)) {
    cat("No results available for", cohort_name, "cohort\n")
    return(NULL)
  }
  
  covariate_df <- l_results[[covariate_name]]
  
  if (is.null(covariate_df) || nrow(covariate_df) == 0) {
    cat("No data available for covariate", covariate_name, "\n")
    return(NULL)
  }
  
  # Prepare data for volcano plot
  plot_data <- covariate_df
  plot_data$neg_log10_q <- -log10(plot_data$p_value_adj)
  plot_data$significant <- plot_data$p_value_adj < 0.05
  plot_data$lfc_value <- plot_data$lfc
  
  # Create volcano plot
  p <- ggplot(plot_data, aes(x = lfc_value, y = neg_log10_q)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray", alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    labs(
      title = paste(cohort_name, "-", covariate_name),
      subtitle = paste("Significant taxa (q < 0.05):", sum(plot_data$significant, na.rm = TRUE)),
      x = "Log Fold Change",
      y = "-log10(q-value)",
      color = "Significant"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom"
    )
  return(p)
}

# Section: Execute... ----
cat("
╔══════════════════════════════════════════════════════════════╗
║              DIFFERENTIAL ABUNDANCE ANALYSIS                 ║
║                        STARTING NOW                          ║
╚══════════════════════════════════════════════════════════════╝
")
setup()
main()
