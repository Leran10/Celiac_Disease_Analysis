# ORF-based Celiac Disease Phage Analysis - Complete Execution Script
# This script runs the complete ORF analysis pipeline with all models and visualizations

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(edgeR)
library(glmmTMB)
library(lme4)
library(geepack)
library(brms)
library(limma)
library(broom.mixed)
library(ggplot2)
library(pheatmap)
library(viridis)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(vegan)

cat("Starting ORF-based PA and abundance analysis...\n")
cat("Time started:", as.character(Sys.time()), "\n")

# Create output directories
dir.create("../Orf_Contig_Phrog_compositional", showWarnings = FALSE, recursive = TRUE)
dir.create("../Orf_Contig_Phrog_compositional/orf_results", showWarnings = FALSE, recursive = TRUE)
dir.create("../Orf_Contig_Phrog_compositional/orf_figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../Orf_Contig_Phrog_compositional/orf_report_figures", showWarnings = FALSE, recursive = TRUE)

cat("Output directories created\n")

# Load data for all three cohorts
cat("Loading ORF abundance and metadata for all cohorts...\n")

tryCatch({
  # Total cohort
  total.orf.abundance.clean_0.75_0.03 <- read.csv("../Orf/data_correct/total/total.orf.abundance.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  colnames(total.orf.abundance.clean_0.75_0.03) <- gsub("X","",colnames(total.orf.abundance.clean_0.75_0.03))
  
  total.orf.metadata.clean_0.75_0.03 <- read.csv("../Orf/data_correct/total/total.orf.metadata.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  
  # Process metadata factors for total cohort
  total.orf.metadata.clean_0.75_0.03$feeding_first_year <- factor(total.orf.metadata.clean_0.75_0.03$feeding_first_year,levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
  total.orf.metadata.clean_0.75_0.03$HLA.Category <- factor(total.orf.metadata.clean_0.75_0.03$HLA.Category,levels = c("Standard Risk","High Risk","Low/No Risk"))
  total.orf.metadata.clean_0.75_0.03$Sex <- factor(total.orf.metadata.clean_0.75_0.03$Sex,levels = c("Female","Male"))
  total.orf.metadata.clean_0.75_0.03$Delivery.Mode <- factor(total.orf.metadata.clean_0.75_0.03$Delivery.Mode,levels = c("Vaginal","C-Section"))
  total.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(total.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
  total.orf.metadata.clean_0.75_0.03$Dx.Status <- factor(total.orf.metadata.clean_0.75_0.03$Dx.Status,levels = c("CONTROL","CELIAC"))
  total.orf.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(total.orf.metadata.clean_0.75_0.03$onset_timeline_combined,levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))
  
  cat("Total cohort loaded: ORFs =", nrow(total.orf.abundance.clean_0.75_0.03), ", Samples =", ncol(total.orf.abundance.clean_0.75_0.03), "\n")
  
}, error = function(e) {
  cat("Error loading total cohort data:", e$message, "\n")
  stop("Cannot proceed without total cohort data")
})

tryCatch({
  # US cohort
  US.orf.abundance.clean_0.75_0.03 <- read.csv("../Orf/data_correct/US/US.orf.abundance.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  colnames(US.orf.abundance.clean_0.75_0.03) <- gsub("X","",colnames(US.orf.abundance.clean_0.75_0.03))
  
  US.orf.metadata.clean_0.75_0.03 <- read.csv("../Orf/data_correct/US/US.orf.metadata.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  
  # Process metadata factors for US cohort
  US.orf.metadata.clean_0.75_0.03$feeding_first_year <- factor(US.orf.metadata.clean_0.75_0.03$feeding_first_year,levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
  US.orf.metadata.clean_0.75_0.03$HLA.Category <- factor(US.orf.metadata.clean_0.75_0.03$HLA.Category,levels = c("Standard Risk","High Risk","Low/No Risk"))
  US.orf.metadata.clean_0.75_0.03$Sex <- factor(US.orf.metadata.clean_0.75_0.03$Sex,levels = c("Female","Male"))
  US.orf.metadata.clean_0.75_0.03$Delivery.Mode <- factor(US.orf.metadata.clean_0.75_0.03$Delivery.Mode,levels = c("Vaginal","C-Section"))
  US.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(US.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
  US.orf.metadata.clean_0.75_0.03$Dx.Status <- factor(US.orf.metadata.clean_0.75_0.03$Dx.Status,levels = c("CONTROL","CELIAC"))
  US.orf.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(US.orf.metadata.clean_0.75_0.03$onset_timeline_combined,levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))
  
  cat("US cohort loaded: ORFs =", nrow(US.orf.abundance.clean_0.75_0.03), ", Samples =", ncol(US.orf.abundance.clean_0.75_0.03), "\n")
  
}, error = function(e) {
  cat("Error loading US cohort data:", e$message, "\n")
})

tryCatch({
  # Italy cohort
  Italy.orf.abundance.clean_0.75_0.03 <- read.csv("../Orf/data_correct/Italy/Italy.orf.abundance.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  colnames(Italy.orf.abundance.clean_0.75_0.03) <- gsub("X","",colnames(Italy.orf.abundance.clean_0.75_0.03))
  
  Italy.orf.metadata.clean_0.75_0.03 <- read.csv("../Orf/data_correct/Italy/Italy.orf.metadata.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  
  # Process metadata factors for Italy cohort
  Italy.orf.metadata.clean_0.75_0.03$feeding_first_year <- factor(Italy.orf.metadata.clean_0.75_0.03$feeding_first_year,levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
  Italy.orf.metadata.clean_0.75_0.03$HLA.Category <- factor(Italy.orf.metadata.clean_0.75_0.03$HLA.Category,levels = c("Standard Risk","High Risk","Low/No Risk"))
  Italy.orf.metadata.clean_0.75_0.03$Sex <- factor(Italy.orf.metadata.clean_0.75_0.03$Sex,levels = c("Female","Male"))
  Italy.orf.metadata.clean_0.75_0.03$Delivery.Mode <- factor(Italy.orf.metadata.clean_0.75_0.03$Delivery.Mode,levels = c("Vaginal","C-Section"))
  Italy.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(Italy.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
  Italy.orf.metadata.clean_0.75_0.03$Dx.Status <- factor(Italy.orf.metadata.clean_0.75_0.03$Dx.Status,levels = c("CONTROL","CELIAC"))
  Italy.orf.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(Italy.orf.metadata.clean_0.75_0.03$onset_timeline_combined,levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))
  
  cat("Italy cohort loaded: ORFs =", nrow(Italy.orf.abundance.clean_0.75_0.03), ", Samples =", ncol(Italy.orf.abundance.clean_0.75_0.03), "\n")
  
}, error = function(e) {
  cat("Error loading Italy cohort data:", e$message, "\n")
})

# Generate PA tables for all cohorts
cat("Generating Presence/Absence tables...\n")

# Function to generate PA table
generate_PA_table <- function(abundance_matrix, target_cpm = 0.5, min_reads_floor = 3L, M = 6L) {
  # Build DGEList & normalize (TMMwsp)
  y <- DGEList(counts = abundance_matrix)
  y <- y[rowSums(y$counts) > 0, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMMwsp")
  
  # CPM-anchored presence/absence
  cpm_mat <- edgeR::cpm(y, normalized.lib.sizes = TRUE)
  present_logical <- (cpm_mat >= target_cpm) & (y$counts >= min_reads_floor)
  
  # Convert to 0/1 matrix
  PA <- matrix(
    as.integer(present_logical),
    nrow = nrow(present_logical),
    ncol = ncol(present_logical),
    dimnames = dimnames(present_logical)
  )
  
  # Keep genes present in at least M samples
  keep <- rowSums(PA) >= M
  PA_kept <- PA[keep, , drop = FALSE]
  
  return(PA_kept)
}

# Generate PA tables
total.PA <- generate_PA_table(total.orf.abundance.clean_0.75_0.03)
cat("Total PA generated:", dim(total.PA), "\n")

if(exists("US.orf.abundance.clean_0.75_0.03")) {
  US.PA <- generate_PA_table(US.orf.abundance.clean_0.75_0.03)
  cat("US PA generated:", dim(US.PA), "\n")
}

if(exists("Italy.orf.abundance.clean_0.75_0.03")) {
  Italy.PA <- generate_PA_table(Italy.orf.abundance.clean_0.75_0.03)
  cat("Italy PA generated:", dim(Italy.PA), "\n")
}

# Run PA Models - Model 1: glmmTMB logit (main model)
cat("Running PA Model 1: glmmTMB logit...\n")

run_PA_glmmTMB_logit <- function(PA_matrix, metadata, cohort_name) {
  # Create long-format data
  pa_long <- PA_matrix %>%
    as.data.frame() %>%
    rownames_to_column("orf_id") %>%
    pivot_longer(cols = -orf_id, names_to = "sample_id", values_to = "PA") %>%
    left_join(metadata %>% rownames_to_column("sample_id"), by = "sample_id")
  
  # Fit model for each ORF
  results <- pa_long %>%
    group_by(orf_id) %>%
    nest() %>%
    mutate(
      model = map(data, ~ {
        tryCatch({
          glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                  HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                  family = binomial(link = "logit"), data = .x)
        }, error = function(e) NULL)
      }),
      tidy_results = map(model, ~ {
        if (!is.null(.x)) {
          tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
        } else NULL
      })
    ) %>%
    select(orf_id, tidy_results) %>%
    unnest(tidy_results)
  
  # Save results
  write.csv(results, paste0("../Orf_Contig_Phrog_compositional/orf_results/", cohort_name, "_PA_model1_glmmTMB_logit.csv"), row.names = FALSE)
  
  cat(cohort_name, "PA Model 1 completed:", nrow(results), "results\n")
  return(results)
}

# Run PA model 1 for all cohorts
total_pa_results <- run_PA_glmmTMB_logit(total.PA, total.orf.metadata.clean_0.75_0.03, "total")

if(exists("US.PA")) {
  US_pa_results <- run_PA_glmmTMB_logit(US.PA, US.orf.metadata.clean_0.75_0.03, "US")
}

if(exists("Italy.PA")) {
  Italy_pa_results <- run_PA_glmmTMB_logit(Italy.PA, Italy.orf.metadata.clean_0.75_0.03, "Italy")
}

# Run Abundance Models - Model 1: Negative Binomial GLMM (main model)
cat("Running Abundance Model 1: Negative Binomial GLMM...\n")

run_abundance_nbinom <- function(abundance_matrix, metadata, cohort_name) {
  # Prepare DGEList and normalization
  y <- DGEList(counts = abundance_matrix, samples = metadata)
  y <- calcNormFactors(y, method = "TMMwsp")
  lib_eff <- with(y$samples, lib.size * norm.factors)
  
  # Create long-format data
  abundance_long <- abundance_matrix %>%
    as.data.frame() %>%
    rownames_to_column("orf_id") %>%
    pivot_longer(cols = -orf_id, names_to = "sample_id", values_to = "count") %>%
    left_join(metadata %>% rownames_to_column("sample_id"), by = "sample_id") %>%
    mutate(offset = log(lib_eff[sample_id]))
  
  # Fit NB model for each ORF
  results <- abundance_long %>%
    group_by(orf_id) %>%
    nest() %>%
    mutate(
      model = map(data, ~ {
        tryCatch({
          glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                  HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                  family = nbinom2, data = .x)
        }, error = function(e) NULL)
      }),
      tidy_results = map(model, ~ {
        if (!is.null(.x)) {
          tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
        } else NULL
      })
    ) %>%
    select(orf_id, tidy_results) %>%
    unnest(tidy_results)
  
  # Save results
  write.csv(results, paste0("../Orf_Contig_Phrog_compositional/orf_results/", cohort_name, "_abundance_model1_nbinom.csv"), row.names = FALSE)
  
  cat(cohort_name, "Abundance Model 1 completed:", nrow(results), "results\n")
  return(results)
}

# Run abundance model 1 for all cohorts
total_abundance_results <- run_abundance_nbinom(total.orf.abundance.clean_0.75_0.03, total.orf.metadata.clean_0.75_0.03, "total")

if(exists("US.orf.abundance.clean_0.75_0.03")) {
  US_abundance_results <- run_abundance_nbinom(US.orf.abundance.clean_0.75_0.03, US.orf.metadata.clean_0.75_0.03, "US")
}

if(exists("Italy.orf.abundance.clean_0.75_0.03")) {
  Italy_abundance_results <- run_abundance_nbinom(Italy.orf.abundance.clean_0.75_0.03, Italy.orf.metadata.clean_0.75_0.03, "Italy")
}

# Extract timepoint-specific results
cat("Extracting timepoint-specific results...\n")

extract_timepoint_results <- function(model_results, model_name, gene_col = "orf_id") {
  if(is.null(model_results) || nrow(model_results) == 0) {
    return(data.frame())
  }
  
  # Extract interaction terms for timepoints
  timepoint_results <- model_results %>%
    filter(grepl("Dx.StatusCELIAC:onset_timeline_combined", term)) %>%
    mutate(
      timepoint = gsub("Dx.StatusCELIAC:onset_timeline_combined", "", term),
      model = model_name,
      gene = get(gene_col)
    ) %>%
    select(gene, timepoint, estimate, std.error, p.value, conf.low, conf.high, model)
  
  return(timepoint_results)
}

# Extract timepoint results for main models
if(exists("total_pa_results")) {
  total_pa_timepoints <- extract_timepoint_results(total_pa_results, "PA_glmmTMB_logit")
  write.csv(total_pa_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/total_PA_model1_glmmTMB_logit_timepoint_specific_results.csv", row.names = FALSE)
}

if(exists("US_pa_results")) {
  US_pa_timepoints <- extract_timepoint_results(US_pa_results, "PA_glmmTMB_logit")
  write.csv(US_pa_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/US_PA_model1_glmmTMB_logit_timepoint_specific_results.csv", row.names = FALSE)
}

if(exists("Italy_pa_results")) {
  Italy_pa_timepoints <- extract_timepoint_results(Italy_pa_results, "PA_glmmTMB_logit")
  write.csv(Italy_pa_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/Italy_PA_model1_glmmTMB_logit_timepoint_specific_results.csv", row.names = FALSE)
}

if(exists("total_abundance_results")) {
  total_abundance_timepoints <- extract_timepoint_results(total_abundance_results, "abundance_nbinom")
  write.csv(total_abundance_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/total_abundance_model1_nbinom_timepoint_specific_results.csv", row.names = FALSE)
}

if(exists("US_abundance_results")) {
  US_abundance_timepoints <- extract_timepoint_results(US_abundance_results, "abundance_nbinom")
  write.csv(US_abundance_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/US_abundance_model1_nbinom_timepoint_specific_results.csv", row.names = FALSE)
}

if(exists("Italy_abundance_results")) {
  Italy_abundance_timepoints <- extract_timepoint_results(Italy_abundance_results, "abundance_nbinom")
  write.csv(Italy_abundance_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/Italy_abundance_model1_nbinom_timepoint_specific_results.csv", row.names = FALSE)
}

# Generate Analysis Summary
cat("Creating analysis summary...\n")

create_analysis_summary <- function() {
  model_files <- list.files("../Orf_Contig_Phrog_compositional/orf_results/", pattern = "\\.csv$", full.names = TRUE)
  
  summary_data <- data.frame()
  
  for(file in model_files) {
    if(grepl("timepoint_specific", file)) next
    
    tryCatch({
      data <- read.csv(file)
      
      if(nrow(data) > 0 && "p.value" %in% colnames(data)) {
        file_summary <- data.frame(
          file = basename(file),
          cohort = gsub("_.*", "", basename(file)),
          model_type = ifelse(grepl("PA", basename(file)), "PA", "Abundance"),
          total_results = nrow(data),
          significant_results = sum(data$p.value < 0.05, na.rm = TRUE),
          valid_pvalues = sum(!is.na(data$p.value)),
          convergence_rate = round(sum(!is.na(data$p.value)) / nrow(data) * 100, 1)
        )
        
        summary_data <- rbind(summary_data, file_summary)
      }
    }, error = function(e) {
      cat("Error processing file:", file, "\n")
    })
  }
  
  return(summary_data)
}

analysis_summary <- create_analysis_summary()
write.csv(analysis_summary, "../Orf_Contig_Phrog_compositional/orf_results/analysis_summary.csv", row.names = FALSE)

# Print summary
if(exists("analysis_summary") && nrow(analysis_summary) > 0) {
  cat("\n=== ORF ANALYSIS SUMMARY ===\n")
  print(analysis_summary)
}

cat("\nBasic ORF analysis completed!\n")
cat("Time finished:", as.character(Sys.time()), "\n")
cat("Results saved to: ../Orf_Contig_Phrog_compositional/orf_results/\n")
cat("Next: Run visualization scripts to generate figures and reports\n")