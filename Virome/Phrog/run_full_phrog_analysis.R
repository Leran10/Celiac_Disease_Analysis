# Full PHROG Analysis Script
library(edgeR)
library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(lme4)
library(geepack)

# Function to load and prepare data for each cohort
load_cohort_data <- function(cohort_name) {
  # Load abundance data
  abundance_file <- paste0("data/", cohort_name, "/", cohort_name, ".phrog.abundance.clean_0.75_0.03.csv")
  abundance_data <- read.csv(abundance_file) %>%
    column_to_rownames("X")
  colnames(abundance_data) <- gsub("X", "", colnames(abundance_data))
  
  # Load metadata
  metadata_file <- paste0("data/", cohort_name, "/", cohort_name, ".metadata.clean_0.75_0.03.csv")
  metadata <- read.csv(metadata_file) %>%
    column_to_rownames("X")
  
  # Clean metadata
  metadata$feeding_first_year <- factor(metadata$feeding_first_year, levels = c("Breast_fed", "Formula", "Breastmilk_and_formula"))
  metadata$HLA.Category <- factor(metadata$HLA.Category, levels = c("Standard Risk", "High Risk", "Low/No Risk"))
  metadata$Sex <- factor(metadata$Sex, levels = c("Female", "Male"))
  metadata$Delivery.Mode <- factor(metadata$Delivery.Mode, levels = c("Vaginal", "C-Section"))
  metadata$Age.at.Gluten.Introduction..months. <- as.numeric(metadata$Age.at.Gluten.Introduction..months.)
  metadata$Dx.Status <- factor(metadata$Dx.Status, levels = c("CONTROL", "CELIAC"))
  metadata$onset_timeline_combined <- factor(metadata$onset_timeline_combined, levels = c("t0", "t0-6", "t0-12", "t0-18", "t0-24", "t0-30", "t0-36", "t0-over42"))
  
  # Create PA data
  pa_data <- as.data.frame((abundance_data > 0) * 1)
  
  # Convert to long format
  pa_long <- pa_data %>%
    rownames_to_column("phrog_id") %>%
    pivot_longer(cols = -phrog_id, names_to = "sample_id", values_to = "PA") %>%
    left_join(metadata %>% rownames_to_column("sample_id"), by = "sample_id")
  
  abundance_long <- abundance_data %>%
    rownames_to_column("phrog_id") %>%
    pivot_longer(cols = -phrog_id, names_to = "sample_id", values_to = "abundance") %>%
    left_join(metadata %>% rownames_to_column("sample_id"), by = "sample_id") %>%
    group_by(sample_id) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(log_offset = log(total_abundance + 1))
  
  return(list(
    abundance = abundance_data,
    metadata = metadata,
    pa_long = pa_long,
    abundance_long = abundance_long
  ))
}

# Function to run PA models
run_pa_model <- function(data, model_type, cohort_name) {
  cat("Running PA model:", model_type, "for cohort:", cohort_name, "\n")
  
  results <- data$pa_long %>%
    group_by(phrog_id) %>%
    nest() %>%
    mutate(
      model = map(data, ~ {
        tryCatch({
          if (model_type == "glmmTMB_logit") {
            glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                    HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                    family = binomial(link = "logit"), data = .x)
          } else if (model_type == "glmmTMB_cloglog") {
            glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                    HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                    family = binomial(link = "cloglog"), data = .x)
          } else if (model_type == "glmer") {
            glmer(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                  HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                  family = binomial, data = .x)
          } else if (model_type == "GEE") {
            geeglm(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                   HLA.Category + feeding_first_year + Delivery.Mode,
                   id = patientID, family = binomial, corstr = "exchangeable", data = .x)
          }
        }, error = function(e) {
          return(NULL)
        })
      }),
      tidy_results = map(model, ~ {
        if (!is.null(.x)) {
          tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
        } else NULL
      })
    ) %>%
    select(phrog_id, tidy_results) %>%
    unnest(tidy_results)
  
  # Save results
  output_file <- paste0("../Orf_Contig_Phrog_compositional/phrog_results/", cohort_name, "_PA_", model_type, "_results.csv")
  write.csv(results, output_file, row.names = FALSE)
  
  cat("Completed PA model:", model_type, "for cohort:", cohort_name, "- Results rows:", nrow(results), "\n")
  return(results)
}

# Function to run abundance models
run_abundance_model <- function(data, model_type, cohort_name) {
  cat("Running abundance model:", model_type, "for cohort:", cohort_name, "\n")
  
  if (model_type == "limma_voom") {
    # limma-voom analysis
    abundance_matrix <- as.matrix(data$abundance)
    metadata <- data$metadata
    
    # Create design matrix
    design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                          HLA.Category + feeding_first_year + Delivery.Mode, data = metadata)
    
    # Calculate correlation for repeated measures
    corfit <- duplicateCorrelation(abundance_matrix, design, block = metadata$patientID)
    
    # Apply voom transformation
    v <- voom(abundance_matrix, design, block = metadata$patientID, correlation = corfit$consensus)
    
    # Fit linear model
    fit <- lmFit(v, design, block = metadata$patientID, correlation = corfit$consensus)
    fit <- eBayes(fit)
    
    # Extract results for all coefficients
    results <- data.frame()
    for(coef in 2:ncol(design)) {
      coef_results <- topTable(fit, coef = coef, number = Inf, sort.by = "none") %>%
        rownames_to_column("phrog_id") %>%
        mutate(coefficient = colnames(design)[coef])
      results <- rbind(results, coef_results)
    }
  } else {
    # glmmTMB-based models
    results <- data$abundance_long %>%
      group_by(phrog_id) %>%
      nest() %>%
      mutate(
        model = map(data, ~ {
          tryCatch({
            if (model_type == "negbinom") {
              glmmTMB(abundance ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                      HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(log_offset),
                      family = nbinom2, data = .x)
            } else if (model_type == "hurdle") {
              glmmTMB(abundance ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                      HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(log_offset),
                      family = truncated_nbinom2, data = .x)
            } else if (model_type == "zinb") {
              glmmTMB(abundance ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                      HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(log_offset),
                      ziformula = ~ 1, family = nbinom2, data = .x)
            }
          }, error = function(e) {
            return(NULL)
          })
        }),
        tidy_results = map(model, ~ {
          if (!is.null(.x)) {
            tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
          } else NULL
        })
      ) %>%
      select(phrog_id, tidy_results) %>%
      unnest(tidy_results)
  }
  
  # Save results
  output_file <- paste0("../Orf_Contig_Phrog_compositional/phrog_results/", cohort_name, "_abundance_", model_type, "_results.csv")
  write.csv(results, output_file, row.names = FALSE)
  
  cat("Completed abundance model:", model_type, "for cohort:", cohort_name, "- Results rows:", nrow(results), "\n")
  return(results)
}

# Main execution
cat("Starting PHROG Analysis\n")
cat("======================\n")

# Load data for all cohorts
cohorts <- c("total", "US", "Italy")
cohort_data <- list()

for (cohort in cohorts) {
  cat("Loading data for cohort:", cohort, "\n")
  cohort_data[[cohort]] <- load_cohort_data(cohort)
  cat("Cohort:", cohort, "- PHROGs:", nrow(cohort_data[[cohort]]$abundance), "- Samples:", ncol(cohort_data[[cohort]]$abundance), "\n")
}

# Run PA models
pa_models <- c("glmmTMB_logit", "glmmTMB_cloglog", "glmer", "GEE")
for (cohort in cohorts) {
  for (model in pa_models) {
    tryCatch({
      run_pa_model(cohort_data[[cohort]], model, cohort)
    }, error = function(e) {
      cat("ERROR in PA model", model, "for cohort", cohort, ":", e$message, "\n")
    })
  }
}

# Run abundance models
abundance_models <- c("negbinom", "hurdle", "zinb", "limma_voom")
for (cohort in cohorts) {
  for (model in abundance_models) {
    tryCatch({
      run_abundance_model(cohort_data[[cohort]], model, cohort)
    }, error = function(e) {
      cat("ERROR in abundance model", model, "for cohort", cohort, ":", e$message, "\n")
    })
  }
}

cat("PHROG Analysis Completed!\n")