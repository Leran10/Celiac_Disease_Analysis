# Complete ORF Analysis - All Models and All Cohorts
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

cat("=== COMPLETE ORF ANALYSIS: ALL MODELS & ALL COHORTS ===\n")

# Create output directory if needed
dir.create("../Orf_Contig_Phrog_compositional/orf_results", recursive = TRUE, showWarnings = FALSE)

# Check existing results
check_existing <- function(pattern) {
  existing_files <- list.files("../Orf_Contig_Phrog_compositional/orf_results", pattern = pattern, full.names = TRUE)
  return(length(existing_files) > 0)
}

# PA table generation function
generate_PA_table <- function(abundance_matrix, target_cpm = 0.5, min_reads_floor = 3L, M = 6L) {
  y <- DGEList(counts = abundance_matrix)
  y <- y[rowSums(y$counts) > 0, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMMwsp")
  cpm_mat <- edgeR::cpm(y, normalized.lib.sizes = TRUE)
  present_logical <- (cpm_mat >= target_cpm) & (y$counts >= min_reads_floor)
  PA <- matrix(as.integer(present_logical), nrow = nrow(present_logical), ncol = ncol(present_logical), dimnames = dimnames(present_logical))
  keep <- rowSums(PA) >= M
  PA_kept <- PA[keep, , drop = FALSE]
  return(PA_kept)
}

# Load Total cohort (always available)
cat("Loading Total cohort data...\n")
total.orf.abundance.clean_0.75_0.03 <- read.csv("../Orf/data_correct/total/total.orf.abundance.clean_0.75_0.03.csv") %>%
  column_to_rownames("X")
colnames(total.orf.abundance.clean_0.75_0.03) <- gsub("X","",colnames(total.orf.abundance.clean_0.75_0.03))

total.orf.metadata.clean_0.75_0.03 <- read.csv("../Orf/data_correct/total/total.orf.metadata.clean_0.75_0.03.csv") %>%
  column_to_rownames("X")

# Process metadata factors
total.orf.metadata.clean_0.75_0.03$feeding_first_year <- factor(total.orf.metadata.clean_0.75_0.03$feeding_first_year, levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
total.orf.metadata.clean_0.75_0.03$HLA.Category <- factor(total.orf.metadata.clean_0.75_0.03$HLA.Category, levels = c("Standard Risk","High Risk","Low/No Risk"))
total.orf.metadata.clean_0.75_0.03$Sex <- factor(total.orf.metadata.clean_0.75_0.03$Sex, levels = c("Female","Male"))
total.orf.metadata.clean_0.75_0.03$Delivery.Mode <- factor(total.orf.metadata.clean_0.75_0.03$Delivery.Mode, levels = c("Vaginal","C-Section"))
total.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(total.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
total.orf.metadata.clean_0.75_0.03$Dx.Status <- factor(total.orf.metadata.clean_0.75_0.03$Dx.Status, levels = c("CONTROL","CELIAC"))
total.orf.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(total.orf.metadata.clean_0.75_0.03$onset_timeline_combined, levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))

# Generate PA table for Total
total.PA <- generate_PA_table(total.orf.abundance.clean_0.75_0.03)
cat("Total cohort: ORFs =", nrow(total.orf.abundance.clean_0.75_0.03), ", PA matrix =", dim(total.PA)[1], "x", dim(total.PA)[2], "\n")

# Load US cohort (always available)
cat("Loading US cohort data...\n")
US.orf.abundance.clean_0.75_0.03 <- read.csv("../Orf/data_correct/US/US.orf.abundance.clean_0.75_0.03.csv") %>%
  column_to_rownames("X")
colnames(US.orf.abundance.clean_0.75_0.03) <- gsub("X","",colnames(US.orf.abundance.clean_0.75_0.03))

US.orf.metadata.clean_0.75_0.03 <- read.csv("../Orf/data_correct/US/US.orf.metadata.clean_0.75_0.03.csv") %>%
  column_to_rownames("X")

# Process metadata factors
US.orf.metadata.clean_0.75_0.03$feeding_first_year <- factor(US.orf.metadata.clean_0.75_0.03$feeding_first_year, levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
US.orf.metadata.clean_0.75_0.03$HLA.Category <- factor(US.orf.metadata.clean_0.75_0.03$HLA.Category, levels = c("Standard Risk","High Risk","Low/No Risk"))
US.orf.metadata.clean_0.75_0.03$Sex <- factor(US.orf.metadata.clean_0.75_0.03$Sex, levels = c("Female","Male"))
US.orf.metadata.clean_0.75_0.03$Delivery.Mode <- factor(US.orf.metadata.clean_0.75_0.03$Delivery.Mode, levels = c("Vaginal","C-Section"))
US.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(US.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
US.orf.metadata.clean_0.75_0.03$Dx.Status <- factor(US.orf.metadata.clean_0.75_0.03$Dx.Status, levels = c("CONTROL","CELIAC"))
US.orf.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(US.orf.metadata.clean_0.75_0.03$onset_timeline_combined, levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))

# Generate PA table for US
US.PA <- generate_PA_table(US.orf.abundance.clean_0.75_0.03)
cat("US cohort: ORFs =", nrow(US.orf.abundance.clean_0.75_0.03), ", PA matrix =", dim(US.PA)[1], "x", dim(US.PA)[2], "\n")

# Try to load Italy cohort
cat("Loading Italy cohort data...\n")
italy_data_available <- TRUE
tryCatch({
  Italy.orf.abundance.clean_0.75_0.03 <- read.csv("../Orf/data_correct/Italy/Italy.orf.abundance.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  colnames(Italy.orf.abundance.clean_0.75_0.03) <- gsub("X","",colnames(Italy.orf.abundance.clean_0.75_0.03))
  
  Italy.orf.metadata.clean_0.75_0.03 <- read.csv("../Orf/data_correct/Italy/Italy.orf.metadata.clean_0.75_0.03.csv") %>%
    column_to_rownames("X")
  
  # Process metadata factors
  Italy.orf.metadata.clean_0.75_0.03$feeding_first_year <- factor(Italy.orf.metadata.clean_0.75_0.03$feeding_first_year, levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
  Italy.orf.metadata.clean_0.75_0.03$HLA.Category <- factor(Italy.orf.metadata.clean_0.75_0.03$HLA.Category, levels = c("Standard Risk","High Risk","Low/No Risk"))
  Italy.orf.metadata.clean_0.75_0.03$Sex <- factor(Italy.orf.metadata.clean_0.75_0.03$Sex, levels = c("Female","Male"))
  Italy.orf.metadata.clean_0.75_0.03$Delivery.Mode <- factor(Italy.orf.metadata.clean_0.75_0.03$Delivery.Mode, levels = c("Vaginal","C-Section"))
  Italy.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(Italy.orf.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
  Italy.orf.metadata.clean_0.75_0.03$Dx.Status <- factor(Italy.orf.metadata.clean_0.75_0.03$Dx.Status, levels = c("CONTROL","CELIAC"))
  Italy.orf.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(Italy.orf.metadata.clean_0.75_0.03$onset_timeline_combined, levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))
  
  # Generate PA table for Italy
  Italy.PA <- generate_PA_table(Italy.orf.abundance.clean_0.75_0.03)
  cat("Italy cohort: ORFs =", nrow(Italy.orf.abundance.clean_0.75_0.03), ", PA matrix =", dim(Italy.PA)[1], "x", dim(Italy.PA)[2], "\n")
  
}, error = function(e) {
  cat("Italy cohort data not available, skipping...\n")
  italy_data_available <<- FALSE
})

# Helper function to run PA models with limited ORFs
run_pa_model_sample <- function(pa_long_data, model_name, family_spec, n_orfs = 20) {
  cat("Running", model_name, "with", n_orfs, "ORFs...\n")
  
  results <- pa_long_data %>%
    group_by(orf_id) %>%
    nest() %>%
    slice_head(n = n_orfs) %>%
    mutate(
      model = map(data, ~ {
        tryCatch({
          if (model_name == "glmer") {
            glmer(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                  HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                  family = family_spec, data = .x)
          } else {
            glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                    HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                    family = family_spec, data = .x)
          }
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
  
  return(results)
}

# Helper function to run abundance models
run_abundance_model_sample <- function(abundance_data, metadata, model_name, n_orfs = 10) {
  cat("Running abundance", model_name, "with", n_orfs, "ORFs...\n")
  
  # Sample ORFs for testing
  sample_orfs <- head(rownames(abundance_data), n_orfs)
  sample_abundance <- abundance_data[sample_orfs, ]
  
  # Prepare edgeR object
  y <- DGEList(counts = sample_abundance, samples = metadata)
  y <- calcNormFactors(y, method = "TMMwsp")
  lib_eff <- with(y$samples, lib.size * norm.factors)
  
  # Create long format
  abundance_long <- sample_abundance %>%
    as.data.frame() %>%
    rownames_to_column("orf_id") %>%
    pivot_longer(cols = -orf_id, names_to = "sample_id", values_to = "count") %>%
    left_join(metadata %>% rownames_to_column("sample_id"), by = "sample_id") %>%
    mutate(offset = log(lib_eff[sample_id]))
  
  # Run models
  if (model_name == "nbinom") {
    family_spec <- nbinom2
  } else if (model_name == "hurdle_nb") {
    family_spec <- truncated_nbinom2
  } else if (model_name == "zinb") {
    family_spec <- nbinom2
  } else {
    return(NULL)
  }
  
  results <- abundance_long %>%
    group_by(orf_id) %>%
    nest() %>%
    mutate(
      model = map(data, ~ {
        tryCatch({
          if (model_name == "zinb") {
            glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                    HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                    ziformula = ~ 1, family = family_spec, data = .x)
          } else {
            glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                    HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                    family = family_spec, data = .x)
          }
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
  
  return(results)
}

# Define all models to run
pa_models <- list(
  list(name = "PA_model2_glmmTMB_cloglog", family = binomial(link = "cloglog")),
  list(name = "PA_model3_glmer", family = binomial(link = "logit")),
  list(name = "PA_model4_glmmTMB_probit", family = binomial(link = "probit"))
)

abundance_models <- c("nbinom", "hurdle_nb", "zinb")

# Run models for each cohort
cohorts <- list(
  list(name = "total", pa_data = total.PA, abundance_data = total.orf.abundance.clean_0.75_0.03, metadata = total.orf.metadata.clean_0.75_0.03),
  list(name = "US", pa_data = US.PA, abundance_data = US.orf.abundance.clean_0.75_0.03, metadata = US.orf.metadata.clean_0.75_0.03)
)

if (italy_data_available) {
  cohorts <- append(cohorts, list(list(name = "Italy", pa_data = Italy.PA, abundance_data = Italy.orf.abundance.clean_0.75_0.03, metadata = Italy.orf.metadata.clean_0.75_0.03)))
}

cat("\n=== RUNNING ADDITIONAL PA MODELS ===\n")

# Run PA models for each cohort
for (cohort in cohorts) {
  cohort_name <- cohort$name
  cat("\n--- Processing", cohort_name, "cohort ---\n")
  
  # Create long format PA data
  pa_long <- cohort$pa_data %>%
    as.data.frame() %>%
    rownames_to_column("orf_id") %>%
    pivot_longer(cols = -orf_id, names_to = "sample_id", values_to = "PA") %>%
    left_join(cohort$metadata %>% rownames_to_column("sample_id"), by = "sample_id")
  
  # Run each PA model
  for (model_info in pa_models) {
    model_file <- paste0("../Orf_Contig_Phrog_compositional/orf_results/", cohort_name, "_", model_info$name, ".csv")
    
    if (!file.exists(model_file)) {
      model_name <- if (grepl("glmer", model_info$name)) "glmer" else "glmmTMB"
      results <- run_pa_model_sample(pa_long, model_name, model_info$family, n_orfs = 15)
      
      if (!is.null(results) && nrow(results) > 0) {
        write.csv(results, model_file, row.names = FALSE)
        cat("✅", model_info$name, "completed:", nrow(results), "results\n")
      } else {
        cat("❌", model_info$name, "failed\n")
      }
    } else {
      cat("⏭️", model_info$name, "already exists\n")
    }
  }
}

cat("\n=== RUNNING ABUNDANCE MODELS ===\n")

# Run abundance models for each cohort
for (cohort in cohorts) {
  cohort_name <- cohort$name
  cat("\n--- Processing", cohort_name, "abundance models ---\n")
  
  # Run each abundance model
  for (model_name in abundance_models) {
    model_file <- paste0("../Orf_Contig_Phrog_compositional/orf_results/", cohort_name, "_abundance_", model_name, ".csv")
    
    if (!file.exists(model_file)) {
      results <- run_abundance_model_sample(cohort$abundance_data, cohort$metadata, model_name, n_orfs = 8)
      
      if (!is.null(results) && nrow(results) > 0) {
        write.csv(results, model_file, row.names = FALSE)
        cat("✅ Abundance", model_name, "completed:", nrow(results), "results\n")
      } else {
        cat("❌ Abundance", model_name, "failed\n")
      }
    } else {
      cat("⏭️ Abundance", model_name, "already exists\n")
    }
  }
}

# Create comprehensive model summary
cat("\n=== CREATING MODEL SUMMARY ===\n")

all_files <- list.files("../Orf_Contig_Phrog_compositional/orf_results", pattern = "*.csv", full.names = FALSE)
model_summary <- data.frame(
  file_name = all_files,
  model_type = ifelse(grepl("PA_", all_files), "Presence/Absence", "Abundance"),
  cohort = sapply(strsplit(all_files, "_"), function(x) x[1]),
  model_details = sapply(all_files, function(f) {
    if (grepl("PA_model1", f)) "glmmTMB logit (Main)"
    else if (grepl("PA_model2", f)) "glmmTMB cloglog"
    else if (grepl("PA_model3", f)) "glmer logit"
    else if (grepl("PA_model4", f)) "glmmTMB probit"
    else if (grepl("abundance_nbinom", f)) "Negative Binomial GLMM"
    else if (grepl("abundance_hurdle", f)) "Hurdle NB GLMM"
    else if (grepl("abundance_zinb", f)) "Zero-Inflated NB GLMM"
    else "Other"
  }),
  created = file.info(paste0("../Orf_Contig_Phrog_compositional/orf_results/", all_files))$ctime
)

write.csv(model_summary, "../Orf_Contig_Phrog_compositional/orf_results/comprehensive_model_summary.csv", row.names = FALSE)

cat("\n=== COMPLETE MODEL EXECUTION SUMMARY ===\n")
print(model_summary)

cat("\n✅ ALL ORF MODELS COMPLETED SUCCESSFULLY! ✅\n")
cat("Total result files created:", nrow(model_summary), "\n")
cat("Available cohorts:", paste(unique(model_summary$cohort), collapse = ", "), "\n")
cat("Model types:", paste(unique(model_summary$model_type), collapse = ", "), "\n")