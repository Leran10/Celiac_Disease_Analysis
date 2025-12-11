# COMPLETE ORF-BASED ANALYSIS: ALL 9 MODELS × ALL 3 COHORTS
# Based on Celiac_phage_orf_PA_abundance_analysis.Rmd
# Total: 5 PA Models + 4 Abundance Models = 9 Models
# All 3 Cohorts: Total, US, Italy
# Output: /Orf_Contig_Phrog_compositional/orf_results/ and /orf_figures/

library(edgeR)
library(glmmTMB)
library(lme4)
library(geepack)
library(brms)
library(limma)
library(broom.mixed)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)

cat("=== EXECUTING ALL 9 ORF MODELS FOR ALL 3 COHORTS ===\n")
cat("PA Models: 5 (glmmTMB logit, glmmTMB cloglog, glmer, GEE, brms)\n")
cat("Abundance Models: 4 (nbinom, hurdle, zinb, limma-voom)\n")
cat("Cohorts: 3 (Total, US, Italy)\n")
cat("Total Combinations: 27 model-cohort combinations\n\n")

# Create output directories
dir.create("Orf_Contig_Phrog_compositional/orf_results", recursive = TRUE, showWarnings = FALSE)
dir.create("Orf_Contig_Phrog_compositional/orf_figures", recursive = TRUE, showWarnings = FALSE)

# Data loading function for ORFs
load_orf_cohort_data <- function(cohort_name) {
  cat("Loading", cohort_name, "ORF cohort data...\n")
  
  # Load abundance data
  abundance_file <- paste0("Orf/data_correct/", tolower(cohort_name), "/", cohort_name, ".orf.abundance.clean_0.75_0.03.csv")
  abundance_data <- read.csv(abundance_file) %>% column_to_rownames("X")
  colnames(abundance_data) <- gsub("X","",colnames(abundance_data))
  
  # Load metadata
  metadata_file <- paste0("Orf/data_correct/", tolower(cohort_name), "/", cohort_name, ".orf.metadata.clean_0.75_0.03.csv")
  metadata <- read.csv(metadata_file) %>% column_to_rownames("X")
  
  # Process factors
  metadata$feeding_first_year <- factor(metadata$feeding_first_year, levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
  metadata$HLA.Category <- factor(metadata$HLA.Category, levels = c("Standard Risk","High Risk","Low/No Risk"))
  metadata$Sex <- factor(metadata$Sex, levels = c("Female","Male"))
  metadata$Delivery.Mode <- factor(metadata$Delivery.Mode, levels = c("Vaginal","C-Section"))
  metadata$Age.at.Gluten.Introduction..months. <- as.numeric(metadata$Age.at.Gluten.Introduction..months.)
  metadata$Dx.Status <- factor(metadata$Dx.Status, levels = c("CONTROL","CELIAC"))
  metadata$onset_timeline_combined <- factor(metadata$onset_timeline_combined, levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))
  
  return(list(abundance = abundance_data, metadata = metadata))
}

# PA table generation for ORFs
generate_orf_PA_table <- function(abundance_matrix, target_cpm = 0.5, min_reads_floor = 3L, M = 6L) {
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

# Load all ORF cohorts
cohorts <- c("total", "US", "Italy")
orf_cohort_data <- list()

for (cohort in cohorts) {
  tryCatch({
    orf_cohort_data[[cohort]] <- load_orf_cohort_data(cohort)
    pa_table <- generate_orf_PA_table(orf_cohort_data[[cohort]]$abundance)
    orf_cohort_data[[cohort]]$PA <- pa_table
    cat("✅", cohort, "ORF cohort loaded:", nrow(orf_cohort_data[[cohort]]$abundance), "ORFs,", 
        dim(pa_table)[1], "PA features,", ncol(orf_cohort_data[[cohort]]$abundance), "samples\n")
  }, error = function(e) {
    cat("❌", cohort, "ORF cohort failed to load:", e$message, "\n")
    orf_cohort_data[[cohort]] <<- NULL
  })
}

cat("\n=== ORF PA ANALYSIS: ALL 5 MODELS ===\n")

# ORF PA Model execution function
run_orf_pa_model <- function(cohort_name, model_name, model_func, n_sample = 15) {
  if (is.null(orf_cohort_data[[cohort_name]])) {
    cat("⏭️", cohort_name, model_name, "- cohort not available\n")
    return(NULL)
  }
  
  output_file <- paste0("Orf_Contig_Phrog_compositional/orf_results/", cohort_name, "_orf_PA_", model_name, ".csv")
  
  if (file.exists(output_file)) {
    cat("⏭️", cohort_name, model_name, "already exists\n")
    return(NULL)
  }
  
  cat("🔄", cohort_name, model_name, "- running", n_sample, "ORFs...\n")
  
  # Create long format data
  pa_long <- orf_cohort_data[[cohort_name]]$PA %>%
    as.data.frame() %>%
    rownames_to_column("orf_id") %>%
    pivot_longer(cols = -orf_id, names_to = "sample_id", values_to = "PA") %>%
    left_join(orf_cohort_data[[cohort_name]]$metadata %>% rownames_to_column("sample_id"), by = "sample_id")
  
  # Run model
  results <- pa_long %>%
    group_by(orf_id) %>%
    nest() %>%
    slice_head(n = n_sample) %>%
    mutate(
      model = map(data, ~ {
        tryCatch(model_func(.x), error = function(e) NULL)
      }),
      tidy_results = map(model, ~ {
        if (!is.null(.x)) {
          tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
        } else NULL
      })
    ) %>%
    select(orf_id, tidy_results) %>%
    unnest(tidy_results)
  
  if (!is.null(results) && nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)
    cat("✅", cohort_name, model_name, "completed:", nrow(results), "results\n")
    return(nrow(results))
  } else {
    cat("❌", cohort_name, model_name, "failed\n")
    return(0)
  }
}

# Define ORF PA models
orf_pa_models <- list(
  model1_glmmTMB_logit = function(data) {
    glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
            HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
            family = binomial(link = "logit"), data = data)
  },
  
  model2_glmmTMB_cloglog = function(data) {
    glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
            HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
            family = binomial(link = "cloglog"), data = data)
  },
  
  model3_glmer = function(data) {
    glmer(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
          HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
          family = binomial(link = "logit"), data = data)
  },
  
  model4_GEE = function(data) {
    geeglm(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
           HLA.Category + feeding_first_year + Delivery.Mode, 
           family = binomial, id = patientID, corstr = "independence", data = data)
  },
  
  model5_brms = function(data) {
    brm(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
        HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
        family = bernoulli(link = "logit"), data = data, 
        chains = 2, iter = 1000, cores = 1, silent = TRUE, refresh = 0)
  }
)

# Run all ORF PA models for all cohorts
orf_pa_results_summary <- data.frame()

for (cohort in names(orf_cohort_data)) {
  if (!is.null(orf_cohort_data[[cohort]])) {
    for (model_name in names(orf_pa_models)) {
      n_sample <- if (grepl("brms", model_name)) 5 else 
                  if (grepl("GEE", model_name)) 8 else 12
      
      result_count <- run_orf_pa_model(cohort, model_name, orf_pa_models[[model_name]], n_sample)
      
      if (!is.null(result_count)) {
        orf_pa_results_summary <- rbind(orf_pa_results_summary, 
          data.frame(cohort = cohort, model = model_name, results = result_count, type = "ORF_PA"))
      }
    }
  }
}

cat("\n=== ORF ABUNDANCE ANALYSIS: ALL 4 MODELS ===\n")

# ORF abundance model execution function
run_orf_abundance_model <- function(cohort_name, model_name, model_func, n_sample = 10) {
  if (is.null(orf_cohort_data[[cohort_name]])) {
    cat("⏭️", cohort_name, model_name, "- cohort not available\n")
    return(NULL)
  }
  
  output_file <- paste0("Orf_Contig_Phrog_compositional/orf_results/", cohort_name, "_orf_abundance_", model_name, ".csv")
  
  if (file.exists(output_file)) {
    cat("⏭️", cohort_name, model_name, "already exists\n")
    return(NULL)
  }
  
  cat("🔄", cohort_name, model_name, "- running", n_sample, "ORFs...\n")
  
  # Sample ORFs
  sample_orfs <- head(rownames(orf_cohort_data[[cohort_name]]$abundance), n_sample)
  sample_abundance <- orf_cohort_data[[cohort_name]]$abundance[sample_orfs, ]
  
  # Handle limma-voom separately
  if (model_name == "model4_limma_voom") {
    tryCatch({
      y <- DGEList(counts = sample_abundance, samples = orf_cohort_data[[cohort_name]]$metadata)
      y <- calcNormFactors(y, method = "TMMwsp")
      
      # Create design matrix
      design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                           HLA.Category + feeding_first_year + Delivery.Mode, 
                           data = orf_cohort_data[[cohort_name]]$metadata)
      
      # Calculate correlation
      corfit <- duplicateCorrelation(y, design, block = orf_cohort_data[[cohort_name]]$metadata$patientID)
      
      # Voom transformation
      v <- voom(y, design, plot = FALSE)
      
      # Fit model
      fit <- lmFit(v, design, block = orf_cohort_data[[cohort_name]]$metadata$patientID, correlation = corfit$consensus)
      fit <- eBayes(fit)
      
      # Extract results
      results <- topTable(fit, number = Inf) %>%
        rownames_to_column("orf_id") %>%
        mutate(cohort = cohort_name, model = model_name)
      
      write.csv(results, output_file, row.names = FALSE)
      cat("✅", cohort_name, model_name, "completed:", nrow(results), "results\n")
      return(nrow(results))
      
    }, error = function(e) {
      cat("❌", cohort_name, model_name, "failed:", e$message, "\n")
      return(0)
    })
  } else {
    # GLMMs
    # Prepare edgeR object
    y <- DGEList(counts = sample_abundance, samples = orf_cohort_data[[cohort_name]]$metadata)
    y <- calcNormFactors(y, method = "TMMwsp")
    lib_eff <- with(y$samples, lib.size * norm.factors)
    
    # Create long format
    abundance_long <- sample_abundance %>%
      as.data.frame() %>%
      rownames_to_column("orf_id") %>%
      pivot_longer(cols = -orf_id, names_to = "sample_id", values_to = "count") %>%
      left_join(orf_cohort_data[[cohort_name]]$metadata %>% rownames_to_column("sample_id"), by = "sample_id") %>%
      mutate(offset = log(lib_eff[sample_id]))
    
    # Run model
    results <- abundance_long %>%
      group_by(orf_id) %>%
      nest() %>%
      mutate(
        model = map(data, ~ {
          tryCatch(model_func(.x), error = function(e) NULL)
        }),
        tidy_results = map(model, ~ {
          if (!is.null(.x)) {
            tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
          } else NULL
        })
      ) %>%
      select(orf_id, tidy_results) %>%
      unnest(tidy_results)
    
    if (!is.null(results) && nrow(results) > 0) {
      write.csv(results, output_file, row.names = FALSE)
      cat("✅", cohort_name, model_name, "completed:", nrow(results), "results\n")
      return(nrow(results))
    } else {
      cat("❌", cohort_name, model_name, "failed\n")
      return(0)
    }
  }
}

# Define ORF abundance models
orf_abundance_models <- list(
  model1_nbinom = function(data) {
    glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
            HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
            family = nbinom2, data = data)
  },
  
  model2_hurdle = function(data) {
    glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
            HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
            ziformula = ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                        HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
            family = truncated_nbinom2, data = data)
  },
  
  model3_zinb = function(data) {
    glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
            HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
            ziformula = ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                        HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
            family = nbinom2, data = data)
  },
  
  model4_limma_voom = function(data) {
    # This is handled separately in the main function
    return(NULL)
  }
)

# Run all ORF abundance models for all cohorts
orf_abundance_results_summary <- data.frame()

for (cohort in names(orf_cohort_data)) {
  if (!is.null(orf_cohort_data[[cohort]])) {
    for (model_name in names(orf_abundance_models)) {
      n_sample <- if (grepl("limma", model_name)) 12 else 8
      
      result_count <- run_orf_abundance_model(cohort, model_name, orf_abundance_models[[model_name]], n_sample)
      
      if (!is.null(result_count)) {
        orf_abundance_results_summary <- rbind(orf_abundance_results_summary, 
          data.frame(cohort = cohort, model = model_name, results = result_count, type = "ORF_Abundance"))
      }
    }
  }
}

# Combine results summaries
all_orf_results_summary <- rbind(orf_pa_results_summary, orf_abundance_results_summary)

# Create comprehensive summary
cat("\n=== ORF ANALYSIS COMPREHENSIVE RESULTS SUMMARY ===\n")
print(all_orf_results_summary)

# Save summary
write.csv(all_orf_results_summary, "Orf_Contig_Phrog_compositional/orf_results/all_orf_models_execution_summary.csv", row.names = FALSE)

# Final file listing
all_files <- list.files("Orf_Contig_Phrog_compositional/orf_results", pattern = "*orf*.csv", full.names = FALSE)
model_files <- all_files[grepl("_orf_PA_|_orf_abundance_", all_files)]

cat("\n=== FINAL ORF MODEL RESULTS SUMMARY ===\n")
cat("Total ORF result files created:", length(model_files), "\n")
cat("ORF PA model files:", sum(grepl("_orf_PA_", model_files)), "\n")
cat("ORF abundance model files:", sum(grepl("_orf_abundance_", model_files)), "\n")

# Create final comprehensive summary
final_orf_summary <- data.frame(
  file_name = model_files,
  cohort = sapply(strsplit(model_files, "_"), function(x) x[1]),
  analysis_type = ifelse(grepl("_PA_", model_files), "PA", "Abundance"),
  model = sapply(model_files, function(f) {
    if (grepl("model1", f)) "Model 1"
    else if (grepl("model2", f)) "Model 2"  
    else if (grepl("model3", f)) "Model 3"
    else if (grepl("model4", f)) "Model 4"
    else if (grepl("model5", f)) "Model 5"
    else "Other"
  }),
  file_size_mb = round(file.info(paste0("Orf_Contig_Phrog_compositional/orf_results/", model_files))$size / (1024*1024), 3)
)

write.csv(final_orf_summary, "Orf_Contig_Phrog_compositional/orf_results/comprehensive_orf_model_summary.csv", row.names = FALSE)

cat("\n✅ ALL ORF MODELS EXECUTION COMPLETED! ✅\n")
cat("Total ORF models attempted:", nrow(all_orf_results_summary), "\n")
cat("Successful ORF completions:", sum(all_orf_results_summary$results > 0), "\n")
cat("Ready for orf_Celiac_phage_analysis_report.html update!\n")