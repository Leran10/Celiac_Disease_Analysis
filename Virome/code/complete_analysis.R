# Complete Celiac Phage Analysis Script
# Load required libraries
suppressMessages({
  library(edgeR)
  library(glmmTMB)
  library(broom.mixed)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(lme4)
  library(geepack)
  library(brms)
  library(limma)
})

cat("Starting Complete Celiac Phage Analysis...\n")

# Function to load and process data for each cohort
load_and_process_cohort <- function(cohort_name) {
  cat("Loading", cohort_name, "cohort data...\n")
  
  # Load abundance data
  abundance_file <- paste0("../Contig/data/", cohort_name, "/", cohort_name, ".contig.abundance.clean_0.75_0.03.csv")
  abundance_data <- read.csv(abundance_file) %>% column_to_rownames("X")
  colnames(abundance_data) <- gsub("X","",colnames(abundance_data))
  
  # Load metadata
  metadata_file <- paste0("../Contig/data/", cohort_name, "/", cohort_name, ".contig.metadata.clean_0.75_0.03.csv")
  metadata_data <- read.csv(metadata_file) %>% column_to_rownames("X")
  
  # Process metadata factors
  metadata_data$feeding_first_year <- factor(metadata_data$feeding_first_year, 
                                           levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
  metadata_data$HLA.Category <- factor(metadata_data$HLA.Category, 
                                      levels = c("Standard Risk","High Risk","Low/No Risk"))
  metadata_data$Sex <- factor(metadata_data$Sex, levels = c("Female","Male"))
  metadata_data$Delivery.Mode <- factor(metadata_data$Delivery.Mode, 
                                       levels = c("Vaginal","C-Section"))
  metadata_data$Age.at.Gluten.Introduction..months. <- as.numeric(metadata_data$Age.at.Gluten.Introduction..months.)
  metadata_data$Dx.Status <- factor(metadata_data$Dx.Status, levels = c("CONTROL","CELIAC"))
  metadata_data$onset_timeline_combined <- factor(metadata_data$onset_timeline_combined, 
                                                 levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))
  
  # Create PA table
  cpm_abundance <- cpm(abundance_data)
  pa_data <- as.data.frame(ifelse(cpm_abundance >= 1, 1, 0))
  
  return(list(abundance = abundance_data, metadata = metadata_data, pa = pa_data))
}

# Load all cohorts
cohorts <- list()
for(cohort in c("total", "US", "Italy")) {
  cohorts[[cohort]] <- load_and_process_cohort(cohort)
}

cat("All cohort data loaded successfully!\n")

# Function to run PA models
run_pa_models <- function(cohort_name, cohort_data) {
  cat("Running PA models for", cohort_name, "cohort...\n")
  
  pa_data <- cohort_data$pa
  metadata <- cohort_data$metadata
  results <- list()
  
  # PA Model 1: glmmTMB logit
  cat("  PA Model 1: glmmTMB logit...\n")
  model1_results <- list()
  for(i in 1:nrow(pa_data)) {
    gene_name <- rownames(pa_data)[i]
    gene_data <- data.frame(PA = as.numeric(pa_data[i, ]), metadata)
    
    tryCatch({
      model <- glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                       HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                       data = gene_data, family = binomial(link = "logit"))
      model_summary <- tidy(model, effects = "fixed")
      model_summary$gene <- gene_name
      model1_results[[gene_name]] <- model_summary
    }, error = function(e) {})
  }
  if(length(model1_results) > 0) {
    results$model1 <- do.call(rbind, model1_results)
    write.csv(results$model1, paste0("Orf_Contig_Phrog_compositional/results/", cohort_name, "_PA_model1_glmmTMB_logit.csv"), row.names = FALSE)
  }
  
  # PA Model 2: glmmTMB cloglog
  cat("  PA Model 2: glmmTMB cloglog...\n")
  model2_results <- list()
  for(i in 1:nrow(pa_data)) {
    gene_name <- rownames(pa_data)[i]
    gene_data <- data.frame(PA = as.numeric(pa_data[i, ]), metadata)
    
    tryCatch({
      model <- glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                       HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                       data = gene_data, family = binomial(link = "cloglog"))
      model_summary <- tidy(model, effects = "fixed")
      model_summary$gene <- gene_name
      model2_results[[gene_name]] <- model_summary
    }, error = function(e) {})
  }
  if(length(model2_results) > 0) {
    results$model2 <- do.call(rbind, model2_results)
    write.csv(results$model2, paste0("Orf_Contig_Phrog_compositional/results/", cohort_name, "_PA_model2_glmmTMB_cloglog.csv"), row.names = FALSE)
  }
  
  # PA Model 3: glmer
  cat("  PA Model 3: glmer...\n")
  model3_results <- list()
  for(i in 1:nrow(pa_data)) {
    gene_name <- rownames(pa_data)[i]
    gene_data <- data.frame(PA = as.numeric(pa_data[i, ]), metadata)
    
    tryCatch({
      model <- glmer(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                     HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                     data = gene_data, family = binomial)
      model_summary <- tidy(model, effects = "fixed")
      model_summary$gene <- gene_name
      model3_results[[gene_name]] <- model_summary
    }, error = function(e) {})
  }
  if(length(model3_results) > 0) {
    results$model3 <- do.call(rbind, model3_results)
    write.csv(results$model3, paste0("Orf_Contig_Phrog_compositional/results/", cohort_name, "_PA_model3_glmer.csv"), row.names = FALSE)
  }
  
  cat("PA models completed for", cohort_name, "cohort\n")
  return(results)
}

# Run PA models for all cohorts
pa_results <- list()
for(cohort_name in names(cohorts)) {
  pa_results[[cohort_name]] <- run_pa_models(cohort_name, cohorts[[cohort_name]])
}

cat("All PA models completed!\n")