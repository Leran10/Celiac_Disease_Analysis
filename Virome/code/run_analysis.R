# Load required libraries
suppressMessages({
  library(edgeR)
  library(glmmTMB)
  library(broom.mixed)
  library(dplyr)
  library(purrr)
  library(tibble)
})

cat("Starting Celiac Phage Analysis...\n")

# Load and process data
total.contig.abundance.clean_0.75_0.03 <- read.csv("../Contig/data/total/total.contig.abundance.clean_0.75_0.03.csv") %>% 
  column_to_rownames("X")
colnames(total.contig.abundance.clean_0.75_0.03) <- gsub("X","",colnames(total.contig.abundance.clean_0.75_0.03))

total.contig.metadata.clean_0.75_0.03 <- read.csv("../Contig/data/total/total.contig.metadata.clean_0.75_0.03.csv") %>% 
  column_to_rownames("X")

# Process metadata factors
total.contig.metadata.clean_0.75_0.03$feeding_first_year <- factor(total.contig.metadata.clean_0.75_0.03$feeding_first_year, 
                                                                   levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
total.contig.metadata.clean_0.75_0.03$HLA.Category <- factor(total.contig.metadata.clean_0.75_0.03$HLA.Category, 
                                                             levels = c("Standard Risk","High Risk","Low/No Risk"))
total.contig.metadata.clean_0.75_0.03$Sex <- factor(total.contig.metadata.clean_0.75_0.03$Sex, levels = c("Female","Male"))
total.contig.metadata.clean_0.75_0.03$Delivery.Mode <- factor(total.contig.metadata.clean_0.75_0.03$Delivery.Mode, 
                                                              levels = c("Vaginal","C-Section"))
total.contig.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(total.contig.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
total.contig.metadata.clean_0.75_0.03$Dx.Status <- factor(total.contig.metadata.clean_0.75_0.03$Dx.Status, 
                                                          levels = c("CONTROL","CELIAC"))
total.contig.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(total.contig.metadata.clean_0.75_0.03$onset_timeline_combined, 
                                                                        levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))

# Create PA table
cat("Creating PA table...\n")
cpm_abundance <- cpm(total.contig.abundance.clean_0.75_0.03)
total.PA <- as.data.frame(ifelse(cpm_abundance >= 1, 1, 0))

cat("PA table created with dimensions:", dim(total.PA), "\n")

# PA Model 1: glmmTMB with logit link
cat("Running PA Model 1: glmmTMB logit link...\n")
total_pa_model1_results <- list()

for(i in 1:nrow(total.PA)) {
  gene_name <- rownames(total.PA)[i]
  gene_data <- data.frame(
    PA = as.numeric(total.PA[i, ]),
    total.contig.metadata.clean_0.75_0.03
  )
  
  tryCatch({
    model <- glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                     HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                     data = gene_data,
                     family = binomial(link = "logit"))
    
    model_summary <- tidy(model, effects = "fixed")
    model_summary$gene <- gene_name
    total_pa_model1_results[[gene_name]] <- model_summary
    
    if(i %% 20 == 0) cat("Processed", i, "genes\n")
  }, error = function(e) {
    cat("Error with gene", gene_name, ":", e$message, "\n")
  })
}

# Combine and save results
if(length(total_pa_model1_results) > 0) {
  total_pa_model1_combined <- do.call(rbind, total_pa_model1_results)
  write.csv(total_pa_model1_combined, "Orf_Contig_Phrog_compositional/results/total_PA_model1_glmmTMB_logit.csv", row.names = FALSE)
  cat("PA Model 1 completed. Results saved.\n")
  cat("Results dimensions:", dim(total_pa_model1_combined), "\n")
} else {
  cat("No successful models for PA Model 1\n")
}

cat("Analysis completed!\n")