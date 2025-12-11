# PHROG Abundance Analysis
library(edgeR)
library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(limma)

# Load and prepare total cohort data for abundance analysis
total.phrog.abundance.clean_0.75_0.03 <- read.csv("data/total/total.phrog.abundance.clean_0.75_0.03.csv") %>%
  column_to_rownames("X")
colnames(total.phrog.abundance.clean_0.75_0.03) <- gsub("X","",colnames(total.phrog.abundance.clean_0.75_0.03))

total.phrog.metadata.clean_0.75_0.03 <- read.csv("data/total/total.metadata.clean_0.75_0.03.csv") %>%
  column_to_rownames("X")

# Clean metadata
total.phrog.metadata.clean_0.75_0.03$feeding_first_year <- factor(total.phrog.metadata.clean_0.75_0.03$feeding_first_year,levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
total.phrog.metadata.clean_0.75_0.03$HLA.Category <- factor(total.phrog.metadata.clean_0.75_0.03$HLA.Category,levels = c("Standard Risk","High Risk","Low/No Risk"))
total.phrog.metadata.clean_0.75_0.03$Sex <- factor(total.phrog.metadata.clean_0.75_0.03$Sex,levels = c("Female","Male"))
total.phrog.metadata.clean_0.75_0.03$Delivery.Mode <- factor(total.phrog.metadata.clean_0.75_0.03$Delivery.Mode,levels = c("Vaginal","C-Section"))
total.phrog.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(total.phrog.metadata.clean_0.75_0.03$Age.at.Gluten.Introduction..months.)
total.phrog.metadata.clean_0.75_0.03$Dx.Status <- factor(total.phrog.metadata.clean_0.75_0.03$Dx.Status,levels = c("CONTROL","CELIAC"))
total.phrog.metadata.clean_0.75_0.03$onset_timeline_combined <- factor(total.phrog.metadata.clean_0.75_0.03$onset_timeline_combined,levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))

# Prepare abundance data for limma-voom analysis
abundance_matrix <- as.matrix(total.phrog.abundance.clean_0.75_0.03)
metadata <- total.phrog.metadata.clean_0.75_0.03

cat("Running limma-voom abundance analysis for total cohort...\n")
cat("PHROGs:", nrow(abundance_matrix), "Samples:", ncol(abundance_matrix), "\n")

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
total_limma_results <- data.frame()
for(coef in 2:ncol(design)) {
  coef_results <- topTable(fit, coef = coef, number = Inf, sort.by = "none") %>%
    rownames_to_column("phrog_id") %>%
    mutate(coefficient = colnames(design)[coef])
  total_limma_results <- rbind(total_limma_results, coef_results)
}

# Save results
write.csv(total_limma_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_limma_results.csv", row.names = FALSE)

cat("Completed limma-voom abundance analysis\n")
cat("Results rows:", nrow(total_limma_results), "\n")

# Summary of Dx.Status abundance effects
dx_abundance_effects <- total_limma_results %>% 
  filter(coefficient == "Dx.StatusCELIAC") %>%
  filter(!is.na(P.Value))

dx_abundance_sig <- dx_abundance_effects %>% filter(P.Value < 0.05)

cat("Dx.Status abundance effects:\n")
cat("- Total PHROGs with results:", nrow(dx_abundance_effects), "\n")
cat("- Significant PHROGs (p < 0.05):", nrow(dx_abundance_sig), "\n")

# Save top significant abundance results
if (nrow(dx_abundance_sig) > 0) {
  top_abundance_phrogs <- dx_abundance_sig %>%
    arrange(P.Value) %>%
    select(phrog_id, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
    head(20)
  
  write.csv(top_abundance_phrogs, "../Orf_Contig_Phrog_compositional/phrog_results/top_significant_phrogs_abundance_dx_status.csv", row.names = FALSE)
  cat("Top 20 significant abundance PHROGs saved\n")
}

cat("Abundance analysis completed!\n")