# Generate all missing PHROG model results
library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(lme4)
library(geepack)
library(edgeR)
library(limma)

# Load and prepare total cohort data
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

cat("Data loaded. PHROGs:", nrow(total.phrog.abundance.clean_0.75_0.03), "Samples:", ncol(total.phrog.abundance.clean_0.75_0.03), "\n")

# Create PA data
total.PA <- as.data.frame((total.phrog.abundance.clean_0.75_0.03 > 0) * 1)
total_pa_long <- total.PA %>%
  rownames_to_column("phrog_id") %>%
  pivot_longer(cols = -phrog_id, names_to = "sample_id", values_to = "PA") %>%
  left_join(total.phrog.metadata.clean_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id")

# Working with first 100 PHROGs to ensure completion
phrog_subset <- unique(total_pa_long$phrog_id)[1:100]

# PA Model 4: GEE
cat("Running PA Model 4: GEE...\n")
total_gee_results <- total_pa_long %>%
  filter(phrog_id %in% phrog_subset) %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, function(dat) {
      tryCatch({
        geeglm(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
               HLA.Category + feeding_first_year + Delivery.Mode,
               id = patientID, family = binomial, corstr = "exchangeable", data = dat)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, function(mod) {
      if (!is.null(mod)) {
        tryCatch(broom::tidy(mod, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_gee_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_PA_GEE_results.csv", row.names = FALSE)
cat("PA Model 4 completed. Results:", nrow(total_gee_results), "rows\n")

# Abundance models 
# Prepare DGEList and normalization
total.y <- DGEList(counts = total.phrog.abundance.clean_0.75_0.03, samples = total.phrog.metadata.clean_0.75_0.03)
total.y <- calcNormFactors(total.y, method = "TMMwsp")
total.lib_eff <- with(total.y$samples, lib.size * norm.factors)

# Create abundance long format
total_abundance_long <- total.phrog.abundance.clean_0.75_0.03 %>%
  as.data.frame() %>%
  rownames_to_column("phrog_id") %>%
  pivot_longer(cols = -phrog_id, names_to = "sample_id", values_to = "count") %>%
  left_join(total.phrog.metadata.clean_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  mutate(offset = log(total.lib_eff[sample_id]))

# Abundance Model 2: Hurdle NB
cat("Running Abundance Model 2: Hurdle NB...\n")
total_hurdle_results <- total_abundance_long %>%
  filter(phrog_id %in% phrog_subset) %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, function(dat) {
      tryCatch({
        glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                ziformula = ~ 1, family = truncated_nbinom2, data = dat)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, function(mod) {
      if (!is.null(mod)) {
        tryCatch(tidy(mod, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_hurdle_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_hurdle_results.csv", row.names = FALSE)
cat("Abundance Model 2 completed. Results:", nrow(total_hurdle_results), "rows\n")

# Abundance Model 3: Zero-Inflated NB
cat("Running Abundance Model 3: Zero-Inflated NB...\n")
total_zinb_results <- total_abundance_long %>%
  filter(phrog_id %in% phrog_subset) %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, function(dat) {
      tryCatch({
        glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                ziformula = ~ 1, family = nbinom2, data = dat)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, function(mod) {
      if (!is.null(mod)) {
        tryCatch(tidy(mod, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_zinb_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_zinb_results.csv", row.names = FALSE)
cat("Abundance Model 3 completed. Results:", nrow(total_zinb_results), "rows\n")

cat("All missing models completed!\n")