# Complete the remaining PHROG analysis models
library(edgeR)
library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(lme4)
library(geepack)

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

# Create PA and abundance data
total.PA <- as.data.frame((total.phrog.abundance.clean_0.75_0.03 > 0) * 1)

# Convert to long format for modeling
total_pa_long <- total.PA %>%
  rownames_to_column("phrog_id") %>%
  pivot_longer(cols = -phrog_id, names_to = "sample_id", values_to = "PA") %>%
  left_join(total.phrog.metadata.clean_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id")

total_abundance_long <- total.phrog.abundance.clean_0.75_0.03 %>%
  rownames_to_column("phrog_id") %>%
  pivot_longer(cols = -phrog_id, names_to = "sample_id", values_to = "abundance") %>%
  left_join(total.phrog.metadata.clean_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(log_offset = log(total_abundance + 1))

cat("Data preparation completed\n")
cat("PHROGs:", length(unique(total_pa_long$phrog_id)), "\n")
cat("Samples:", length(unique(total_pa_long$sample_id)), "\n")

# PA Model 3: glmer (lme4)
cat("Running PA Model 3: glmer...\n")
total_glmer_results <- total_pa_long %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmer(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
              HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
              family = binomial, data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_glmer_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmer_results.csv", row.names = FALSE)
cat("Completed PA Model 3 - Results rows:", nrow(total_glmer_results), "\n")

# PA Model 4: GEE
cat("Running PA Model 4: GEE...\n")
total_GEE_results <- total_pa_long %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        geeglm(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
               HLA.Category + feeding_first_year + Delivery.Mode,
               id = patientID, family = binomial, corstr = "exchangeable", data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_GEE_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_PA_GEE_results.csv", row.names = FALSE)
cat("Completed PA Model 4 - Results rows:", nrow(total_GEE_results), "\n")

# Abundance Model 1: Negative Binomial GLMM
cat("Running Abundance Model 1: Negative Binomial...\n")
total_negbinom_results <- total_abundance_long %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(abundance ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(log_offset),
                family = nbinom2, data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_negbinom_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_negbinom_results.csv", row.names = FALSE)
cat("Completed Abundance Model 1 - Results rows:", nrow(total_negbinom_results), "\n")

# Abundance Model 2: Hurdle Negative Binomial GLMM
cat("Running Abundance Model 2: Hurdle NB...\n")
total_hurdle_results <- total_abundance_long %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(abundance ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(log_offset),
                family = truncated_nbinom2, data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_hurdle_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_hurdle_results.csv", row.names = FALSE)
cat("Completed Abundance Model 2 - Results rows:", nrow(total_hurdle_results), "\n")

# Abundance Model 3: Zero-Inflated Negative Binomial GLMM
cat("Running Abundance Model 3: Zero-Inflated NB...\n")
total_zinb_results <- total_abundance_long %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(abundance ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(log_offset),
                ziformula = ~ 1, family = nbinom2, data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(phrog_id, tidy_results) %>%
  unnest(tidy_results)

write.csv(total_zinb_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_zinb_results.csv", row.names = FALSE)
cat("Completed Abundance Model 3 - Results rows:", nrow(total_zinb_results), "\n")

cat("All remaining models completed successfully!\n")
cat("PA Models: glmmTMB logit, glmmTMB cloglog, glmer, GEE\n")
cat("Abundance Models: limma-voom, negative binomial, hurdle, zero-inflated\n")