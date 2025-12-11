# PHROG Analysis Script
library(edgeR)
library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

# Load and prepare data
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

# Create PA data
total.PA <- as.data.frame((total.phrog.abundance.clean_0.75_0.03 > 0) * 1)

# Convert to long format for modeling
total_pa_long <- total.PA %>%
  rownames_to_column("phrog_id") %>%
  pivot_longer(cols = -phrog_id, names_to = "sample_id", values_to = "PA") %>%
  left_join(total.phrog.metadata.clean_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id")

print("Data preparation completed")
print(paste("Total PHROGs:", length(unique(total_pa_long$phrog_id))))
print(paste("Total observations:", nrow(total_pa_long)))

# Test with first 10 PHROGs for PA Model 1
phrog_subset <- unique(total_pa_long$phrog_id)[1:10]
print(paste("Testing with PHROGs:", paste(phrog_subset, collapse=", ")))

# PA Model 1: glmmTMB logit for total cohort (subset)
total_glmmTMB_logit_results <- total_pa_long %>%
  filter(phrog_id %in% phrog_subset) %>%
  group_by(phrog_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                family = binomial(link = "logit"), data = .x)
      }, error = function(e) {
        cat("Error for PHROG:", .x$phrog_id[1], "-", e$message, "\n")
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

print("PA Model 1 test completed")
print(paste("Results rows:", nrow(total_glmmTMB_logit_results)))

# Save test results
write.csv(total_glmmTMB_logit_results, "../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_logit_test_results.csv", row.names = FALSE)

print("Test results saved successfully")