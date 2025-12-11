# Run remaining PA models only
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
  
  return(list(pa_long = pa_long))
}

# Load total cohort data
cat("Loading total cohort data...\n")
total_data <- load_cohort_data("total")

# Run remaining PA models for total cohort
cat("Running glmer for total cohort...\n")
total_glmer_results <- total_data$pa_long %>%
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
cat("Completed glmer for total cohort - Results rows:", nrow(total_glmer_results), "\n")

cat("Running GEE for total cohort...\n")
total_GEE_results <- total_data$pa_long %>%
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
cat("Completed GEE for total cohort - Results rows:", nrow(total_GEE_results), "\n")

cat("PA models for total cohort completed!\n")