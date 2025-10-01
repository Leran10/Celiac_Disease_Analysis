#!/usr/bin/env Rscript
# Compositional Analysis for Corrected Viral ORF Data
# Based on the original CLAUDE.md template but using corrected data
# Author: Claude AI
# Date: 2025-07-25

library(dplyr)
library(vegan)
library(ggplot2)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(ggrepel)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

cat("=== COMPOSITIONAL ANALYSIS FOR CORRECTED VIRAL ORF DATA ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD CORRECTED DATA
# ============================================================================

cat("1. Loading corrected data...\n")

# Load corrected ORF abundance data
orf_data_raw <- read.csv("total_orf.abundance.clean.csv")
rownames(orf_data_raw) <- orf_data_raw$X
orf_data <- orf_data_raw[, -1]  # Remove X column

# Clean sample names by removing X prefix
orf_samples <- colnames(orf_data)
orf_samples_clean <- gsub("^X", "", orf_samples)
colnames(orf_data) <- orf_samples_clean

# Load corrected metadata
metadata <- read.csv("total_metadata.clean.csv")

# Match samples between datasets using X column (contains sample IDs)
sample_ids_from_cols <- colnames(orf_data)
sample_ids_from_metadata <- metadata$X

# Find common samples
common_samples <- intersect(sample_ids_from_cols, sample_ids_from_metadata)
metadata_matched <- metadata[metadata$X %in% common_samples, ]
orf_data_matched <- orf_data[, common_samples]

# Set row names for metadata to match sample names
rownames(metadata_matched) <- metadata_matched$X

cat("Data loaded:\n")
cat("  ORF count:", nrow(orf_data_matched), "ORFs\n")
cat("  Sample count:", ncol(orf_data_matched), "samples\n")
cat("  Matched samples:", length(common_samples), "\n")

# Prepare metadata with required columns
metadata_clean <- metadata_matched %>%
  mutate(
    sample_id = X,
    Dx.Status = factor(Dx.Status, levels = c("CELIAC", "CONTROL")),
    patientID = factor(patientID)
  )

# Add other factors if they exist
if("Country" %in% colnames(metadata_clean)) {
  metadata_clean$Country <- factor(metadata_clean$Country)
}
if("Sex" %in% colnames(metadata_clean)) {
  metadata_clean$Sex <- factor(metadata_clean$Sex)
}
if("HLA.Category" %in% colnames(metadata_clean)) {
  metadata_clean$HLA.Category <- factor(metadata_clean$HLA.Category)
}
if("feeding_first_year" %in% colnames(metadata_clean)) {
  metadata_clean$feeding_first_year <- factor(metadata_clean$feeding_first_year)
}
if("Delivery.Mode" %in% colnames(metadata_clean)) {
  metadata_clean$Delivery.Mode <- factor(metadata_clean$Delivery.Mode)
}

# Remove rows with missing critical data
metadata_clean <- metadata_clean[complete.cases(metadata_clean[, c("sample_id", "patientID", "Dx.Status", "onset_timeline_numeric")]), ]

cat("Clean metadata samples:", nrow(metadata_clean), "\n")
cat("Patients:", length(unique(metadata_clean$patientID)), "\n")
cat("Disease groups:\n")
print(table(metadata_clean$Dx.Status))

# ============================================================================
# PART 2: CALCULATE DIVERSITY METRICS
# ============================================================================

cat("\n2. Calculating diversity metrics...\n")

# Function to calculate comprehensive diversity metrics
calculate_diversity_metrics <- function(abundance_matrix, metadata) {
  
  diversity_data <- data.frame(
    sample_id = colnames(abundance_matrix),
    stringsAsFactors = FALSE
  )
  
  # Transpose for vegan functions (samples as rows)
  abundance_t <- t(abundance_matrix)
  
  # Calculate diversity metrics
  diversity_data$richness <- specnumber(abundance_t)  # Species richness
  diversity_data$shannon <- diversity(abundance_t, index = "shannon")  # Shannon diversity
  diversity_data$simpson <- diversity(abundance_t, index = "simpson")  # Simpson diversity
  diversity_data$evenness <- diversity_data$shannon / log(diversity_data$richness)  # Pielou's evenness
  diversity_data$total_abundance <- rowSums(abundance_t)  # Total viral abundance
  
  # Calculate dominance (relative abundance of most abundant ORF)
  diversity_data$dominance <- apply(abundance_t, 1, function(x) max(x) / sum(x))
  
  # Calculate viral load coefficient of variation
  diversity_data$viral_load_cv <- apply(abundance_t, 1, function(x) {
    if(mean(x) == 0) return(0)
    sd(x) / mean(x)
  })
  
  # Handle infinite/NaN values
  diversity_data$evenness[is.infinite(diversity_data$evenness) | is.nan(diversity_data$evenness)] <- 0
  diversity_data$viral_load_cv[is.infinite(diversity_data$viral_load_cv) | is.nan(diversity_data$viral_load_cv)] <- 0
  
  # Merge with metadata
  diversity_full <- merge(diversity_data, metadata, by = "sample_id")
  
  return(diversity_full)
}

# Calculate diversity metrics
diversity_data_full <- calculate_diversity_metrics(orf_data_matched, metadata_clean)

cat("Diversity metrics calculated for", nrow(diversity_data_full), "samples\n")
cat("Metrics calculated: richness, shannon, simpson, evenness, total_abundance, dominance, viral_load_cv\n")

# Save diversity data
write.csv(diversity_data_full, "diversity_data_full.csv", row.names = FALSE)

# ============================================================================
# PART 3: DIVERSITY TRAJECTORY ANALYSIS
# ============================================================================

cat("\n3. Performing diversity trajectory analysis...\n")

# Define diversity metrics to analyze
diversity_metrics <- c("richness", "shannon", "simpson", "evenness", "total_abundance", "dominance", "viral_load_cv")

# Function to perform limma analysis for diversity trajectories
perform_diversity_trajectory_analysis <- function(diversity_data, metrics) {
  
  results_list <- list()
  
  for(metric in metrics) {
    cat("  Analyzing", metric, "...\n")
    
    # Create design matrix
    design_formula <- as.formula("~ Dx.Status * onset_timeline_numeric + Country + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode")
    
    # Handle missing values
    complete_data <- diversity_data[complete.cases(diversity_data[, c(metric, "Dx.Status", "onset_timeline_numeric", "Country", "Sex", "Age.at.Gluten.Introduction..months.", "HLA.Category", "feeding_first_year", "Delivery.Mode")]), ]
    
    if(nrow(complete_data) < 10) {
      cat("    Skipping", metric, "- insufficient data\n")
      next
    }
    
    # Create design matrix
    design <- model.matrix(design_formula, data = complete_data)
    
    # Create response vector
    response <- complete_data[[metric]]
    
    # Fit limma model with blocking by patient
    if("patientID" %in% colnames(complete_data)) {
      patient_factor <- factor(complete_data$patientID)
      # Use duplicateCorrelation for repeated measures
      corfit <- duplicateCorrelation(matrix(response, nrow = 1), design, block = patient_factor)
      fit <- lmFit(matrix(response, nrow = 1), design, block = patient_factor, correlation = corfit$consensus)
    } else {
      fit <- lmFit(matrix(response, nrow = 1), design)
    }
    
    # Apply empirical Bayes
    fit <- eBayes(fit)
    
    # Extract results for interaction term
    coef_names <- colnames(design)
    interaction_coef <- grep("Dx.Status.*onset_timeline_numeric", coef_names, value = TRUE)
    
    if(length(interaction_coef) > 0) {
      results <- topTable(fit, coef = interaction_coef, number = 1, sort.by = "none")
      results$metric <- metric
      results$n_samples <- nrow(complete_data)
      results$n_patients <- length(unique(complete_data$patientID))
      
      results_list[[metric]] <- results
    }
  }
  
  # Combine results
  if(length(results_list) > 0) {
    combined_results <- do.call(rbind, results_list)
    return(combined_results)
  } else {
    return(NULL)
  }
}

# Perform trajectory analysis
trajectory_results <- perform_diversity_trajectory_analysis(diversity_data_full, diversity_metrics)

if(!is.null(trajectory_results)) {
  cat("Trajectory analysis completed for", nrow(trajectory_results), "metrics\n")
  write.csv(trajectory_results, "diversity_trajectory_results.csv", row.names = FALSE)
  
  # Print significant results
  sig_results <- trajectory_results[trajectory_results$adj.P.Val < 0.05, ]
  if(nrow(sig_results) > 0) {
    cat("Significant trajectory interactions (adj.p < 0.05):\n")
    print(sig_results[, c("metric", "logFC", "P.Value", "adj.P.Val")])
  } else {
    cat("No significant trajectory interactions found\n")
  }
} else {
  cat("No trajectory results could be generated\n")
}

# ============================================================================
# PART 4: SLOPE ANALYSIS
# ============================================================================

cat("\n4. Performing slope analysis...\n")

# Function to calculate patient-level slopes
calculate_patient_slopes <- function(diversity_data, metrics) {
  
  slope_data <- data.frame()
  
  # Get unique patients
  patients <- unique(diversity_data$patientID)
  
  for(patient in patients) {
    patient_data <- diversity_data[diversity_data$patientID == patient, ]
    
    # Need at least 3 timepoints for slope calculation
    if(nrow(patient_data) < 3) next
    
    # Get patient metadata (constant across timepoints)
    patient_info <- patient_data[1, c("patientID", "Dx.Status", "Country", "Sex", "Age.at.Gluten.Introduction..months.", "HLA.Category", "feeding_first_year", "Delivery.Mode")]
    
    # Calculate slopes for each metric
    patient_slopes <- patient_info
    
    for(metric in metrics) {
      # Fit linear model: metric ~ onset_timeline_numeric
      if(sum(!is.na(patient_data[[metric]])) >= 3) {
        lm_result <- lm(as.formula(paste(metric, "~ onset_timeline_numeric")), data = patient_data)
        patient_slopes[[paste0(metric, "_slope")]] <- coef(lm_result)["onset_timeline_numeric"]
      } else {
        patient_slopes[[paste0(metric, "_slope")]] <- NA
      }
    }
    
    slope_data <- rbind(slope_data, patient_slopes)
  }
  
  return(slope_data)
}

# Calculate patient slopes
slope_data <- calculate_patient_slopes(diversity_data_full, diversity_metrics)

cat("Calculated slopes for", nrow(slope_data), "patients\n")
write.csv(slope_data, "slope_data.csv", row.names = FALSE)

# Function to analyze slope differences between groups
analyze_slope_differences <- function(slope_data, metrics) {
  
  results_list <- list()
  
  for(metric in metrics) {
    slope_col <- paste0(metric, "_slope")
    
    if(!slope_col %in% colnames(slope_data)) next
    
    cat("  Analyzing", metric, "slope differences...\n")
    
    # Remove missing values
    complete_data <- slope_data[!is.na(slope_data[[slope_col]]), ]
    
    if(nrow(complete_data) < 10) {
      cat("    Skipping", metric, "- insufficient data\n")
      next
    }
    
    # Create design matrix for slope analysis
    design_formula <- as.formula(paste("~", "Dx.Status + Country + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode"))
    design <- model.matrix(design_formula, data = complete_data)
    
    # Create response vector
    response <- complete_data[[slope_col]]
    
    # Fit limma model
    fit <- lmFit(matrix(response, nrow = 1), design)
    fit <- eBayes(fit)
    
    # Extract results for disease status
    dx_coef <- grep("Dx.Status", colnames(design), value = TRUE)
    
    if(length(dx_coef) > 0) {
      results <- topTable(fit, coef = dx_coef[1], number = 1, sort.by = "none")
      results$metric <- metric
      results$n_patients <- nrow(complete_data)
      
      # Add group means
      group_means <- aggregate(complete_data[[slope_col]], by = list(complete_data$Dx.Status), mean, na.rm = TRUE)
      celiac_mean <- group_means$x[group_means$Group.1 == "CELIAC"]
      control_mean <- group_means$x[group_means$Group.1 == "CONTROL"]
      
      results$celiac_mean_slope <- ifelse(length(celiac_mean) > 0, celiac_mean, NA)
      results$control_mean_slope <- ifelse(length(control_mean) > 0, control_mean, NA)
      
      results_list[[metric]] <- results
    }
  }
  
  # Combine results
  if(length(results_list) > 0) {
    combined_results <- do.call(rbind, results_list)
    return(combined_results)
  } else {
    return(NULL)
  }
}

# Analyze slope differences
slope_results <- analyze_slope_differences(slope_data, diversity_metrics)

if(!is.null(slope_results)) {
  cat("Slope analysis completed for", nrow(slope_results), "metrics\n")
  write.csv(slope_results, "slope_analysis_results.csv", row.names = FALSE)
  
  # Print significant results
  sig_slope_results <- slope_results[slope_results$adj.P.Val < 0.05, ]
  if(nrow(sig_slope_results) > 0) {
    cat("Significant slope differences (adj.p < 0.05):\n")
    print(sig_slope_results[, c("metric", "logFC", "P.Value", "adj.P.Val", "celiac_mean_slope", "control_mean_slope")])
  } else {
    cat("No significant slope differences found\n")
  }
} else {
  cat("No slope results could be generated\n")
}

cat("\n=== COMPOSITIONAL ANALYSIS PART 1 COMPLETED ===\n")
cat("Generated files:\n")
cat("- diversity_data_full.csv: Sample-level diversity metrics\n")
cat("- diversity_trajectory_results.csv: Trajectory interaction results\n")
cat("- slope_data.csv: Patient-level slope data\n")
cat("- slope_analysis_results.csv: Slope difference analysis\n")

cat("\nNext: Run stability and turnover analysis, then visualizations\n")
cat("Analysis completed at:", format(Sys.time()), "\n")