#!/usr/bin/env Rscript
# Compositional Analysis Part 2 - Stability, Turnover, and Visualizations
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

cat("=== COMPOSITIONAL ANALYSIS PART 2 ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD PREVIOUS RESULTS AND DATA
# ============================================================================

cat("1. Loading previous results and data...\n")

# Load data from part 1
diversity_data_full <- read.csv("diversity_data_full.csv")
slope_data <- read.csv("slope_data.csv")
trajectory_results <- read.csv("diversity_trajectory_results.csv")
slope_results <- read.csv("slope_analysis_results.csv")

# Load original data
orf_data_raw <- read.csv("total_orf.abundance.clean.csv")
rownames(orf_data_raw) <- orf_data_raw$X
orf_data <- orf_data_raw[, -1]
orf_samples_clean <- gsub("^X", "", colnames(orf_data))
colnames(orf_data) <- orf_samples_clean

metadata <- read.csv("total_metadata.clean.csv")

cat("Data loaded successfully\n")

# ============================================================================
# PART 2: STABILITY ANALYSIS
# ============================================================================

cat("\n2. Performing stability analysis...\n")

# Function to calculate patient-level stability metrics
calculate_patient_stability <- function(diversity_data, metrics) {
  
  stability_data <- data.frame()
  
  # Get unique patients
  patients <- unique(diversity_data$patientID)
  
  for(patient in patients) {
    patient_data <- diversity_data[diversity_data$patientID == patient, ]
    
    # Need at least 3 timepoints for stability calculation
    if(nrow(patient_data) < 3) next
    
    # Get patient metadata
    patient_info <- patient_data[1, c("patientID", "Dx.Status", "Country", "Sex", "Age.at.Gluten.Introduction..months.", "HLA.Category", "feeding_first_year", "Delivery.Mode")]
    
    # Calculate stability (coefficient of variation) for each metric
    patient_stability <- patient_info
    
    for(metric in metrics) {
      metric_values <- patient_data[[metric]]
      metric_values <- metric_values[!is.na(metric_values)]
      
      if(length(metric_values) >= 3 && mean(metric_values) != 0) {
        # Coefficient of variation (lower = more stable)
        cv <- sd(metric_values) / mean(metric_values)
        patient_stability[[paste0(metric, "_stability")]] <- cv
      } else {
        patient_stability[[paste0(metric, "_stability")]] <- NA
      }
    }
    
    stability_data <- rbind(stability_data, patient_stability)
  }
  
  return(stability_data)
}

# Calculate stability metrics
diversity_metrics <- c("richness", "shannon", "simpson", "evenness", "total_abundance", "dominance", "viral_load_cv")
stability_data <- calculate_patient_stability(diversity_data_full, diversity_metrics)

cat("Calculated stability for", nrow(stability_data), "patients\n")
write.csv(stability_data, "stability_data.csv", row.names = FALSE)

# Function to analyze stability differences between groups
analyze_stability_differences <- function(stability_data, metrics) {
  
  results_list <- list()
  
  for(metric in metrics) {
    stability_col <- paste0(metric, "_stability")
    
    if(!stability_col %in% colnames(stability_data)) next
    
    cat("  Analyzing", metric, "stability differences...\n")
    
    # Remove missing/infinite values
    complete_data <- stability_data[!is.na(stability_data[[stability_col]]) & 
                                   is.finite(stability_data[[stability_col]]), ]
    
    if(nrow(complete_data) < 10) {
      cat("    Skipping", metric, "- insufficient data\n")
      next
    }
    
    # Create design matrix
    design_formula <- as.formula(paste("~", "Dx.Status + Country + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode"))
    design <- model.matrix(design_formula, data = complete_data)
    
    # Create response vector
    response <- complete_data[[stability_col]]
    
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
      group_means <- aggregate(complete_data[[stability_col]], by = list(complete_data$Dx.Status), mean, na.rm = TRUE)
      celiac_mean <- group_means$x[group_means$Group.1 == "CELIAC"]
      control_mean <- group_means$x[group_means$Group.1 == "CONTROL"]
      
      results$celiac_mean_stability <- ifelse(length(celiac_mean) > 0, celiac_mean, NA)
      results$control_mean_stability <- ifelse(length(control_mean) > 0, control_mean, NA)
      
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

# Analyze stability differences
stability_results <- analyze_stability_differences(stability_data, diversity_metrics)

if(!is.null(stability_results)) {
  cat("Stability analysis completed for", nrow(stability_results), "metrics\n")
  write.csv(stability_results, "stability_analysis_results.csv", row.names = FALSE)
  
  # Print significant results
  sig_stability_results <- stability_results[stability_results$adj.P.Val < 0.05, ]
  if(nrow(sig_stability_results) > 0) {
    cat("Significant stability differences (adj.p < 0.05):\n")
    print(sig_stability_results[, c("metric", "logFC", "P.Value", "adj.P.Val", "celiac_mean_stability", "control_mean_stability")])
  } else {
    cat("No significant stability differences found\n")
  }
} else {
  cat("No stability results could be generated\n")
}

# ============================================================================
# PART 3: TURNOVER ANALYSIS
# ============================================================================

cat("\n3. Performing turnover analysis...\n")

# Function to calculate patient-level turnover metrics
calculate_patient_turnover <- function(orf_data, metadata, min_timepoints = 3) {
  
  turnover_data <- data.frame()
  
  # Get unique patients
  patients <- unique(metadata$patientID)
  
  for(patient in patients) {
    patient_samples <- metadata[metadata$patientID == patient, ]
    
    # Need at least min_timepoints for turnover calculation
    if(nrow(patient_samples) < min_timepoints) next
    
    # Order by timeline
    patient_samples <- patient_samples[order(patient_samples$onset_timeline_numeric), ]
    
    # Get abundance data for this patient
    patient_orf_data <- orf_data[, patient_samples$sample_id, drop = FALSE]
    
    # Calculate turnover between consecutive timepoints
    turnover_values <- c()
    
    for(i in 2:ncol(patient_orf_data)) {
      # Calculate Bray-Curtis dissimilarity between consecutive timepoints
      timepoint_data <- patient_orf_data[, c(i-1, i)]
      # Transpose for vegdist (samples as rows)
      timepoint_data_t <- data.frame(t(timepoint_data))
      
      # Calculate dissimilarity
      dissim <- vegdist(timepoint_data_t, method = "bray")
      turnover_values <- c(turnover_values, as.numeric(dissim))
    }
    
    # Calculate mean turnover for this patient
    if(length(turnover_values) > 0) {
      patient_info <- patient_samples[1, c("patientID", "Dx.Status", "Country", "Sex", "Age.at.Gluten.Introduction..months.", "HLA.Category", "feeding_first_year", "Delivery.Mode")]
      patient_info$mean_turnover <- mean(turnover_values, na.rm = TRUE)
      patient_info$median_turnover <- median(turnover_values, na.rm = TRUE)
      patient_info$n_transitions <- length(turnover_values)
      
      turnover_data <- rbind(turnover_data, patient_info)
    }
  }
  
  return(turnover_data)
}

# Calculate turnover metrics
turnover_data <- calculate_patient_turnover(orf_data, diversity_data_full)

cat("Calculated turnover for", nrow(turnover_data), "patients\n")
write.csv(turnover_data, "turnover_data.csv", row.names = FALSE)

# Function to analyze turnover differences between groups
analyze_turnover_differences <- function(turnover_data) {
  
  results_list <- list()
  turnover_metrics <- c("mean_turnover", "median_turnover")
  
  for(metric in turnover_metrics) {
    cat("  Analyzing", metric, "differences...\n")
    
    # Remove missing values
    complete_data <- turnover_data[!is.na(turnover_data[[metric]]), ]
    
    if(nrow(complete_data) < 10) {
      cat("    Skipping", metric, "- insufficient data\n")
      next
    }
    
    # Create design matrix
    design_formula <- as.formula(paste("~", "Dx.Status + Country + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode"))
    design <- model.matrix(design_formula, data = complete_data)
    
    # Create response vector
    response <- complete_data[[metric]]
    
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
      group_means <- aggregate(complete_data[[metric]], by = list(complete_data$Dx.Status), mean, na.rm = TRUE)
      celiac_mean <- group_means$x[group_means$Group.1 == "CELIAC"]
      control_mean <- group_means$x[group_means$Group.1 == "CONTROL"]
      
      results$celiac_mean_turnover <- ifelse(length(celiac_mean) > 0, celiac_mean, NA)
      results$control_mean_turnover <- ifelse(length(control_mean) > 0, control_mean, NA)
      
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

# Analyze turnover differences
turnover_results <- analyze_turnover_differences(turnover_data)

if(!is.null(turnover_results)) {
  cat("Turnover analysis completed for", nrow(turnover_results), "metrics\n")
  write.csv(turnover_results, "turnover_analysis_results.csv", row.names = FALSE)
  
  # Print significant results
  sig_turnover_results <- turnover_results[turnover_results$adj.P.Val < 0.05, ]
  if(nrow(sig_turnover_results) > 0) {
    cat("Significant turnover differences (adj.p < 0.05):\n")
    print(sig_turnover_results[, c("metric", "logFC", "P.Value", "adj.P.Val", "celiac_mean_turnover", "control_mean_turnover")])
  } else {
    cat("No significant turnover differences found\n")
  }
} else {
  cat("No turnover results could be generated\n")
}

cat("\n=== COMPOSITIONAL ANALYSIS PART 2 COMPLETED ===\n")
cat("Generated files:\n")
cat("- stability_data.csv: Patient-level stability metrics\n")
cat("- stability_analysis_results.csv: Stability difference analysis\n")
cat("- turnover_data.csv: Patient-level turnover metrics\n")
cat("- turnover_analysis_results.csv: Turnover difference analysis\n")

cat("\nNext: Create comprehensive results summary and visualizations\n")
cat("Analysis completed at:", format(Sys.time()), "\n")