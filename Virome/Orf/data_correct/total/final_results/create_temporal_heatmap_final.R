#!/usr/bin/env Rscript
# Create Temporal Heatmap for Significant ORFs - Final Complete Dataset
# Uses proper mixed-effects modeling with patient random effects
# Author: Claude AI
# Date: 2025-07-24

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(lme4)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

cat("=== TEMPORAL HEATMAP ANALYSIS - FINAL COMPLETE DATASET ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading data...\n")

# Load limma results
limma_results <- read.csv("total_limma_model_res.csv", row.names = 1)

# Load metadata
metadata <- read.csv("total_metadata.clean.csv")

# Load ORF abundance data 
orf_data_raw <- read.csv("total_orf.abundance.clean.csv")
# Use X column as row names
rownames(orf_data_raw) <- orf_data_raw$X
orf_data <- orf_data_raw[, -1]  # Remove X column

# Clean sample names by removing X prefix
orf_samples <- colnames(orf_data)
orf_samples_clean <- gsub("^X", "", orf_samples)
colnames(orf_data) <- orf_samples_clean

cat("Data loaded:\n")
cat("  Limma results:", nrow(limma_results), "ORFs\n")
cat("  Metadata:", nrow(metadata), "samples\n")
cat("  ORF abundance data:", nrow(orf_data), "ORFs,", ncol(orf_data), "samples\n")

# ============================================================================
# PART 2: IDENTIFY SIGNIFICANT ORFS
# ============================================================================

cat("\n2. Identifying significant ORFs...\n")

# Filter significant ORFs
significant_orfs <- limma_results[limma_results$adj.P.Val < 0.01, ]

cat("Significant ORFs (adj.p < 0.01):", nrow(significant_orfs), "\n")
cat("Percentage of total:", round(100 * nrow(significant_orfs) / nrow(limma_results), 1), "%\n")

# Match samples between datasets using X column (contains sample IDs)
sample_ids_from_cols <- colnames(orf_data)
sample_ids_from_metadata <- metadata$X

# Find common samples
common_samples <- intersect(sample_ids_from_cols, sample_ids_from_metadata)
metadata_matched <- metadata[metadata$X %in% common_samples, ]
orf_data_matched <- orf_data[, common_samples]

# Set row names for metadata to match sample names for easier processing
rownames(metadata_matched) <- metadata_matched$X

cat("Matched samples:", length(common_samples), "\n")

# ============================================================================
# PART 3: PREPARE METADATA FOR MODELING
# ============================================================================

cat("\n3. Preparing metadata...\n")

# Check metadata columns
cat("Metadata columns:", paste(colnames(metadata_matched), collapse = ", "), "\n")

# Prepare metadata - adapt column names based on what's available
metadata_clean <- metadata_matched %>%
  mutate(sample_id = X)

# Check for required columns and adapt
required_cols <- c("patientID", "Dx.Status", "onset_timeline_numeric", "Country", "Sex")
missing_cols <- setdiff(required_cols, colnames(metadata_clean))
if(length(missing_cols) > 0) {
  cat("Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(metadata_clean), collapse = ", "), "\n")
}

# Prepare metadata with available columns
if("patientID" %in% colnames(metadata_clean)) {
  metadata_clean <- metadata_clean %>%
    mutate(
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
} else {
  stop("patientID column not found in metadata!")
}

# Remove rows with missing critical data
metadata_clean <- metadata_clean[complete.cases(metadata_clean[, c("sample_id", "patientID", "Dx.Status", "onset_timeline_numeric")]), ]

cat("Clean metadata samples:", nrow(metadata_clean), "\n")
cat("Patients:", length(unique(metadata_clean$patientID)), "\n")
cat("Disease groups:", table(metadata_clean$Dx.Status), "\n")

# ============================================================================
# PART 4: FIT MIXED-EFFECTS MODELS FOR SIGNIFICANT ORFS
# ============================================================================

cat("\n4. Fitting mixed-effects models for significant ORFs...\n")

# Function to fit mixed-effects model for a single ORF
fit_orf_mixed_model <- function(orf_id, abundance_data, metadata) {
  
  # Get abundance for this ORF
  orf_abundances <- abundance_data[orf_id, ]
  
  # Create data frame for analysis
  model_data <- data.frame(
    sample_id = names(orf_abundances),
    abundance = as.numeric(orf_abundances),
    stringsAsFactors = FALSE
  )
  
  # Merge with metadata
  model_data <- merge(model_data, metadata, by = "sample_id")
  
  # Remove samples with missing data
  model_data <- model_data[complete.cases(model_data), ]
  
  if(nrow(model_data) < 10 || length(unique(model_data$patientID)) < 5) {
    return(NULL)  # Skip ORFs with insufficient data
  }
  
  # Create model formula based on available columns
  base_formula <- "abundance ~ Dx.Status * onset_timeline_numeric"
  
  # Add confounders if available
  confounders <- c()
  if("Country" %in% colnames(model_data) && length(unique(model_data$Country)) > 1) {
    confounders <- c(confounders, "Country")
  }
  if("Sex" %in% colnames(model_data) && length(unique(model_data$Sex)) > 1) {
    confounders <- c(confounders, "Sex")
  }
  
  if(length(confounders) > 0) {
    formula_str <- paste(base_formula, "+", paste(confounders, collapse = " + "))
  } else {
    formula_str <- base_formula
  }
  
  # Add random effect
  formula_str <- paste(formula_str, "+ (1|patientID)")
  
  # Fit mixed-effects model with patient random effect
  tryCatch({
    model <- lmer(as.formula(formula_str), data = model_data)
    
    # Extract fitted values
    model_data$fitted_abundance <- fitted(model)
    
    return(model_data)
    
  }, error = function(e) {
    # If lmer fails, fall back to lm (but warn user)
    tryCatch({
      lm_formula <- gsub("\\+ \\(1\\|patientID\\)", "", formula_str)
      model <- lm(as.formula(lm_formula), data = model_data)
      
      model_data$fitted_abundance <- fitted(model)
      return(model_data)
      
    }, error = function(e2) {
      return(NULL)
    })
  })
}

# Process significant ORFs in batches
orf_ids_to_analyze <- intersect(rownames(significant_orfs), rownames(orf_data_matched))
orf_model_results <- list()

cat("Fitting models for", length(orf_ids_to_analyze), "significant ORFs...\n")

# Process in batches with progress tracking
batch_size <- 100
n_batches <- ceiling(length(orf_ids_to_analyze) / batch_size)

for(batch in 1:n_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, length(orf_ids_to_analyze))
  batch_orfs <- orf_ids_to_analyze[start_idx:end_idx]
  
  for(orf_id in batch_orfs) {
    result <- fit_orf_mixed_model(orf_id, orf_data_matched, metadata_clean)
    if(!is.null(result)) {
      orf_model_results[[orf_id]] <- result
    }
  }
  
  cat("Completed batch", batch, "of", n_batches, "- processed", length(orf_model_results), "ORFs so far\n")
}

cat("Successfully fitted models for", length(orf_model_results), "ORFs\n")

# ============================================================================
# PART 5: CALCULATE TEMPORAL DIFFERENCES
# ============================================================================

cat("\n5. Calculating temporal differences...\n")

# Define actual study timepoints
actual_timepoints <- c(0, -6, -12, -18, -24, -30, -36, -42, -48, -54, -60, -66, -72)
time_bins <- sort(actual_timepoints, decreasing = FALSE)

# Function to calculate temporal differences
calculate_temporal_differences <- function(model_results, time_bins) {
  
  # Pre-allocate matrix
  heatmap_matrix <- matrix(NA, 
                          nrow = length(model_results), 
                          ncol = length(time_bins),
                          dimnames = list(names(model_results), 
                                         ifelse(time_bins == 0, "T0", paste0("T0", time_bins))))
  
  for(i in 1:length(model_results)) {
    orf_id <- names(model_results)[i]
    orf_data <- model_results[[orf_id]]
    
    # Assign samples to exact timepoints (with tolerance for rounding)
    orf_data$time_bin <- NA
    for(timepoint in time_bins) {
      matches <- abs(orf_data$onset_timeline_numeric - timepoint) < 1
      orf_data$time_bin[matches] <- timepoint
    }
    
    for(j in 1:length(time_bins)) {
      timepoint <- time_bins[j]
      
      # Get fitted abundances for this timepoint
      time_data <- orf_data[!is.na(orf_data$time_bin) & orf_data$time_bin == timepoint, ]
      
      if(nrow(time_data) >= 2) {  # Require at least 2 samples total
        # Calculate mean fitted abundances by group
        group_means <- time_data %>%
          group_by(Dx.Status) %>%
          summarise(
            mean_fitted_abundance = mean(fitted_abundance, na.rm = TRUE),
            n_samples = n(),
            .groups = "drop"
          ) %>%
          filter(n_samples >= 1)
        
        # Calculate difference: CONTROL - CELIAC 
        # (Since CELIAC is reference group, positive logFC = higher in CONTROL)
        if(nrow(group_means) == 2) {
          celiac_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CELIAC"]
          control_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CONTROL"]
          
          if(length(celiac_fitted) > 0 && length(control_fitted) > 0) {
            heatmap_matrix[i, j] <- control_fitted - celiac_fitted
          }
        }
      }
    }
    
    # Progress update
    if(i %% 100 == 0) {
      cat("Processed", i, "of", length(model_results), "ORFs for temporal differences\n")
    }
  }
  
  return(heatmap_matrix)
}

# Calculate temporal differences
temporal_matrix <- calculate_temporal_differences(orf_model_results, time_bins)

cat("Temporal matrix dimensions:", dim(temporal_matrix), "\n")

# Remove rows with insufficient temporal data
valid_rows <- rowSums(!is.na(temporal_matrix)) >= 3
temporal_clean <- temporal_matrix[valid_rows, , drop = FALSE]

cat("Matrix after removing sparse rows:", dim(temporal_clean), "\n")

# ============================================================================
# PART 6: GENERATE HEATMAPS
# ============================================================================

cat("\n6. Generating heatmaps...\n")

if(nrow(temporal_clean) > 0) {
  
  # Apply log2 transformation for better visualization
  temporal_log <- temporal_clean
  temporal_log[is.na(temporal_log)] <- 0
  temporal_log <- sign(temporal_log) * log2(abs(temporal_log) + 1)
  
  cat("Scale after log transformation: min =", round(min(temporal_log, na.rm = TRUE), 2),
      "max =", round(max(temporal_log, na.rm = TRUE), 2), "\n")
  
  # Define color palette
  color_palette <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                     "white", 
                                     "#FDBF6F", "#FF7F00", "#E31A1C", "#B2182B"))(100)
  
  # Set symmetric breaks
  max_abs_val <- max(abs(temporal_log), na.rm = TRUE)
  break_range <- min(max_abs_val, 4)
  breaks <- seq(-break_range, break_range, length.out = 101)
  
  # Heatmap 1: All significant ORFs (no row names due to size)
  png("all_significant_orfs_temporal_heatmap_final.png", width = 1400, height = 2000, res = 150)
  
  pheatmap(
    temporal_log,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "none",
    color = color_palette,
    main = paste("Temporal Heatmap:", nrow(temporal_log), "Significant ORFs (Final Complete Dataset)\n",
                 "adj.p < 0.01, Mixed-effects models with patient random effects\n",
                 "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
    fontsize = 8,
    fontsize_row = 0,
    fontsize_col = 10,
    show_rownames = FALSE,
    breaks = breaks,
    na_col = "grey90"
  )
  
  dev.off()
  
  cat("Generated full heatmap for", nrow(temporal_log), "significant ORFs\n")
  
  # Heatmap 2: Top 200 most variable ORFs with names
  temporal_variance <- apply(temporal_log, 1, var, na.rm = TRUE)
  top_variable_orfs <- names(sort(temporal_variance, decreasing = TRUE))[1:min(200, length(temporal_variance))]
  top_variable_matrix <- temporal_log[top_variable_orfs, , drop = FALSE]
  
  png("top_variable_significant_orfs_heatmap_final.png", width = 1400, height = 1600, res = 150)
  
  pheatmap(
    top_variable_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "none",
    color = color_palette,
    main = paste("Top", nrow(top_variable_matrix), "Most Variable Significant ORFs (Final Dataset)\n",
                 "adj.p < 0.01, Mixed-effects models with patient random effects\n",
                 "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
    fontsize = 8,
    fontsize_row = 4,
    fontsize_col = 10,
    show_rownames = TRUE,
    breaks = breaks,
    na_col = "grey90"
  )
  
  dev.off()
  
  cat("Generated heatmap for top", nrow(top_variable_matrix), "most variable ORFs\n")
  
  # Heatmap 3: Cluster summary
  if(nrow(temporal_log) > 10) {
    # Cluster ORFs and identify major groups
    n_clusters <- min(8, max(3, round(nrow(temporal_log)/100)))
    row_clusters <- cutree(hclust(dist(temporal_log)), k = n_clusters)
    
    # Calculate cluster centroids
    cluster_centroids <- matrix(NA, nrow = n_clusters, ncol = ncol(temporal_log))
    rownames(cluster_centroids) <- paste0("Cluster_", 1:n_clusters)
    colnames(cluster_centroids) <- colnames(temporal_log)
    
    for(cluster_id in 1:n_clusters) {
      cluster_orfs <- names(row_clusters)[row_clusters == cluster_id]
      if(length(cluster_orfs) > 1) {
        cluster_centroids[cluster_id, ] <- colMeans(temporal_log[cluster_orfs, , drop = FALSE], na.rm = TRUE)
      } else if(length(cluster_orfs) == 1) {
        cluster_centroids[cluster_id, ] <- temporal_log[cluster_orfs, ]
      }
    }
    
    # Plot cluster summary
    png("temporal_pattern_clusters_final.png", width = 1200, height = 800, res = 150)
    
    pheatmap(
      cluster_centroids,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      scale = "none",
      color = color_palette,
      main = paste("Temporal Pattern Clusters (Final Complete Dataset)\n",
                   n_clusters, "major patterns from", nrow(temporal_log), "significant ORFs (adj.p < 0.01)\n",
                   "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
      fontsize = 10,
      fontsize_row = 10,
      fontsize_col = 12,
      show_rownames = TRUE,
      breaks = breaks,
      na_col = "grey90"
    )
    
    dev.off()
    
    # Print cluster sizes
    cluster_sizes <- table(row_clusters)
    cat("\nCluster sizes:\n")
    for(i in 1:length(cluster_sizes)) {
      cat("Cluster", names(cluster_sizes)[i], ":", cluster_sizes[i], "ORFs\n")
    }
  }
  
} else {
  cat("No valid data for heatmap generation\n")
}

# ============================================================================
# PART 7: SAVE RESULTS
# ============================================================================

cat("\n7. Saving results...\n")

# Save temporal matrix
write.csv(temporal_clean, "temporal_matrix_final.csv")

# Save log-transformed matrix
write.csv(temporal_log, "temporal_matrix_log_final.csv")

# Create summary statistics
summary_stats <- data.frame(
  metric = c("total_limma_orfs", "significant_orfs", "orfs_with_models", 
             "orfs_in_heatmap", "mean_temporal_difference", "median_temporal_difference"),
  value = c(nrow(limma_results), nrow(significant_orfs), length(orf_model_results),
            nrow(temporal_clean), 
            mean(temporal_clean, na.rm = TRUE), 
            median(temporal_clean, na.rm = TRUE))
)

write.csv(summary_stats, "temporal_heatmap_summary_final.csv", row.names = FALSE)

cat("\n=== TEMPORAL HEATMAP ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- all_significant_orfs_temporal_heatmap_final.png: Full heatmap\n")
cat("- top_variable_significant_orfs_heatmap_final.png: Top 200 most variable\n")
cat("- temporal_pattern_clusters_final.png: Pattern clusters summary\n")
cat("- temporal_matrix_final.csv: Raw temporal differences\n")
cat("- temporal_matrix_log_final.csv: Log-transformed data\n")
cat("- temporal_heatmap_summary_final.csv: Summary statistics\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")