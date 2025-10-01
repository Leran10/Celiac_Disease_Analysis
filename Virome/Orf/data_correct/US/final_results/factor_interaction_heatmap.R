#!/usr/bin/env Rscript
# Factor Model Interaction Heatmap - US Cohort
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Create heatmap using overlapping significant ORFs from Dx.Status * factor(timeline) interactions

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== FACTOR MODEL INTERACTION HEATMAP ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD EXISTING INTERACTION RESULTS
# ============================================================================

cat("1. Loading existing interaction results...\n")

# Load interaction results from the 3 timepoints with significant ORFs
t24_results <- read.csv("significant_orfs_Dx.StatusCONTROL:t0-24.csv", row.names = 1)
t36_results <- read.csv("significant_orfs_Dx.StatusCONTROL:t0-36.csv", row.names = 1)
t42_results <- read.csv("significant_orfs_Dx.StatusCONTROL:t0-over42.csv", row.names = 1)

cat("Loaded interaction results:\n")
cat("- t0-24:", nrow(t24_results), "significant ORFs\n")
cat("- t0-36:", nrow(t36_results), "significant ORFs\n") 
cat("- t0-over42:", nrow(t42_results), "significant ORFs\n")

# Get ORF lists
t24_orfs <- rownames(t24_results)
t36_orfs <- rownames(t36_results)
t42_orfs <- rownames(t42_results)

# Find overlapping ORFs across all 3 timepoints
overlap_all <- intersect(intersect(t24_orfs, t36_orfs), t42_orfs)

# Also find pairwise overlaps for comparison
overlap_24_36 <- intersect(t24_orfs, t36_orfs)
overlap_24_42 <- intersect(t24_orfs, t42_orfs)
overlap_36_42 <- intersect(t36_orfs, t42_orfs)

cat("\nOverlap analysis:\n")
cat("- t0-24 ∩ t0-36:", length(overlap_24_36), "ORFs\n")
cat("- t0-24 ∩ t0-over42:", length(overlap_24_42), "ORFs\n")
cat("- t0-36 ∩ t0-over42:", length(overlap_36_42), "ORFs\n")
cat("- All three timepoints:", length(overlap_all), "ORFs\n")

# ============================================================================
# PART 2: FIT FACTOR MODEL AND GET FITTED VALUES
# ============================================================================

cat("\n2. Fitting factor model to get fitted values...\n")

# Load abundance data
abundance_data_raw <- read.csv("US.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names (remove X prefix if present)
sample_names <- colnames(abundance_data)
sample_names_clean <- gsub("^X", "", sample_names)
colnames(abundance_data) <- sample_names_clean

# Load metadata
metadata_raw <- read.csv("US.metadata.clean.csv")
rownames(metadata_raw) <- metadata_raw$X
metadata <- metadata_raw

# Match samples
matching_samples <- intersect(colnames(abundance_data), rownames(metadata))
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")

# Set reference levels (same as original analysis)
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))
metadata_filtered$onset_timeline_combined <- factor(metadata_filtered$onset_timeline_combined, 
                                                   levels = c("t0", "t0-6", "t0-12", "t0-18", "t0-24", "t0-30", "t0-36", "t0-over42"))

# Create design matrix with FACTOR timeline
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

cat("Design matrix dimensions:", nrow(linear.model.design), "x", ncol(linear.model.design), "\n")

# Convert abundance data to matrix for limma
abundance_matrix <- as.matrix(abundance_filtered)

# Fit limma model with voomLmFit and patient blocking
cat("Running factor model with patient blocking...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get fitted values from the factor model
fitted_values <- fitted(linear.model.fit)

cat("Factor model fitted successfully\n")

# ============================================================================
# PART 3: CREATE HEATMAP MATRICES
# ============================================================================

cat("\n3. Creating heatmap matrices for different ORF sets...\n")

# Function to calculate temporal differences using fitted values
calculate_factor_temporal_differences <- function(fitted_values, metadata, orf_list) {
  
  # Get unique timepoints
  timepoints <- levels(metadata$onset_timeline_combined)
  
  # Create matrix
  heatmap_matrix <- matrix(NA, 
                          nrow = length(orf_list), 
                          ncol = length(timepoints),
                          dimnames = list(orf_list, timepoints))
  
  for(i in 1:length(orf_list)) {
    orf_id <- orf_list[i]
    
    if(orf_id %in% rownames(fitted_values)) {
      orf_fitted <- fitted_values[orf_id, ]
      
      for(j in 1:length(timepoints)) {
        timepoint <- timepoints[j]
        
        # Get samples at this timepoint
        timepoint_samples <- metadata$onset_timeline_combined == timepoint
        
        if(sum(timepoint_samples) > 0) {
          timepoint_data <- data.frame(
            fitted_abundance = as.numeric(orf_fitted[timepoint_samples]),
            Dx.Status = metadata$Dx.Status[timepoint_samples]
          )
          
          # Calculate group means
          group_means <- timepoint_data %>%
            group_by(Dx.Status) %>%
            summarise(
              mean_fitted_abundance = mean(fitted_abundance, na.rm = TRUE),
              n_samples = n(),
              .groups = "drop"
            ) %>%
            filter(n_samples >= 1)
          
          if(nrow(group_means) == 2) {
            celiac_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CELIAC"]
            control_fitted <- group_means$mean_fitted_abundance[group_means$Dx.Status == "CONTROL"]
            
            if(length(celiac_fitted) > 0 && length(control_fitted) > 0) {
              # Same calculation as previous heatmaps: control - celiac
              heatmap_matrix[i, j] <- control_fitted - celiac_fitted
            }
          }
        }
      }
    }
    
    if(i %% 100 == 0) {
      cat("Processed", i, "of", length(orf_list), "ORFs\n")
    }
  }
  
  return(heatmap_matrix)
}

# Create heatmap matrices for different ORF sets
cat("Creating heatmap for all 3-timepoint overlapping ORFs...\n")
if(length(overlap_all) > 0) {
  heatmap_all_overlap <- calculate_factor_temporal_differences(fitted_values, metadata_filtered, overlap_all)
  cat("Created matrix for", length(overlap_all), "overlapping ORFs\n")
}

cat("Creating heatmap for t0-36 & t0-over42 overlapping ORFs...\n")
if(length(overlap_36_42) > 10) {  # Only if reasonable number
  heatmap_36_42 <- calculate_factor_temporal_differences(fitted_values, metadata_filtered, overlap_36_42)
  cat("Created matrix for", length(overlap_36_42), "t0-36 & t0-over42 overlapping ORFs\n")
}

# ============================================================================
# PART 4: CREATE HEATMAPS
# ============================================================================

cat("\n4. Creating factor model interaction heatmaps...\n")

# Set up color scheme
create_heatmap <- function(matrix_data, title_text, filename) {
  if(nrow(matrix_data) > 0) {
    # Remove rows with all NA
    valid_rows <- rowSums(!is.na(matrix_data)) > 0
    clean_matrix <- matrix_data[valid_rows, ]
    
    if(nrow(clean_matrix) > 0) {
      data_range <- range(clean_matrix, na.rm = TRUE)
      max_abs <- max(abs(data_range), na.rm = TRUE)
      color_limits <- c(-max_abs, max_abs)
      colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
      
      png(filename, width = 2200, height = max(800, nrow(clean_matrix) * 15), res = 150)
      pheatmap(clean_matrix,
               cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
               color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
               main = paste(title_text, "\nRed = Higher in CONTROL, Blue = Higher in CELIAC"),
               fontsize = 18, 
               fontsize_row = max(8, min(12, 1000/nrow(clean_matrix))), 
               fontsize_col = 20,
               cex_main = 1.5,
               show_rownames = ifelse(nrow(clean_matrix) <= 50, TRUE, FALSE),
               show_colnames = TRUE, 
               border_color = NA)
      dev.off()
      
      cat("Created heatmap:", filename, "with", nrow(clean_matrix), "ORFs\n")
      return(clean_matrix)
    }
  }
  return(NULL)
}

# Create heatmaps
final_matrices <- list()

if(length(overlap_all) > 0) {
  final_matrices$all_overlap <- create_heatmap(heatmap_all_overlap, 
    paste("US Cohort: Factor Model -", length(overlap_all), "ORFs Significant in All 3 Interaction Timepoints"),
    "factor_interaction_all_overlap_heatmap.png")
}

if(exists("heatmap_36_42") && nrow(heatmap_36_42) > 0) {
  final_matrices$overlap_36_42 <- create_heatmap(heatmap_36_42, 
    paste("US Cohort: Factor Model -", length(overlap_36_42), "ORFs Significant in t0-36 & t0-over42"),
    "factor_interaction_36_42_overlap_heatmap.png")
}

# ============================================================================
# PART 5: SAVE RESULTS AND SUMMARY
# ============================================================================

cat("\n5. Saving results and creating summary...\n")

# Save overlap information
overlap_summary <- data.frame(
  Comparison = c("t0-24 only", "t0-36 only", "t0-over42 only", 
                 "t0-24 ∩ t0-36", "t0-24 ∩ t0-over42", "t0-36 ∩ t0-over42",
                 "All three timepoints"),
  Number_of_ORFs = c(length(setdiff(setdiff(t24_orfs, t36_orfs), t42_orfs)),
                      length(setdiff(setdiff(t36_orfs, t24_orfs), t42_orfs)),
                      length(setdiff(setdiff(t42_orfs, t24_orfs), t36_orfs)),
                      length(overlap_24_36),
                      length(overlap_24_42), 
                      length(overlap_36_42),
                      length(overlap_all)),
  stringsAsFactors = FALSE
)

write.csv(overlap_summary, "factor_interaction_overlap_summary.csv", row.names = FALSE)

# Save overlapping ORF lists
if(length(overlap_all) > 0) {
  write.csv(data.frame(ORF = overlap_all), "all_timepoints_overlapping_orfs.csv", row.names = FALSE)
}

if(length(overlap_36_42) > 0) {
  write.csv(data.frame(ORF = overlap_36_42), "t36_t42_overlapping_orfs.csv", row.names = FALSE)
}

# Save heatmap matrices
if(!is.null(final_matrices$all_overlap)) {
  write.csv(final_matrices$all_overlap, "factor_interaction_all_overlap_matrix.csv")
}

if(!is.null(final_matrices$overlap_36_42)) {
  write.csv(final_matrices$overlap_36_42, "factor_interaction_36_42_matrix.csv")
}

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== FACTOR MODEL INTERACTION HEATMAP ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")

if(length(overlap_all) > 0) {
  cat("- factor_interaction_all_overlap_heatmap.png: Heatmap of", length(overlap_all), "ORFs significant in all 3 timepoints\n")
  cat("- all_timepoints_overlapping_orfs.csv: List of overlapping ORFs\n")
  cat("- factor_interaction_all_overlap_matrix.csv: Heatmap matrix data\n")
}

if(length(overlap_36_42) > 0) {
  cat("- factor_interaction_36_42_overlap_heatmap.png: Heatmap of", length(overlap_36_42), "ORFs in t0-36 & t0-over42\n")
  cat("- t36_t42_overlapping_orfs.csv: List of t0-36 & t0-over42 overlapping ORFs\n")
  cat("- factor_interaction_36_42_matrix.csv: Heatmap matrix data\n")
}

cat("- factor_interaction_overlap_summary.csv: Summary of all overlaps\n")

cat("\nOverlap Summary:\n")
print(overlap_summary)

cat("\nKey Insights:\n")
if(length(overlap_all) > 0) {
  cat("- Found", length(overlap_all), "ORFs with consistent interaction effects across all 3 significant timepoints\n")
  cat("- These represent the most robust temporal biomarkers with interaction patterns\n")
}

if(length(overlap_36_42) > length(overlap_all)) {
  cat("- Found", length(overlap_36_42), "ORFs significant in both t0-36 and t0-over42 (later timepoints)\n")
  cat("- These may represent late-stage interaction effects\n")
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")