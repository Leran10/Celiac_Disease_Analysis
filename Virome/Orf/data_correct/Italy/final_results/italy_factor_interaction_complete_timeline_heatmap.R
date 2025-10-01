#!/usr/bin/env Rscript
# Italy Complete Timeline Factor Model Interaction Heatmap
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Create heatmap showing all timepoints for interaction-significant ORFs (Italy cohort)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== ITALY COMPLETE TIMELINE FACTOR MODEL INTERACTION HEATMAP ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD OVERLAPPING ORFS FROM ITALY ANALYSIS
# ============================================================================

cat("1. Loading overlapping ORFs from Italy interaction analysis...\n")

# Load the union ORFs from Italy factor analysis
union_orfs <- read.csv("italy_factor_interaction_union_orfs.csv")
orf_list <- union_orfs$ORF

cat("Loaded", length(orf_list), "ORFs significant in Italy interaction timepoints\n")

# ============================================================================
# PART 2: LOAD DATA AND RECREATE FACTOR MODEL
# ============================================================================

cat("2. Loading data and fitting factor model...\n")

# Load abundance data
abundance_data_raw <- read.csv("Italy.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names (remove X prefix if present)
sample_names <- colnames(abundance_data)
sample_names_clean <- gsub("^X", "", sample_names)
colnames(abundance_data) <- sample_names_clean

# Load metadata
metadata_raw <- read.csv("Italy.metadata.clean.csv")
rownames(metadata_raw) <- metadata_raw[, 1]
metadata <- metadata_raw

# Match samples
matching_samples <- intersect(colnames(abundance_data), rownames(metadata))
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")

# Set reference levels with reversed timeline order (same as US analysis)
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))
metadata_filtered$onset_timeline_combined <- factor(metadata_filtered$onset_timeline_combined, 
                                                   levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0"))

# Set other factor levels
metadata_filtered$feeding_first_year <- factor(metadata_filtered$feeding_first_year,
                                               levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
metadata_filtered$HLA.Category <- factor(metadata_filtered$HLA.Category,
                                         levels = c("Standard Risk","High Risk","Low/No Risk"))
metadata_filtered$Sex <- factor(metadata_filtered$Sex,
                                levels = c("Female","Male"))
metadata_filtered$Delivery.Mode <- factor(metadata_filtered$Delivery.Mode,
                                          levels = c("Vaginal","C-Section"))
metadata_filtered$Age.at.Gluten.Introduction..months. <- as.numeric(metadata_filtered$Age.at.Gluten.Introduction..months.)

cat("Timeline levels:", levels(metadata_filtered$onset_timeline_combined), "\n")

# Create design matrix with FACTOR timeline
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

cat("Design matrix dimensions:", nrow(linear.model.design), "x", ncol(linear.model.design), "\n")

# Convert abundance data to matrix for limma
abundance_matrix <- as.matrix(abundance_filtered)

# Fit limma model with voomLmFit and patient blocking
cat("Fitting factor model with patient blocking...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get fitted values from the factor model
fitted_values <- fitted(linear.model.fit)

cat("Factor model fitted successfully\n")

# ============================================================================
# PART 3: CREATE COMPLETE TIMELINE HEATMAP MATRIX
# ============================================================================

cat("3. Creating complete timeline heatmap matrix...\n")

# Get all timepoints
all_timepoints <- levels(metadata_filtered$onset_timeline_combined)
cat("All timepoints:", paste(all_timepoints, collapse = ", "), "\n")

# Create matrix for all timepoints
heatmap_matrix_complete <- matrix(NA, 
                                 nrow = length(orf_list), 
                                 ncol = length(all_timepoints),
                                 dimnames = list(orf_list, all_timepoints))

cat("Processing", length(orf_list), "ORFs across", length(all_timepoints), "timepoints...\n")

for(i in 1:length(orf_list)) {
  orf_id <- orf_list[i]
  
  if(orf_id %in% rownames(fitted_values)) {
    orf_fitted <- fitted_values[orf_id, ]
    
    for(j in 1:length(all_timepoints)) {
      timepoint <- all_timepoints[j]
      
      # Get samples at this timepoint
      timepoint_samples <- metadata_filtered$onset_timeline_combined == timepoint
      
      if(sum(timepoint_samples) > 0) {
        timepoint_data <- data.frame(
          fitted_abundance = as.numeric(orf_fitted[timepoint_samples]),
          Dx.Status = metadata_filtered$Dx.Status[timepoint_samples]
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
            # Calculate abundance difference: control - celiac
            heatmap_matrix_complete[i, j] <- control_fitted - celiac_fitted
          }
        }
      }
    }
  }
  
  if(i %% 500 == 0) {
    cat("Processed", i, "of", length(orf_list), "ORFs\n")
  }
}

# Remove ORFs with too many missing values
valid_orfs <- rowSums(!is.na(heatmap_matrix_complete)) >= (ncol(heatmap_matrix_complete) * 0.5)
heatmap_matrix_clean <- heatmap_matrix_complete[valid_orfs, ]

cat("Final matrix:", nrow(heatmap_matrix_clean), "ORFs x", ncol(heatmap_matrix_clean), "timepoints\n")
cat("Data range:", round(range(heatmap_matrix_clean, na.rm = TRUE), 2), "\n")

# ============================================================================
# PART 4: CREATE COMPLETE TIMELINE HEATMAP
# ============================================================================

cat("4. Creating complete timeline heatmap...\n")

if(nrow(heatmap_matrix_clean) > 0) {
  # Set up color scheme
  data_range <- range(heatmap_matrix_clean, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # Create the complete timeline heatmap
  png("italy_factor_interaction_complete_timeline_heatmap.png", 
      width = 2400, height = max(800, nrow(heatmap_matrix_clean) * 12), res = 150)
  
  pheatmap(heatmap_matrix_clean,
           cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
           color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Italy Cohort: Factor Model -", nrow(heatmap_matrix_clean), 
                       "Interaction-Significant ORFs Across Complete Timeline\n",
                       "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 18, 
           fontsize_row = max(6, min(10, 800/nrow(heatmap_matrix_clean))), 
           fontsize_col = 20,
           cex_main = 1.4,
           show_rownames = FALSE,  # Too many to show individual names clearly
           show_colnames = TRUE, 
           border_color = NA)
  
  dev.off()
  
  cat("Complete timeline heatmap created successfully!\n")
  
  # Also create a smaller version with top variable ORFs for better visualization
  if(nrow(heatmap_matrix_clean) > 200) {
    orf_variances <- apply(heatmap_matrix_clean, 1, function(x) var(x, na.rm = TRUE))
    orf_variances <- orf_variances[!is.na(orf_variances)]
    n_top <- min(200, length(orf_variances))
    top_variable_orfs <- names(sort(orf_variances, decreasing = TRUE)[1:n_top])
    top_heatmap_matrix <- heatmap_matrix_clean[top_variable_orfs, ]
    
    png("italy_factor_interaction_complete_timeline_top200_heatmap.png", 
        width = 2400, height = max(800, n_top * 15), res = 150)
    
    pheatmap(top_heatmap_matrix,
             cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
             color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
             main = paste("Italy Cohort: Factor Model - Top", n_top, 
                         "Most Variable Interaction-Significant ORFs\n",
                         "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
             fontsize = 18, 
             fontsize_row = max(8, min(12, 1000/n_top)), 
             fontsize_col = 20,
             cex_main = 1.4,
             show_rownames = TRUE,
             show_colnames = TRUE, 
             border_color = NA)
    
    dev.off()
    
    cat("Top", n_top, "variable ORFs heatmap also created\n")
  }
}

# ============================================================================
# PART 5: SAVE DATA AND CREATE SUMMARY
# ============================================================================

cat("5. Saving data and creating summary...\n")

# Save the complete timeline matrix
write.csv(heatmap_matrix_clean, "italy_factor_interaction_complete_timeline_matrix.csv")

# Create summary statistics for each timepoint
timepoint_summary <- data.frame(
  Timepoint = colnames(heatmap_matrix_clean),
  Valid_ORFs = colSums(!is.na(heatmap_matrix_clean)),
  Mean_Difference = round(colMeans(heatmap_matrix_clean, na.rm = TRUE), 3),
  SD_Difference = round(apply(heatmap_matrix_clean, 2, sd, na.rm = TRUE), 3),
  Range_Min = round(apply(heatmap_matrix_clean, 2, min, na.rm = TRUE), 3),
  Range_Max = round(apply(heatmap_matrix_clean, 2, max, na.rm = TRUE), 3),
  stringsAsFactors = FALSE
)

write.csv(timepoint_summary, "italy_factor_interaction_complete_timeline_summary.csv", row.names = FALSE)

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== ITALY COMPLETE TIMELINE FACTOR MODEL HEATMAP COMPLETED ===\n")
cat("Generated files:\n")
cat("- italy_factor_interaction_complete_timeline_heatmap.png: Complete timeline heatmap with", nrow(heatmap_matrix_clean), "ORFs\n")

if(exists("n_top")) {
  cat("- italy_factor_interaction_complete_timeline_top200_heatmap.png: Top", n_top, "most variable ORFs\n")
}

cat("- italy_factor_interaction_complete_timeline_matrix.csv: Complete heatmap matrix data\n")
cat("- italy_factor_interaction_complete_timeline_summary.csv: Summary statistics by timepoint\n")

cat("\nSummary by timepoint:\n")
print(timepoint_summary)

cat("\nKey Features:\n")
cat("- ORFs included:", nrow(heatmap_matrix_clean), "(interaction-significant from t0-24, t0-30)\n")
cat("- Timepoints shown:", ncol(heatmap_matrix_clean), "(complete timeline)\n")
cat("- Data range:", round(min(heatmap_matrix_clean, na.rm = TRUE), 2), "to", round(max(heatmap_matrix_clean, na.rm = TRUE), 2), "\n")

# Identify timepoints with strongest patterns
strongest_positive <- timepoint_summary$Timepoint[which.max(timepoint_summary$Mean_Difference)]
strongest_negative <- timepoint_summary$Timepoint[which.min(timepoint_summary$Mean_Difference)]

cat("- Strongest CONTROL dominance:", strongest_positive, "(mean difference =", timepoint_summary$Mean_Difference[timepoint_summary$Timepoint == strongest_positive], ")\n")
cat("- Strongest CELIAC dominance:", strongest_negative, "(mean difference =", timepoint_summary$Mean_Difference[timepoint_summary$Timepoint == strongest_negative], ")\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")