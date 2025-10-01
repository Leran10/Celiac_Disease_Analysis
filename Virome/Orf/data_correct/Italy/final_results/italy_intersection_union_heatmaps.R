#!/usr/bin/env Rscript
# Italy Intersection and Union Heatmaps
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Create heatmaps for intersection (AND) and union (OR) of significant timepoints

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== ITALY INTERSECTION AND UNION HEATMAPS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD INDIVIDUAL TIMEPOINT RESULTS
# ============================================================================

cat("1. Loading individual timepoint results...\n")

# Load t0-24 significant ORFs
t0_24_orfs <- read.csv("italy_significant_orfs_Dx_StatusCONTROL_onset_timeline_combinedt0-24.csv")
t0_24_orf_list <- t0_24_orfs$ORF

# Load t0-30 significant ORFs  
t0_30_orfs <- read.csv("italy_significant_orfs_Dx_StatusCONTROL_onset_timeline_combinedt0-30.csv")
t0_30_orf_list <- t0_30_orfs$ORF

cat("t0-24 significant ORFs:", length(t0_24_orf_list), "\n")
cat("t0-30 significant ORFs:", length(t0_30_orf_list), "\n")

# ============================================================================
# PART 2: CALCULATE INTERSECTION AND UNION
# ============================================================================

cat("2. Calculating intersection and union...\n")

# Calculate intersection (ORFs significant in BOTH timepoints)
intersection_orfs <- intersect(t0_24_orf_list, t0_30_orf_list)

# Calculate union (ORFs significant in EITHER timepoint)
union_orfs <- union(t0_24_orf_list, t0_30_orf_list)

cat("Intersection (t0-24 AND t0-30):", length(intersection_orfs), "ORFs\n")
cat("Union (t0-24 OR t0-30):", length(union_orfs), "ORFs\n")

# Verify union matches previous analysis
union_check <- read.csv("italy_factor_interaction_union_orfs.csv")
cat("Previous union count:", nrow(union_check), "ORFs\n")
cat("Union calculation matches previous:", length(union_orfs) == nrow(union_check), "\n")

# Save intersection and union lists
intersection_df <- data.frame(ORF = intersection_orfs, stringsAsFactors = FALSE)
union_df <- data.frame(ORF = union_orfs, stringsAsFactors = FALSE)

write.csv(intersection_df, "italy_intersection_orfs_t0-24_AND_t0-30.csv", row.names = FALSE)
write.csv(union_df, "italy_union_orfs_t0-24_OR_t0-30.csv", row.names = FALSE)

# Create detailed breakdown
breakdown_df <- data.frame(
  Category = c("t0-24 only", "t0-30 only", "Both timepoints (intersection)", "Either timepoint (union)"),
  Number_of_ORFs = c(
    length(setdiff(t0_24_orf_list, t0_30_orf_list)),  # t0-24 only
    length(setdiff(t0_30_orf_list, t0_24_orf_list)),  # t0-30 only  
    length(intersection_orfs),                         # intersection
    length(union_orfs)                                 # union
  ),
  stringsAsFactors = FALSE
)

write.csv(breakdown_df, "italy_intersection_union_breakdown.csv", row.names = FALSE)

cat("\nBreakdown:\n")
print(breakdown_df)

# ============================================================================
# PART 3: LOAD DATA AND RECREATE FACTOR MODEL
# ============================================================================

cat("3. Loading data and fitting factor model...\n")

# Load abundance data
abundance_data_raw <- read.csv("Italy.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names
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

# Set factor levels with reversed timeline order
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))
metadata_filtered$onset_timeline_combined <- factor(metadata_filtered$onset_timeline_combined, 
                                                   levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0"))

# Set other factor levels
metadata_filtered$feeding_first_year <- factor(metadata_filtered$feeding_first_year,
                                               levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
metadata_filtered$HLA.Category <- factor(metadata_filtered$HLA.Category,
                                         levels = c("Standard Risk","High Risk","Low/No Risk"))
metadata_filtered$Sex <- factor(metadata_filtered$Sex, levels = c("Female","Male"))
metadata_filtered$Delivery.Mode <- factor(metadata_filtered$Delivery.Mode, levels = c("Vaginal","C-Section"))
metadata_filtered$Age.at.Gluten.Introduction..months. <- as.numeric(metadata_filtered$Age.at.Gluten.Introduction..months.)

# Create design matrix
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

# Convert abundance data to matrix
abundance_matrix <- as.matrix(abundance_filtered)

# Fit model
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get fitted values
fitted_values <- fitted(linear.model.fit)

cat("Factor model fitted successfully\n")

# ============================================================================
# PART 4: CREATE HEATMAP MATRIX FUNCTION
# ============================================================================

create_heatmap_matrix <- function(orf_list, analysis_name) {
  cat("Creating heatmap matrix for", analysis_name, "with", length(orf_list), "ORFs...\n")
  
  all_timepoints <- levels(metadata_filtered$onset_timeline_combined)
  
  heatmap_matrix <- matrix(NA, 
                          nrow = length(orf_list), 
                          ncol = length(all_timepoints),
                          dimnames = list(orf_list, all_timepoints))
  
  for(i in 1:length(orf_list)) {
    orf_id <- orf_list[i]
    
    if(orf_id %in% rownames(fitted_values)) {
      orf_fitted <- fitted_values[orf_id, ]
      
      for(j in 1:length(all_timepoints)) {
        timepoint <- all_timepoints[j]
        
        timepoint_samples <- metadata_filtered$onset_timeline_combined == timepoint
        
        if(sum(timepoint_samples) > 0) {
          timepoint_data <- data.frame(
            fitted_abundance = as.numeric(orf_fitted[timepoint_samples]),
            Dx.Status = metadata_filtered$Dx.Status[timepoint_samples]
          )
          
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
              heatmap_matrix[i, j] <- control_fitted - celiac_fitted
            }
          }
        }
      }
    }
    
    if(i %% 500 == 0 && i > 0) {
      cat("Processed", i, "of", length(orf_list), "ORFs for", analysis_name, "\n")
    }
  }
  
  # Remove ORFs with too many missing values
  valid_orfs <- rowSums(!is.na(heatmap_matrix)) >= (ncol(heatmap_matrix) * 0.5)
  heatmap_matrix_clean <- heatmap_matrix[valid_orfs, ]
  
  cat("Final", analysis_name, "matrix:", nrow(heatmap_matrix_clean), "x", ncol(heatmap_matrix_clean), "\n")
  
  return(heatmap_matrix_clean)
}

# ============================================================================
# PART 5: CREATE INTERSECTION HEATMAP
# ============================================================================

cat("4. Creating intersection heatmap (t0-24 AND t0-30)...\n")

intersection_matrix <- create_heatmap_matrix(intersection_orfs, "intersection")

if(nrow(intersection_matrix) > 0) {
  # Set up color scheme
  data_range <- range(intersection_matrix, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # Create intersection heatmap
  png("italy_intersection_heatmap_t0-24_AND_t0-30.png", 
      width = 2400, height = max(800, nrow(intersection_matrix) * 15), res = 150)
  
  pheatmap(intersection_matrix,
           cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
           color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Italy Cohort: Intersection ORFs Significant in BOTH t0-24 AND t0-30\n",
                       "n =", nrow(intersection_matrix), "ORFs | Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 18, 
           fontsize_row = max(8, min(12, 1000/nrow(intersection_matrix))), 
           fontsize_col = 20,
           cex_main = 1.2,
           show_rownames = nrow(intersection_matrix) <= 50,  # Show row names only if <= 50 ORFs
           show_colnames = TRUE, 
           border_color = NA)
  
  dev.off()
  
  cat("Intersection heatmap created successfully!\n")
  
  # Save intersection matrix data
  write.csv(intersection_matrix, "italy_intersection_heatmap_matrix.csv")
}

# ============================================================================
# PART 6: CREATE UNION HEATMAP  
# ============================================================================

cat("5. Creating union heatmap (t0-24 OR t0-30)...\n")

union_matrix <- create_heatmap_matrix(union_orfs, "union")

if(nrow(union_matrix) > 0) {
  # Set up color scheme
  data_range <- range(union_matrix, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # Create union heatmap
  png("italy_union_heatmap_t0-24_OR_t0-30.png", 
      width = 2400, height = max(800, nrow(union_matrix) * 12), res = 150)
  
  pheatmap(union_matrix,
           cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
           color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("Italy Cohort: Union ORFs Significant in EITHER t0-24 OR t0-30\n",
                       "n =", nrow(union_matrix), "ORFs | Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 18, 
           fontsize_row = max(6, min(10, 800/nrow(union_matrix))), 
           fontsize_col = 20,
           cex_main = 1.2,
           show_rownames = FALSE,  # Too many to show individual names
           show_colnames = TRUE, 
           border_color = NA)
  
  dev.off()
  
  cat("Union heatmap created successfully!\n")
  
  # Save union matrix data
  write.csv(union_matrix, "italy_union_heatmap_matrix.csv")
  
  # Create top variable ORFs version if there are many ORFs
  if(nrow(union_matrix) > 200) {
    orf_variances <- apply(union_matrix, 1, function(x) var(x, na.rm = TRUE))
    orf_variances <- orf_variances[!is.na(orf_variances)]
    n_top <- min(200, length(orf_variances))
    top_variable_orfs <- names(sort(orf_variances, decreasing = TRUE)[1:n_top])
    top_union_matrix <- union_matrix[top_variable_orfs, ]
    
    png("italy_union_heatmap_top200_t0-24_OR_t0-30.png", 
        width = 2400, height = max(800, n_top * 15), res = 150)
    
    pheatmap(top_union_matrix,
             cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
             color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
             main = paste("Italy Cohort: Top", n_top, "Most Variable Union ORFs (t0-24 OR t0-30)\n",
                         "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
             fontsize = 18, 
             fontsize_row = max(8, min(12, 1000/n_top)), 
             fontsize_col = 20,
             cex_main = 1.2,
             show_rownames = TRUE,
             show_colnames = TRUE, 
             border_color = NA)
    
    dev.off()
    
    cat("Top", n_top, "variable union ORFs heatmap also created\n")
  }
}

# ============================================================================
# PART 7: CREATE SUMMARY STATISTICS
# ============================================================================

cat("6. Creating summary statistics...\n")

# Create summaries for both matrices
if(nrow(intersection_matrix) > 0) {
  intersection_summary <- data.frame(
    Timepoint = colnames(intersection_matrix),
    Valid_ORFs = colSums(!is.na(intersection_matrix)),
    Mean_Difference = round(colMeans(intersection_matrix, na.rm = TRUE), 3),
    SD_Difference = round(apply(intersection_matrix, 2, sd, na.rm = TRUE), 3),
    Range_Min = round(apply(intersection_matrix, 2, min, na.rm = TRUE), 3),
    Range_Max = round(apply(intersection_matrix, 2, max, na.rm = TRUE), 3),
    Analysis = "Intersection",
    stringsAsFactors = FALSE
  )
  
  write.csv(intersection_summary, "italy_intersection_heatmap_summary.csv", row.names = FALSE)
}

if(nrow(union_matrix) > 0) {
  union_summary <- data.frame(
    Timepoint = colnames(union_matrix),
    Valid_ORFs = colSums(!is.na(union_matrix)),
    Mean_Difference = round(colMeans(union_matrix, na.rm = TRUE), 3),
    SD_Difference = round(apply(union_matrix, 2, sd, na.rm = TRUE), 3),
    Range_Min = round(apply(union_matrix, 2, min, na.rm = TRUE), 3),
    Range_Max = round(apply(union_matrix, 2, max, na.rm = TRUE), 3),
    Analysis = "Union",
    stringsAsFactors = FALSE
  )
  
  write.csv(union_summary, "italy_union_heatmap_summary.csv", row.names = FALSE)
}

# ============================================================================
# PART 8: FINAL SUMMARY
# ============================================================================

cat("\n=== ITALY INTERSECTION AND UNION HEATMAPS COMPLETED ===\n")
cat("Generated files:\n")

if(nrow(intersection_matrix) > 0) {
  cat("- italy_intersection_heatmap_t0-24_AND_t0-30.png: Intersection heatmap (", nrow(intersection_matrix), "ORFs)\n")
  cat("- italy_intersection_heatmap_matrix.csv: Intersection matrix data\n")
  cat("- italy_intersection_heatmap_summary.csv: Intersection summary statistics\n")
}

if(nrow(union_matrix) > 0) {
  cat("- italy_union_heatmap_t0-24_OR_t0-30.png: Union heatmap (", nrow(union_matrix), "ORFs)\n")
  cat("- italy_union_heatmap_matrix.csv: Union matrix data\n")
  cat("- italy_union_heatmap_summary.csv: Union summary statistics\n")
  
  if(exists("n_top")) {
    cat("- italy_union_heatmap_top200_t0-24_OR_t0-30.png: Top", n_top, "variable union ORFs\n")
  }
}

cat("- italy_intersection_orfs_t0-24_AND_t0-30.csv: List of intersection ORFs\n")
cat("- italy_union_orfs_t0-24_OR_t0-30.csv: List of union ORFs\n")
cat("- italy_intersection_union_breakdown.csv: Detailed breakdown\n")

cat("\nFinal breakdown:\n")
print(breakdown_df)

if(nrow(intersection_matrix) > 0) {
  cat("\nIntersection summary:\n")
  print(intersection_summary)
}

if(nrow(union_matrix) > 0) {
  cat("\nUnion summary:\n")
  print(union_summary)
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")