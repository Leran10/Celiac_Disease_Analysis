#!/usr/bin/env Rscript
# Factor Model Interaction Heatmap - Union of All Significant ORFs - US Cohort
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Create heatmap using UNION of all significant ORFs from Dx.Status * factor(timeline) interactions

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)

cat("=== FACTOR MODEL INTERACTION UNION HEATMAP ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD EXISTING INTERACTION RESULTS AND CREATE UNION
# ============================================================================

cat("1. Loading existing interaction results and creating union...\n")

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

# Create UNION of all significant ORFs (instead of intersection)
union_all <- unique(c(t24_orfs, t36_orfs, t42_orfs))

# Also calculate overlaps for comparison and summary
overlap_all <- intersect(intersect(t24_orfs, t36_orfs), t42_orfs)
overlap_24_36 <- intersect(t24_orfs, t36_orfs)
overlap_24_42 <- intersect(t24_orfs, t42_orfs)
overlap_36_42 <- intersect(t36_orfs, t42_orfs)

cat("\nORF set analysis:\n")
cat("- Union of all three timepoints:", length(union_all), "ORFs\n")
cat("- t0-24 ∩ t0-36:", length(overlap_24_36), "ORFs\n")
cat("- t0-24 ∩ t0-over42:", length(overlap_24_42), "ORFs\n")
cat("- t0-36 ∩ t0-over42:", length(overlap_36_42), "ORFs\n")
cat("- All three timepoints (intersection):", length(overlap_all), "ORFs\n")

# Create breakdown of which ORFs are in which timepoints
union_breakdown <- data.frame(
  ORF = union_all,
  in_t24 = union_all %in% t24_orfs,
  in_t36 = union_all %in% t36_orfs,
  in_t42 = union_all %in% t42_orfs,
  stringsAsFactors = FALSE
)

# Add summary column showing timepoint membership
union_breakdown$timepoint_pattern <- apply(union_breakdown[, c("in_t24", "in_t36", "in_t42")], 1, function(x) {
  timepoints <- c("t24", "t36", "t42")[x]
  if(length(timepoints) == 3) return("All_three")
  if(length(timepoints) == 2) return(paste(timepoints, collapse = "_"))
  if(length(timepoints) == 1) return(paste(timepoints, "only"))
  return("None")
})

# Summary statistics
timepoint_summary <- table(union_breakdown$timepoint_pattern)
cat("\nBreakdown of union ORFs by timepoint membership:\n")
print(timepoint_summary)

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
# PART 3: CREATE HEATMAP MATRIX FOR UNION ORF SET
# ============================================================================

cat("\n3. Creating heatmap matrix for union ORF set...\n")

# Function to calculate temporal differences using fitted values
calculate_factor_temporal_differences <- function(fitted_values, metadata, orf_list) {
  
  # Get unique timepoints in REVERSED order (t0-over42 to t0)
  timepoints <- rev(levels(metadata$onset_timeline_combined))
  
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

# Create heatmap matrix for union ORF set
cat("Creating heatmap for union of all significant ORFs...\n")
heatmap_union <- calculate_factor_temporal_differences(fitted_values, metadata_filtered, union_all)
cat("Created matrix for", length(union_all), "union ORFs\n")

# ============================================================================
# PART 4: CREATE HEATMAP
# ============================================================================

cat("\n4. Creating factor model interaction union heatmap...\n")

# Set up color scheme and create heatmap
create_union_heatmap <- function(matrix_data, title_text, filename) {
  if(nrow(matrix_data) > 0) {
    # Remove rows with all NA
    valid_rows <- rowSums(!is.na(matrix_data)) > 0
    clean_matrix <- matrix_data[valid_rows, ]
    
    if(nrow(clean_matrix) > 0) {
      data_range <- range(clean_matrix, na.rm = TRUE)
      max_abs <- max(abs(data_range), na.rm = TRUE)
      color_limits <- c(-max_abs, max_abs)
      colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
      
      # For large ORF sets, create a larger heatmap
      png(filename, width = 2400, height = max(1200, nrow(clean_matrix) * 12), res = 150)
      pheatmap(clean_matrix,
               cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
               color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
               main = paste(title_text, "\nRed = Higher in CONTROL, Blue = Higher in CELIAC\n(Reversed Timeline Order: t0-over42 → t0)"),
               fontsize = 16, 
               fontsize_row = max(6, min(10, 800/nrow(clean_matrix))), 
               fontsize_col = 18,
               cex_main = 1.3,
               show_rownames = ifelse(nrow(clean_matrix) <= 100, TRUE, FALSE),
               show_colnames = TRUE, 
               border_color = NA)
      dev.off()
      
      cat("Created heatmap:", filename, "with", nrow(clean_matrix), "ORFs\n")
      return(clean_matrix)
    }
  }
  return(NULL)
}

# Create the union heatmap
final_union_matrix <- create_union_heatmap(heatmap_union, 
  paste("US Cohort: Factor Model - Union of", length(union_all), "ORFs Significant in ANY Interaction Timepoint"),
  "factor_interaction_union_all_heatmap.png")

# ============================================================================
# PART 5: SAVE RESULTS AND SUMMARY
# ============================================================================

cat("\n5. Saving results and creating summary...\n")

# Save union ORF breakdown
write.csv(union_breakdown, "factor_interaction_union_breakdown.csv", row.names = FALSE)

# Save union summary statistics
union_summary <- data.frame(
  Set = c("Union (any timepoint)", "t0-24 only", "t0-36 only", "t0-over42 only", 
          "t0-24 ∩ t0-36", "t0-24 ∩ t0-over42", "t0-36 ∩ t0-over42",
          "All three timepoints (intersection)"),
  Number_of_ORFs = c(length(union_all),
                      length(setdiff(setdiff(t24_orfs, t36_orfs), t42_orfs)),
                      length(setdiff(setdiff(t36_orfs, t24_orfs), t42_orfs)),
                      length(setdiff(setdiff(t42_orfs, t24_orfs), t36_orfs)),
                      length(overlap_24_36),
                      length(overlap_24_42), 
                      length(overlap_36_42),
                      length(overlap_all)),
  Percentage_of_Union = c(100.0,
                          round(100 * length(setdiff(setdiff(t24_orfs, t36_orfs), t42_orfs)) / length(union_all), 1),
                          round(100 * length(setdiff(setdiff(t36_orfs, t24_orfs), t42_orfs)) / length(union_all), 1),
                          round(100 * length(setdiff(setdiff(t42_orfs, t24_orfs), t36_orfs)) / length(union_all), 1),
                          round(100 * length(overlap_24_36) / length(union_all), 1),
                          round(100 * length(overlap_24_42) / length(union_all), 1),
                          round(100 * length(overlap_36_42) / length(union_all), 1),
                          round(100 * length(overlap_all) / length(union_all), 1)),
  stringsAsFactors = FALSE
)

write.csv(union_summary, "factor_interaction_union_summary.csv", row.names = FALSE)

# Save union ORF list
write.csv(data.frame(ORF = union_all), "factor_interaction_union_orfs.csv", row.names = FALSE)

# Save heatmap matrix
if(!is.null(final_union_matrix)) {
  write.csv(final_union_matrix, "factor_interaction_union_matrix.csv")
}

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== FACTOR MODEL INTERACTION UNION HEATMAP ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- factor_interaction_union_all_heatmap.png: Heatmap of", length(union_all), "ORFs significant in ANY interaction timepoint\n")
cat("- factor_interaction_union_orfs.csv: Complete list of union ORFs\n")
cat("- factor_interaction_union_matrix.csv: Heatmap matrix data\n")
cat("- factor_interaction_union_breakdown.csv: Detailed breakdown of which ORFs appear in which timepoints\n")
cat("- factor_interaction_union_summary.csv: Summary statistics\n")

cat("\nUnion Summary:\n")
print(union_summary)

cat("\nTimepoint Membership Breakdown:\n")
print(timepoint_summary)

cat("\nKey Insights:\n")
cat("- Union includes", length(union_all), "ORFs that are significant in at least one interaction timepoint\n")
cat("- This is", round(length(union_all) / length(overlap_all), 1), "times larger than the intersection (", length(overlap_all), "ORFs)\n")

# Calculate coverage statistics
cat("- Coverage by timepoint:\n")
cat("  * t0-24:", sum(union_breakdown$in_t24), "ORFs (", round(100*sum(union_breakdown$in_t24)/length(union_all), 1), "%)\n")
cat("  * t0-36:", sum(union_breakdown$in_t36), "ORFs (", round(100*sum(union_breakdown$in_t36)/length(union_all), 1), "%)\n")
cat("  * t0-over42:", sum(union_breakdown$in_t42), "ORFs (", round(100*sum(union_breakdown$in_t42)/length(union_all), 1), "%)\n")

cat("- Unique to single timepoints:", sum(timepoint_summary[grepl("only", names(timepoint_summary))]), "ORFs\n")
cat("- Present in multiple timepoints:", sum(timepoint_summary[!grepl("only", names(timepoint_summary)) & names(timepoint_summary) != "All_three"]), "ORFs\n")
cat("- Present in all three timepoints:", timepoint_summary["All_three"], "ORFs\n")

cat("\nThis union approach captures the full breadth of temporal interaction effects,\n")
cat("including ORFs that may be specific to individual timepoints as well as those with consistent patterns.\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")