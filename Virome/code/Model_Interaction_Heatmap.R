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
# extract results
US.results_t0_6 <- topTable(US.contrast_fit, coef = "CONTROL_vs_CELIAC_at_t0_6",number=Inf) %>% filter(adj.P.Val < 0.05)
# 3

# Get ORF lists
t6_orfs <- rownames(US.results_t0_6)


US.results_t0_12 <- topTable(US.contrast_fit, coef = "CONTROL_vs_CELIAC_at_t0_12",number=Inf) %>% filter(adj.P.Val < 0.05)
# 43

t12_orfs <- rownames(US.results_t0_12)


US.results_t0_18 <- topTable(US.contrast_fit, coef = "CONTROL_vs_CELIAC_at_t0_18",number=Inf) %>% filter(adj.P.Val < 0.05)
# 2730

t18_orfs <- rownames(US.results_t0_18)

US.results_t0_24 <- topTable(US.contrast_fit, coef = "CONTROL_vs_CELIAC_at_t0_24",number=Inf) %>% filter(adj.P.Val < 0.05)
# 219

t24_orfs <- rownames(US.results_t0_24)

US.results_t0_30 <- topTable(US.contrast_fit, coef = "CONTROL_vs_CELIAC_at_t0_30",number=Inf) %>% filter(adj.P.Val < 0.05)
dim(US.results_t0_30)
# 101

t30_orfs <- rownames(US.results_t0_30)

US.results_t0_36 <- topTable(US.contrast_fit, coef = "CONTROL_vs_CELIAC_at_t0_36",number=Inf) %>% filter(adj.P.Val < 0.05)
# 4

t36_orfs <- rownames(US.results_t0_36)

US.results_t0_over42 <- topTable(US.contrast_fit, coef = "CONTROL_vs_CELIAC_at_t0_over42",number=Inf) %>% filter(adj.P.Val < 0.05)
# 2641

t42_orfs <- rownames(US.results_t0_over42)


# Create UNION of all significant ORFs (instead of intersection)
union_all <- unique(c(t6_orfs,t12_orfs,t18_orfs,t24_orfs, t30_orfs,t36_orfs, t42_orfs))

# 
# # Also calculate overlaps for comparison and summary
# overlap_all <- intersect(intersect(t24_orfs, t36_orfs), t42_orfs)
# overlap_24_36 <- intersect(t24_orfs, t36_orfs)
# overlap_24_42 <- intersect(t24_orfs, t42_orfs)
# overlap_36_42 <- intersect(t36_orfs, t42_orfs)
# 
# cat("\nORF set analysis:\n")
# cat("- Union of all three timepoints:", length(union_all), "ORFs\n")
# cat("- t0-24 ∩ t0-36:", length(overlap_24_36), "ORFs\n")
# cat("- t0-24 ∩ t0-over42:", length(overlap_24_42), "ORFs\n")
# cat("- t0-36 ∩ t0-over42:", length(overlap_36_42), "ORFs\n")
# cat("- All three timepoints (intersection):", length(overlap_all), "ORFs\n")

# Create breakdown of which ORFs are in which timepoints
union_breakdown <- data.frame(
  ORF = union_all,
  in_t6 = union_all %in% t6_orfs,
  in_t12 = union_all %in% t12_orfs,
  in_t18 = union_all %in% t18_orfs,
  in_t24 = union_all %in% t24_orfs,
  in_t30 = union_all %in% t30_orfs,
  in_t36 = union_all %in% t36_orfs,
  in_t42 = union_all %in% t42_orfs,
  stringsAsFactors = FALSE
)

# # Add summary column showing timepoint membership
# union_breakdown$timepoint_pattern <- apply(union_breakdown[, c("in_t24", "in_t36", "in_t42")], 1, function(x) {
#   timepoints <- c("t24", "t36", "t42")[x]
#   if(length(timepoints) == 3) return("All_three")
#   if(length(timepoints) == 2) return(paste(timepoints, collapse = "_"))
#   if(length(timepoints) == 1) return(paste(timepoints, "only"))
#   return("None")
# })
# 
# # Summary statistics
# timepoint_summary <- table(union_breakdown$timepoint_pattern)
# cat("\nBreakdown of union ORFs by timepoint membership:\n")
# print(timepoint_summary)

# ============================================================================
# PART 2: FIT FACTOR MODEL AND GET FITTED VALUES
# ============================================================================

US.categorical.model.design <- model.matrix(
  ~ Dx.Status * onset_timeline_combined +
    Sex + Age.at.Gluten.Introduction..months. +
    HLA.Category + feeding_first_year + Delivery.Mode,
  data = US.metadata.clean
)
colnames(US.categorical.model.design) <- make.names(colnames(US.categorical.model.design))

# DGE object
US.dge <- DGEList(counts = US.orf.abundance.clean)


# 1) Raw library sizes (per sample; don’t drop any samples here)
US.lib_raw <- colSums(US.dge$counts)

# 2) Choose a sane floor for the "smallest library" used ONLY for the CPM cutoff
#    (10th percentile OR 10% of median, whichever is larger)
US.Lfloor  <- max(quantile(US.lib_raw, 0.10), 0.10 * median(US.lib_raw))

# 3) Feature filtering using the floored library sizes (samples are NOT removed)
US.keep <- filterByExpr(US.dge, design = US.categorical.model.design, min.count = 10,
                        lib.size = pmax(US.lib_raw, US.Lfloor))

# 4) Apply the row filter and recompute library sizes for downstream steps
US.dge.filter  <- US.dge[US.keep, , keep.lib.sizes = FALSE]
# 3342  197

# computes normalization factors—one per sample—that scale library sizes so counts become comparable across samples despite different sequencing depths and composition
US.dge.filter <- calcNormFactors(US.dge.filter, method = "TMMwsp")

# 1. voom step (per-observation precision weights): Converts counts to log2-CPM, learns the mean–variance trend, and assigns a weight to each data point 
# 2. sample-quality step (array-style weights): Estimates one extra weight per sample that reflects how globally noisy/odd that sample is. Noisy samples get weights < 1 (down-weighted); very clean ones can be > 1 (up-weighted).
US.v <- voomWithQualityWeights(US.dge.filter, US.categorical.model.design, plot = TRUE)


# Estimate within-patient correlation (repeated measures)
US.corfit <- duplicateCorrelation(US.v, US.categorical.model.design,
                                  block = US.metadata.clean$patientID)
cat("Consensus correlation =", round(US.corfit$consensus, 3), "\n")

# Fit the weighted linear model + robust EB shrinkage
US.categorical.model.fit <- lmFit(US.v, US.categorical.model.design,
                                  block = US.metadata.clean$patientID,
                                  correlation = US.corfit$consensus)
US.categorical.model.fit <- eBayes(US.categorical.model.fit, trend = TRUE, robust = TRUE)

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