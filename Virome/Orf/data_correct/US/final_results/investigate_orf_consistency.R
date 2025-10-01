#!/usr/bin/env Rscript
# Investigate ORF Consistency Pattern Analysis
# Author: Claude AI
# Date: 2025-08-12
# Purpose: Investigate why ORFs show such consistent temporal patterns

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)
library(corrplot)

cat("=== INVESTIGATING ORF CONSISTENCY PATTERNS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND FIT MODELS
# ============================================================================

cat("1. Loading data and fitting factor model...\n")

# Load abundance data
abundance_data_raw <- read.csv("US.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names
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

# Set reference levels with reversed timeline order
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))
metadata_filtered$onset_timeline_combined <- factor(metadata_filtered$onset_timeline_combined, 
                                                   levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0"))

# Create design matrix
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

# Convert abundance data to matrix
abundance_matrix <- as.matrix(abundance_filtered)

# Fit factor model
cat("Fitting factor model with patient blocking...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get fitted values
fitted_values <- fitted(linear.model.fit)

cat("Factor model fitted successfully\n")

# ============================================================================
# PART 2: LOAD SIGNIFICANT AND NON-SIGNIFICANT ORFS
# ============================================================================

cat("2. Identifying significant vs non-significant ORFs...\n")

# Load union significant ORFs from previous analysis
union_orfs <- read.csv("factor_interaction_union_orfs.csv")
significant_orfs <- union_orfs$ORF

cat("Loaded", length(significant_orfs), "interaction-significant ORFs\n")

# Get all ORFs in the dataset
all_orfs <- rownames(abundance_matrix)
cat("Total ORFs in dataset:", length(all_orfs), "\n")

# Identify non-significant ORFs
non_significant_orfs <- setdiff(all_orfs, significant_orfs)
cat("Non-significant ORFs:", length(non_significant_orfs), "\n")

# ============================================================================
# PART 3: CREATE HEATMAP FUNCTION
# ============================================================================

create_heatmap_matrix <- function(orf_list, name_suffix) {
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
    
    if(i %% 500 == 0) {
      cat("Processed", i, "of", length(orf_list), "ORFs for", name_suffix, "\n")
    }
  }
  
  # Remove ORFs with too many missing values
  valid_orfs <- rowSums(!is.na(heatmap_matrix)) >= (ncol(heatmap_matrix) * 0.5)
  heatmap_matrix_clean <- heatmap_matrix[valid_orfs, ]
  
  return(heatmap_matrix_clean)
}

# ============================================================================
# PART 4: POINT 1 - CHECK NON-SIGNIFICANT ORFS PATTERN
# ============================================================================

cat("3. Point 1: Analyzing non-significant ORFs pattern...\n")

# Sample random subset of non-significant ORFs for comparison (to avoid memory issues)
set.seed(123)
sample_size <- min(1000, length(non_significant_orfs))
sampled_non_sig_orfs <- sample(non_significant_orfs, sample_size)

cat("Sampling", sample_size, "non-significant ORFs for comparison\n")

# Create heatmap matrices
significant_matrix <- create_heatmap_matrix(significant_orfs[1:min(1000, length(significant_orfs))], "significant")
non_significant_matrix <- create_heatmap_matrix(sampled_non_sig_orfs, "non-significant")

cat("Significant ORFs matrix:", nrow(significant_matrix), "x", ncol(significant_matrix), "\n")
cat("Non-significant ORFs matrix:", nrow(non_significant_matrix), "x", ncol(non_significant_matrix), "\n")

# Calculate consistency metrics
calculate_consistency <- function(matrix_data, name) {
  # Calculate correlation between timepoints
  timepoint_correlations <- cor(matrix_data, use = "complete.obs")
  mean_correlation <- mean(timepoint_correlations[upper.tri(timepoint_correlations)], na.rm = TRUE)
  
  # Calculate proportion of ORFs following same direction at each timepoint
  timepoint_directions <- apply(matrix_data, 2, function(x) {
    pos_prop <- sum(x > 0, na.rm = TRUE) / sum(!is.na(x))
    neg_prop <- sum(x < 0, na.rm = TRUE) / sum(!is.na(x))
    return(c(positive = pos_prop, negative = neg_prop))
  })
  
  # Overall consistency score (how much majority direction dominates)
  consistency_scores <- apply(timepoint_directions, 2, max)
  overall_consistency <- mean(consistency_scores)
  
  cat("=== Consistency Analysis for", name, "ORFs ===\n")
  cat("Mean inter-timepoint correlation:", round(mean_correlation, 3), "\n")
  cat("Overall consistency score:", round(overall_consistency, 3), "\n")
  cat("Data range:", round(range(matrix_data, na.rm = TRUE), 2), "\n")
  
  return(list(
    mean_correlation = mean_correlation,
    consistency_score = overall_consistency,
    timepoint_directions = timepoint_directions,
    data_range = range(matrix_data, na.rm = TRUE)
  ))
}

significant_consistency <- calculate_consistency(significant_matrix, "Significant")
non_significant_consistency <- calculate_consistency(non_significant_matrix, "Non-significant")

# Create comparison heatmaps
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Significant ORFs heatmap (sample)
png("consistency_significant_sample_heatmap.png", width = 2000, height = 1500, res = 150)
data_range_sig <- range(significant_matrix, na.rm = TRUE)
max_abs_sig <- max(abs(data_range_sig))
pheatmap(significant_matrix,
         cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
         color = colors, breaks = seq(-max_abs_sig, max_abs_sig, length.out = 101),
         main = paste("Sample of Interaction-Significant ORFs (n =", nrow(significant_matrix), ")\nRed = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 14, show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
dev.off()

# Non-significant ORFs heatmap
png("consistency_non_significant_sample_heatmap.png", width = 2000, height = 1500, res = 150)
data_range_nonsig <- range(non_significant_matrix, na.rm = TRUE)
max_abs_nonsig <- max(abs(data_range_nonsig))
pheatmap(non_significant_matrix,
         cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
         color = colors, breaks = seq(-max_abs_nonsig, max_abs_nonsig, length.out = 101),
         main = paste("Sample of Non-Interaction-Significant ORFs (n =", nrow(non_significant_matrix), ")\nRed = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 14, show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
dev.off()

# ============================================================================
# PART 5: POINT 2 - EXAMINE RAW ABUNDANCE DISTRIBUTIONS
# ============================================================================

cat("4. Point 2: Examining raw abundance distributions...\n")

# Sample ORFs for raw abundance analysis
sample_orfs_for_raw <- significant_orfs[1:min(50, length(significant_orfs))]

# Extract raw abundances for sampled ORFs
raw_abundance_data <- abundance_filtered[sample_orfs_for_raw, ]

# Create long format for plotting
raw_data_long <- data.frame()
for(orf in rownames(raw_abundance_data)) {
  for(sample in colnames(raw_abundance_data)) {
    if(sample %in% rownames(metadata_filtered)) {
      raw_data_long <- rbind(raw_data_long, data.frame(
        ORF = orf,
        Sample = sample,
        Raw_Abundance = as.numeric(raw_abundance_data[orf, sample]),
        Dx_Status = metadata_filtered[sample, "Dx.Status"],
        Timeline = metadata_filtered[sample, "onset_timeline_combined"],
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Calculate raw abundance summary statistics
raw_abundance_summary <- raw_data_long %>%
  group_by(ORF, Timeline, Dx_Status) %>%
  summarise(
    mean_abundance = mean(Raw_Abundance, na.rm = TRUE),
    median_abundance = median(Raw_Abundance, na.rm = TRUE),
    sd_abundance = sd(Raw_Abundance, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  filter(n_samples >= 2)

# Create raw abundance comparison plot
cat("Creating raw abundance distribution plots...\n")

# Box plot of raw abundances by timeline and disease status
png("raw_abundance_distribution_by_timeline.png", width = 2400, height = 1200, res = 150)
p1 <- ggplot(raw_data_long, aes(x = Timeline, y = log10(Raw_Abundance + 1), fill = Dx_Status)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("CELIAC" = "#2166AC", "CONTROL" = "#B2182B")) +
  labs(title = "Raw Abundance Distributions by Timeline and Disease Status",
       subtitle = paste("Sample of", length(sample_orfs_for_raw), "interaction-significant ORFs"),
       x = "Timeline", y = "Log10(Raw Abundance + 1)", fill = "Disease Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))
print(p1)
dev.off()

# Summary statistics
raw_stats_summary <- raw_data_long %>%
  group_by(Timeline, Dx_Status) %>%
  summarise(
    mean_log_abundance = mean(log10(Raw_Abundance + 1), na.rm = TRUE),
    median_log_abundance = median(log10(Raw_Abundance + 1), na.rm = TRUE),
    sd_log_abundance = sd(log10(Raw_Abundance + 1), na.rm = TRUE),
    n_observations = n(),
    .groups = "drop"
  )

write.csv(raw_stats_summary, "raw_abundance_summary_stats.csv", row.names = FALSE)

cat("Raw abundance analysis completed\n")

# ============================================================================
# PART 6: POINT 3 - INDIVIDUAL ORF VARIANCE AND OUTLIERS
# ============================================================================

cat("5. Point 3: Analyzing individual ORF variance and identifying outliers...\n")

# Calculate variance metrics for each ORF
orf_variance_analysis <- data.frame()

for(orf in rownames(significant_matrix)) {
  orf_values <- significant_matrix[orf, ]
  orf_values_clean <- orf_values[!is.na(orf_values)]
  
  if(length(orf_values_clean) >= 4) {  # Need at least 4 timepoints
    variance <- var(orf_values_clean)
    range_val <- diff(range(orf_values_clean))
    mean_abs_value <- mean(abs(orf_values_clean))
    
    # Check if ORF follows general pattern (correlation with mean pattern)
    mean_pattern <- colMeans(significant_matrix, na.rm = TRUE)
    mean_pattern_clean <- mean_pattern[!is.na(orf_values)]
    orf_values_for_cor <- orf_values_clean[names(orf_values_clean) %in% names(mean_pattern_clean)]
    mean_pattern_for_cor <- mean_pattern_clean[names(mean_pattern_clean) %in% names(orf_values_for_cor)]
    
    if(length(orf_values_for_cor) >= 3 && length(mean_pattern_for_cor) >= 3) {
      pattern_correlation <- cor(orf_values_for_cor, mean_pattern_for_cor, use = "complete.obs")
    } else {
      pattern_correlation <- NA
    }
    
    orf_variance_analysis <- rbind(orf_variance_analysis, data.frame(
      ORF = orf,
      Variance = variance,
      Range = range_val,
      Mean_Abs_Value = mean_abs_value,
      Pattern_Correlation = pattern_correlation,
      N_Timepoints = length(orf_values_clean),
      stringsAsFactors = FALSE
    ))
  }
}

# Identify outliers
orf_variance_analysis$Variance_Percentile <- rank(orf_variance_analysis$Variance) / nrow(orf_variance_analysis) * 100
orf_variance_analysis$Correlation_Percentile <- rank(-orf_variance_analysis$Pattern_Correlation, na.last = "keep") / sum(!is.na(orf_variance_analysis$Pattern_Correlation)) * 100

# High variance outliers (top 5%)
high_variance_outliers <- orf_variance_analysis[orf_variance_analysis$Variance_Percentile >= 95, ]
# Low correlation outliers (bottom 5%)
low_correlation_outliers <- orf_variance_analysis[orf_variance_analysis$Correlation_Percentile <= 5 & !is.na(orf_variance_analysis$Correlation_Percentile), ]

cat("=== ORF Variance Analysis ===\n")
cat("Total ORFs analyzed:", nrow(orf_variance_analysis), "\n")
cat("High variance outliers (top 5%):", nrow(high_variance_outliers), "\n")
cat("Low correlation outliers (bottom 5%):", nrow(low_correlation_outliers), "\n")
cat("Mean pattern correlation (all ORFs):", round(mean(orf_variance_analysis$Pattern_Correlation, na.rm = TRUE), 3), "\n")

# Create variance analysis plots
png("orf_variance_analysis.png", width = 2400, height = 1600, res = 150)
par(mfrow = c(2, 2))

# Variance distribution
hist(orf_variance_analysis$Variance, breaks = 50, main = "Distribution of ORF Variances", 
     xlab = "Variance", col = "lightblue", border = "darkblue")
abline(v = quantile(orf_variance_analysis$Variance, 0.95), col = "red", lwd = 2, lty = 2)

# Pattern correlation distribution
hist(orf_variance_analysis$Pattern_Correlation, breaks = 50, main = "Distribution of Pattern Correlations", 
     xlab = "Correlation with Mean Pattern", col = "lightgreen", border = "darkgreen")
abline(v = quantile(orf_variance_analysis$Pattern_Correlation, 0.05, na.rm = TRUE), col = "red", lwd = 2, lty = 2)

# Variance vs Correlation scatter
plot(orf_variance_analysis$Pattern_Correlation, orf_variance_analysis$Variance,
     xlab = "Pattern Correlation", ylab = "Variance", 
     main = "Variance vs Pattern Correlation", pch = 16, col = alpha("black", 0.6))
abline(h = quantile(orf_variance_analysis$Variance, 0.95), col = "red", lwd = 2, lty = 2)
abline(v = quantile(orf_variance_analysis$Pattern_Correlation, 0.05, na.rm = TRUE), col = "red", lwd = 2, lty = 2)

# Range distribution
hist(orf_variance_analysis$Range, breaks = 50, main = "Distribution of ORF Ranges", 
     xlab = "Range (Max - Min)", col = "lightyellow", border = "orange")

dev.off()

# Save detailed results
write.csv(orf_variance_analysis, "orf_variance_detailed_analysis.csv", row.names = FALSE)
write.csv(high_variance_outliers, "high_variance_outlier_orfs.csv", row.names = FALSE)
write.csv(low_correlation_outliers, "low_correlation_outlier_orfs.csv", row.names = FALSE)

# Create heatmap of outlier ORFs
if(nrow(high_variance_outliers) > 0) {
  outlier_orfs <- c(high_variance_outliers$ORF, low_correlation_outliers$ORF)
  outlier_orfs <- unique(outlier_orfs)
  outlier_orfs <- outlier_orfs[outlier_orfs %in% rownames(significant_matrix)]
  
  if(length(outlier_orfs) > 0) {
    outlier_matrix <- significant_matrix[outlier_orfs, , drop = FALSE]
    
    png("outlier_orfs_heatmap.png", width = 2000, height = max(600, length(outlier_orfs) * 20), res = 150)
    data_range_outlier <- range(outlier_matrix, na.rm = TRUE)
    max_abs_outlier <- max(abs(data_range_outlier))
    pheatmap(outlier_matrix,
             cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
             color = colors, breaks = seq(-max_abs_outlier, max_abs_outlier, length.out = 101),
             main = paste("Outlier ORFs: High Variance & Low Pattern Correlation (n =", nrow(outlier_matrix), ")"),
             fontsize = 12, show_rownames = TRUE, show_colnames = TRUE, border_color = NA,
             fontsize_row = max(6, min(10, 400/nrow(outlier_matrix))))
    dev.off()
    
    cat("Created outlier ORFs heatmap with", nrow(outlier_matrix), "ORFs\n")
  }
}

# ============================================================================
# PART 7: SUMMARY AND CONCLUSIONS
# ============================================================================

cat("6. Creating comprehensive summary...\n")

# Create summary report
summary_report <- data.frame(
  Metric = c("Total ORFs in dataset", "Interaction-significant ORFs", "Non-significant ORFs sampled",
             "Significant ORFs mean correlation", "Non-significant ORFs mean correlation",
             "Significant ORFs consistency score", "Non-significant ORFs consistency score",
             "High variance outliers", "Low correlation outliers",
             "Mean pattern correlation (all significant ORFs)"),
  Value = c(length(all_orfs), length(significant_orfs), sample_size,
            round(significant_consistency$mean_correlation, 3),
            round(non_significant_consistency$mean_correlation, 3),
            round(significant_consistency$consistency_score, 3),
            round(non_significant_consistency$consistency_score, 3),
            nrow(high_variance_outliers), nrow(low_correlation_outliers),
            round(mean(orf_variance_analysis$Pattern_Correlation, na.rm = TRUE), 3)),
  stringsAsFactors = FALSE
)

write.csv(summary_report, "consistency_investigation_summary.csv", row.names = FALSE)

cat("\n=== CONSISTENCY INVESTIGATION COMPLETED ===\n")
cat("Generated files:\n")
cat("- consistency_significant_sample_heatmap.png: Sample of significant ORFs\n")
cat("- consistency_non_significant_sample_heatmap.png: Sample of non-significant ORFs\n")
cat("- raw_abundance_distribution_by_timeline.png: Raw abundance distributions\n")
cat("- orf_variance_analysis.png: Variance analysis plots\n")
cat("- outlier_orfs_heatmap.png: Heatmap of outlier ORFs\n")
cat("- consistency_investigation_summary.csv: Overall summary\n")
cat("- orf_variance_detailed_analysis.csv: Detailed variance metrics\n")
cat("- raw_abundance_summary_stats.csv: Raw abundance statistics\n")

print(summary_report)

cat("\nAnalysis completed at:", format(Sys.time()), "\n")