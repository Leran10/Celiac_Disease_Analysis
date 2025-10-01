#!/usr/bin/env Rscript
# Comprehensive US Cohort Analysis - Matching Total Cohort Methodology
# Author: Claude AI
# Date: 2025-08-05

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)
library(vegan)
library(gridExtra)
library(ggrepel)
library(stringr)

cat("=== COMPREHENSIVE US COHORT ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND FIT LIMMA MODELS
# ============================================================================

cat("1. Loading US cohort data and fitting limma models...\n")

# Load abundance data
abundance_data_raw <- read.csv("../US.orf.abundance.clean.csv")
rownames(abundance_data_raw) <- abundance_data_raw$X
abundance_data <- abundance_data_raw[, -1]

# Fix sample names (remove X prefix if present)
sample_names <- colnames(abundance_data)
sample_names_clean <- gsub("^X", "", sample_names)
colnames(abundance_data) <- sample_names_clean

# Load metadata
metadata_raw <- read.csv("../US.metadata.clean.csv")
rownames(metadata_raw) <- metadata_raw$X
metadata <- metadata_raw

cat("Loaded abundance data with", nrow(abundance_data), "ORFs and", ncol(abundance_data), "samples\n")
cat("Loaded metadata with", nrow(metadata), "samples\n")

# Match samples
matching_samples <- intersect(colnames(abundance_data), rownames(metadata))
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")

# Ensure Dx.Status is a factor with CELIAC as reference
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))

# Create design matrix (US cohort doesn't have Country)
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_numeric + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

cat("Design matrix dimensions:", nrow(linear.model.design), "x", ncol(linear.model.design), "\n")

# Convert abundance data to matrix for limma
abundance_matrix <- as.matrix(abundance_filtered)

# Fit limma model with voomLmFit and patient blocking
cat("Running voomLmFit with patient blocking...\n")
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

# Get significant ORFs for the INTERACTION term (Dx.StatusCONTROL:onset_timeline_numeric)
topTable_results <- topTable(linear.model.fit, coef = "Dx.StatusCONTROL:onset_timeline_numeric", number = Inf, sort.by = "P")
significant_orfs <- topTable_results[topTable_results$adj.P.Val < 0.05, ]

cat("Identified", nrow(significant_orfs), "significant ORFs at adj.p < 0.05 (interaction term)\n")
cat("Percentage of total ORFs:", round(nrow(significant_orfs)/nrow(abundance_matrix)*100, 1), "%\n")

# Save significant ORFs results
write.csv(significant_orfs, "US_limma_model_Sig_res_final.csv", row.names = TRUE)

# ============================================================================
# PART 2: CREATE TEMPORAL HEATMAPS WITH LIMMA FITTED VALUES
# ============================================================================

cat("\n2. Creating temporal heatmaps with limma fitted values...\n")

# Get fitted values from the model
fitted_values <- fitted(linear.model.fit)

# Create model_results list for significant ORFs
significant_orf_names <- rownames(significant_orfs)
available_orfs <- rownames(fitted_values)
matching_sig_orfs <- intersect(significant_orf_names, available_orfs)

model_results <- list()
for(orf_id in matching_sig_orfs) {
  orf_fitted <- fitted_values[orf_id, ]
  model_results[[orf_id]] <- data.frame(
    onset_timeline_numeric = metadata_filtered$onset_timeline_numeric,
    Dx.Status = metadata_filtered$Dx.Status,
    fitted_abundance = as.numeric(orf_fitted),
    stringsAsFactors = FALSE
  )
}

cat("Created model_results list with", length(model_results), "ORFs\n")

# Define the calculate_temporal_differences function
calculate_temporal_differences <- function(model_results, time_bins) {
  heatmap_matrix <- matrix(NA, 
                          nrow = length(model_results), 
                          ncol = length(time_bins),
                          dimnames = list(names(model_results), 
                                         ifelse(time_bins == 0, "T0", paste0("T", time_bins))))
  
  for(i in 1:length(model_results)) {
    orf_id <- names(model_results)[i]
    orf_data <- model_results[[orf_id]]
    
    orf_data$time_bin <- NA
    for(timepoint in time_bins) {
      matches <- abs(orf_data$onset_timeline_numeric - timepoint) < 1
      orf_data$time_bin[matches] <- timepoint
    }
    
    for(j in 1:length(time_bins)) {
      timepoint <- time_bins[j]
      time_data <- orf_data[!is.na(orf_data$time_bin) & orf_data$time_bin == timepoint, ]
      
      if(nrow(time_data) >= 2) {
        group_means <- time_data %>%
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
    
    if(i %% 100 == 0) {
      cat("Processed", i, "of", length(model_results), "ORFs for temporal differences\n")
    }
  }
  
  return(heatmap_matrix)
}

# Calculate temporal differences
time_bins <- sort(unique(metadata_filtered$onset_timeline_numeric))
heatmap_matrix <- calculate_temporal_differences(model_results, time_bins)

# Filter out ORFs with too many missing values
valid_orfs <- rowSums(!is.na(heatmap_matrix)) >= (ncol(heatmap_matrix) * 0.3)
heatmap_matrix_clean <- heatmap_matrix[valid_orfs, ]

cat("Final heatmap matrix:", nrow(heatmap_matrix_clean), "ORFs x", ncol(heatmap_matrix_clean), "timepoints\n")

# Create heatmaps
if(nrow(heatmap_matrix_clean) > 0) {
  data_range <- range(heatmap_matrix_clean, na.rm = TRUE)
  max_abs <- max(abs(data_range), na.rm = TRUE)
  color_limits <- c(-max_abs, max_abs)
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # All significant ORFs heatmap
  png("all_significant_orfs_temporal_heatmap_US_limma_fitted.png", width = 1400, height = 1000, res = 150)
  pheatmap(heatmap_matrix_clean,
           cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
           color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("US Cohort: Temporal Heatmap -", nrow(heatmap_matrix_clean), "Significant ORFs (Limma Fitted Values)\n",
                        "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 18, fontsize_row = 12, fontsize_col = 20,
           cex_main = 1.5,
           show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
  dev.off()
  
  # Top variable ORFs heatmap
  orf_variances <- apply(heatmap_matrix_clean, 1, function(x) var(x, na.rm = TRUE))
  orf_variances <- orf_variances[!is.na(orf_variances)]
  n_top <- min(200, length(orf_variances))
  top_variable_orfs <- names(sort(orf_variances, decreasing = TRUE)[1:n_top])
  top_heatmap_matrix <- heatmap_matrix_clean[top_variable_orfs, ]
  
  png("top_variable_significant_orfs_heatmap_US_limma_fitted.png", 
      width = 1800, height = max(800, n_top * 15), res = 150)
  pheatmap(top_heatmap_matrix,
           cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
           color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
           main = paste("US Cohort: Top", n_top, "Most Variable Significant ORFs (Limma Fitted Values)\n",
                        "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
           fontsize = 18, fontsize_row = max(10, min(14, 1200/n_top)), fontsize_col = 20,
           cex_main = 1.5,
           show_rownames = TRUE, show_colnames = TRUE, border_color = NA)
  dev.off()
}

cat("Temporal heatmaps created\n")

# Save data
write.csv(heatmap_matrix_clean, "temporal_heatmap_matrix_US_limma_fitted.csv")
saveRDS(model_results, "limma_model_results_fitted_values_US.rds")

# ============================================================================
# PART 3: PATTERN CLUSTERING ANALYSIS
# ============================================================================

cat("\n3. Creating pattern clustering analysis...\n")

# Function to calculate pattern features
calculate_pattern_features <- function(orf_profile) {
  valid_values <- orf_profile[!is.na(orf_profile)]
  if(length(valid_values) < 3) {
    return(list(slope = NA, variance = NA, early_late_diff = NA, 
                peak_timing = NA, monotonicity = NA, pattern_type = "insufficient_data"))
  }
  
  timepoints <- 1:length(valid_values)
  slope <- coef(lm(valid_values ~ timepoints))[2]
  variance <- var(valid_values)
  
  n <- length(valid_values)
  early_mean <- mean(valid_values[1:min(3, n)])
  late_mean <- mean(valid_values[max(1, n-2):n])
  early_late_diff <- late_mean - early_mean
  
  peak_timing <- which.max(abs(valid_values)) / length(valid_values)
  diffs <- diff(valid_values)
  monotonicity <- abs(sum(sign(diffs))) / length(diffs)
  
  if(abs(slope) < 0.1 && variance < 1) {
    pattern_type <- "stable"
  } else if(slope > 0.1) {
    pattern_type <- "increasing"
  } else if(slope < -0.1) {
    pattern_type <- "decreasing"
  } else {
    pattern_type <- "complex"
  }
  
  return(list(slope = slope, variance = variance, early_late_diff = early_late_diff,
              peak_timing = peak_timing, monotonicity = monotonicity, pattern_type = pattern_type))
}

# Calculate features for all ORFs
pattern_features <- list()
for(i in 1:nrow(heatmap_matrix_clean)) {
  orf_name <- rownames(heatmap_matrix_clean)[i]
  orf_profile <- as.numeric(heatmap_matrix_clean[i, ])
  pattern_features[[orf_name]] <- calculate_pattern_features(orf_profile)
}

# Convert to data frame
features_df <- data.frame(
  orf_name = names(pattern_features),
  slope = sapply(pattern_features, function(x) x$slope),
  variance = sapply(pattern_features, function(x) x$variance),
  early_late_diff = sapply(pattern_features, function(x) x$early_late_diff),
  peak_timing = sapply(pattern_features, function(x) x$peak_timing),
  monotonicity = sapply(pattern_features, function(x) x$monotonicity),
  pattern_type = sapply(pattern_features, function(x) x$pattern_type),
  stringsAsFactors = FALSE
)

features_clean <- features_df[!is.na(features_df$slope), ]

# Perform clustering
feature_matrix <- features_clean[, c("slope", "variance", "early_late_diff", "peak_timing", "monotonicity")]
feature_matrix_scaled <- scale(feature_matrix)
complete_rows <- complete.cases(feature_matrix_scaled)
feature_matrix_clean <- feature_matrix_scaled[complete_rows, ]
features_for_clustering <- features_clean[complete_rows, ]

# Determine optimal clusters
wss <- sapply(1:10, function(k) {
  set.seed(123)
  kmeans(feature_matrix_clean, centers = k, nstart = 10)$tot.withinss
})
optimal_k <- which.min(diff(diff(wss))) + 1
if(optimal_k < 3) optimal_k <- 4

set.seed(123)
kmeans_result <- kmeans(feature_matrix_clean, centers = optimal_k, nstart = 25)
features_for_clustering$cluster <- kmeans_result$cluster

# Create clustered heatmap
clustered_orfs <- features_for_clustering %>%
  arrange(cluster, pattern_type, desc(abs(slope))) %>%
  pull(orf_name)

heatmap_clustered <- heatmap_matrix_clean[clustered_orfs, ]

cluster_annotation <- data.frame(
  Cluster = paste0("Cluster_", features_for_clustering$cluster[match(clustered_orfs, features_for_clustering$orf_name)]),
  Pattern_Type = features_for_clustering$pattern_type[match(clustered_orfs, features_for_clustering$orf_name)]
)
rownames(cluster_annotation) <- clustered_orfs

cluster_colors <- rainbow(optimal_k, alpha = 0.8)
names(cluster_colors) <- paste0("Cluster_", 1:optimal_k)
pattern_colors <- c("increasing" = "#E31A1C", "decreasing" = "#1F78B4", 
                   "stable" = "#33A02C", "complex" = "#FF7F00")
annotation_colors <- list(Cluster = cluster_colors, Pattern_Type = pattern_colors)

png("temporal_pattern_clusters_US.png", width = 1400, height = 1200, res = 150)
pheatmap(heatmap_clustered,
         cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
         color = colors, breaks = seq(color_limits[1], color_limits[2], length.out = 101),
         annotation_row = cluster_annotation, annotation_colors = annotation_colors,
         main = paste("US Cohort: Pattern-Based Clustering -", nrow(heatmap_clustered), "ORFs in", optimal_k, "Temporal Programs\n",
                      "Red = Higher in CONTROL, Blue = Higher in CELIAC"),
         fontsize = 8, fontsize_row = 4, fontsize_col = 10,
         show_rownames = FALSE, show_colnames = TRUE, border_color = NA)
dev.off()

cat("Pattern clustering completed with", optimal_k, "clusters\n")

# ============================================================================
# PART 4: DIVERSITY ANALYSIS
# ============================================================================

cat("\n4. Running diversity analysis...\n")

# Calculate diversity metrics
diversity_metrics <- data.frame(sample_id = matching_samples, stringsAsFactors = FALSE)
diversity_metrics <- merge(diversity_metrics, metadata_filtered, by.x = "sample_id", by.y = "row.names")

abundance_matrix_t <- t(abundance_filtered)
for(i in 1:nrow(diversity_metrics)) {
  sample_abundance <- abundance_matrix_t[i, ]
  diversity_metrics$richness[i] <- sum(sample_abundance > 0)
  
  shannon_val <- vegan::diversity(sample_abundance, index = "shannon")
  diversity_metrics$shannon[i] <- ifelse(is.finite(shannon_val), shannon_val, 0)
  
  simpson_val <- vegan::diversity(sample_abundance, index = "simpson")
  diversity_metrics$simpson[i] <- ifelse(is.finite(simpson_val), simpson_val, 0)
  
  if(diversity_metrics$richness[i] > 1) {
    diversity_metrics$evenness[i] <- diversity_metrics$shannon[i] / log(diversity_metrics$richness[i])
  } else {
    diversity_metrics$evenness[i] <- 0
  }
  
  diversity_metrics$total_abundance[i] <- sum(sample_abundance)
  diversity_metrics$dominance[i] <- 1 - diversity_metrics$simpson[i]
}

diversity_metrics <- diversity_metrics %>%
  group_by(patientID) %>%
  mutate(viral_load_cv = ifelse(mean(total_abundance, na.rm = TRUE) > 0, 
                               sd(total_abundance, na.rm = TRUE) / mean(total_abundance, na.rm = TRUE), 
                               0)) %>%
  ungroup()

diversity_metrics$viral_load_cv[!is.finite(diversity_metrics$viral_load_cv)] <- 0
diversity_metrics$evenness[!is.finite(diversity_metrics$evenness)] <- 0

# Diversity trajectory analysis
trajectory_design <- model.matrix(~ Dx.Status * onset_timeline_numeric + Sex + 
                                  Age.at.Gluten.Introduction..months. + HLA.Category + 
                                  feeding_first_year + Delivery.Mode, 
                                  data = diversity_metrics)

diversity_matrix <- as.matrix(t(diversity_metrics[, c("richness", "shannon", "simpson", 
                                                     "evenness", "total_abundance", 
                                                     "dominance", "viral_load_cv")]))
diversity_matrix[!is.finite(diversity_matrix)] <- 0

diversity_fit <- lmFit(diversity_matrix, trajectory_design)
diversity_fit <- eBayes(diversity_fit)

trajectory_results <- topTable(diversity_fit, coef = "Dx.StatusCONTROL:onset_timeline_numeric", 
                               number = Inf, sort.by = "P")

# Create diversity trajectory heatmap
png("diversity_trajectory_heatmap_US.png", width = 800, height = 600, res = 150)
heatmap_matrix_div <- as.matrix(trajectory_results$logFC)
rownames(heatmap_matrix_div) <- rownames(trajectory_results)
colnames(heatmap_matrix_div) <- "Trajectory_Interaction"

pheatmap(heatmap_matrix_div,
         cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "US Cohort: Diversity Trajectory Interactions (Dx.Status × Timeline)",
         fontsize = 12, fontsize_row = 10, show_colnames = FALSE, border_color = "white")
dev.off()

# Slope analysis
slope_data <- diversity_metrics %>%
  group_by(patientID) %>%
  summarise(
    Dx.Status = first(Dx.Status),
    Sex = first(Sex),
    Age.at.Gluten.Introduction..months. = first(Age.at.Gluten.Introduction..months.),
    HLA.Category = first(HLA.Category),
    feeding_first_year = first(feeding_first_year),
    Delivery.Mode = first(Delivery.Mode),
    richness_slope = if(n() >= 3) coef(lm(richness ~ onset_timeline_numeric))[2] else NA,
    shannon_slope = if(n() >= 3) coef(lm(shannon ~ onset_timeline_numeric))[2] else NA,
    simpson_slope = if(n() >= 3) coef(lm(simpson ~ onset_timeline_numeric))[2] else NA,
    evenness_slope = if(n() >= 3) coef(lm(evenness ~ onset_timeline_numeric))[2] else NA,
    total_abundance_slope = if(n() >= 3) coef(lm(total_abundance ~ onset_timeline_numeric))[2] else NA,
    dominance_slope = if(n() >= 3) coef(lm(dominance ~ onset_timeline_numeric))[2] else NA,
    viral_load_cv_slope = if(n() >= 3) coef(lm(viral_load_cv ~ onset_timeline_numeric))[2] else NA,
    .groups = "drop"
  ) %>%
  filter(!is.na(richness_slope))

slope_design <- model.matrix(~ Dx.Status + Sex + Age.at.Gluten.Introduction..months. + 
                            HLA.Category + feeding_first_year + Delivery.Mode, 
                            data = slope_data)

slope_matrix <- as.matrix(t(slope_data[, c("richness_slope", "shannon_slope", "simpson_slope",
                                          "evenness_slope", "total_abundance_slope", 
                                          "dominance_slope", "viral_load_cv_slope")]))

slope_fit <- lmFit(slope_matrix, slope_design)
slope_fit <- eBayes(slope_fit)
slope_results <- topTable(slope_fit, coef = "Dx.StatusCONTROL", number = Inf, sort.by = "P")

# Create slope analysis plots
slope_plot_data <- slope_results
slope_plot_data$metric <- rownames(slope_plot_data)
slope_plot_data$significant <- slope_plot_data$adj.P.Val < 0.05

p_slope_bar <- ggplot(slope_plot_data, aes(x = reorder(metric, logFC), y = logFC)) +
  geom_col(aes(fill = significant), alpha = 0.8) +
  geom_text(aes(label = paste0("p=", round(adj.P.Val, 3))), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
  labs(title = "US Cohort: Slope Analysis Results",
       subtitle = "Effect sizes (logFC) for CONTROL vs CELIAC slope differences",
       x = "Diversity Metric", y = "Effect Size (logFC)") +
  theme_minimal() + theme(legend.position = "none")

ggsave("slope_analysis_results_US.png", p_slope_bar, width = 10, height = 6, dpi = 300)

# Create mean trajectories plot
trajectory_summary <- diversity_metrics %>%
  group_by(Dx.Status, onset_timeline_numeric) %>%
  summarise(
    richness_mean = mean(richness, na.rm = TRUE),
    shannon_mean = mean(shannon, na.rm = TRUE),
    simpson_mean = mean(simpson, na.rm = TRUE),
    evenness_mean = mean(evenness, na.rm = TRUE),
    .groups = "drop"
  )

metrics_to_plot <- c("richness_mean", "shannon_mean", "simpson_mean", "evenness_mean")
plot_list <- list()

for(metric in metrics_to_plot) {
  metric_name <- gsub("_mean", "", metric)
  p <- ggplot(trajectory_summary, aes(x = onset_timeline_numeric, y = .data[[metric]], color = Dx.Status)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = c("CELIAC" = "#E31A1C", "CONTROL" = "#1F78B4")) +
    labs(title = paste("US Cohort: Mean", str_to_title(metric_name), "Trajectory"),
         x = "Months Relative to Diagnosis", y = str_to_title(metric_name)) +
    theme_minimal() + theme(legend.position = "bottom")
  plot_list[[metric]] <- p
}

combined_trajectories <- do.call(grid.arrange, c(plot_list, ncol = 2))
ggsave("diversity_trajectories_by_group_US.png", combined_trajectories, width = 12, height = 8, dpi = 300)

# Create summary statistics
summary_stats <- data.frame(
  Metric = c("Total Samples", "Total Patients", "Significant ORFs (Interaction)", "Pattern Clusters", 
             "Trajectory Interactions", "Significant Slopes"),
  Count = c(nrow(diversity_metrics), length(unique(diversity_metrics$patientID)), 
            nrow(significant_orfs), optimal_k, sum(trajectory_results$adj.P.Val < 0.05), 
            sum(slope_results$adj.P.Val < 0.05)),
  Category = c("Dataset", "Dataset", "ORF Analysis", "Pattern Analysis", 
               "Trajectory Analysis", "Slope Analysis")
)

p_summary <- ggplot(summary_stats, aes(x = reorder(Metric, Count), y = Count, fill = Category)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.1, size = 4) +
  coord_flip() +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(title = "US Cohort: Analysis Summary Statistics",
       x = "Analysis Component", y = "Count") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("analysis_summary_statistics_US.png", p_summary, width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 5: SAVE ALL DATA FILES
# ============================================================================

cat("\n5. Saving all analysis data files...\n")

write.csv(diversity_metrics, "diversity_data_full_US.csv", row.names = FALSE)
write.csv(trajectory_results, "diversity_trajectory_results_US.csv", row.names = TRUE)
write.csv(slope_results, "slope_analysis_results_US.csv", row.names = TRUE)
write.csv(slope_data, "slope_data_US.csv", row.names = FALSE)
write.csv(features_for_clustering, "pattern_features_US.csv", row.names = FALSE)
write.csv(summary_stats, "analysis_summary_US.csv", row.names = FALSE)

# ============================================================================
# PART 6: FINAL SUMMARY
# ============================================================================

cat("\n=== US COHORT COMPREHENSIVE ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- US_limma_model_Sig_res_final.csv: Significant ORFs results\n")
cat("- all_significant_orfs_temporal_heatmap_US_limma_fitted.png: All significant ORFs heatmap\n")
cat("- top_variable_significant_orfs_heatmap_US_limma_fitted.png: Top variable ORFs heatmap\n")
cat("- temporal_pattern_clusters_US.png: Pattern clustering analysis\n")
cat("- diversity_trajectory_heatmap_US.png: Diversity trajectory interactions\n")
cat("- slope_analysis_results_US.png: Slope analysis results\n")
cat("- diversity_trajectories_by_group_US.png: Mean diversity trajectories\n")
cat("- analysis_summary_statistics_US.png: Summary statistics\n")
cat("- All supporting data files (.csv)\n")

cat("\nKey Results Summary:\n")
cat("- Total samples:", nrow(diversity_metrics), "\n")
cat("- Total patients:", length(unique(diversity_metrics$patientID)), "\n")
cat("- Significant ORFs:", nrow(significant_orfs), "(", round(nrow(significant_orfs)/nrow(abundance_matrix)*100, 1), "% of total)\n")
cat("- Pattern clusters identified:", optimal_k, "\n")
cat("- Significant trajectory interactions:", sum(trajectory_results$adj.P.Val < 0.05), "\n")
cat("- Significant slope differences:", sum(slope_results$adj.P.Val < 0.05), "\n")
cat("- Data range (heatmap):", round(min(heatmap_matrix_clean, na.rm = TRUE), 2), "to", round(max(heatmap_matrix_clean, na.rm = TRUE), 2), "\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")