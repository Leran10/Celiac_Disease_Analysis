#!/usr/bin/env Rscript
# Categorical Timeline Analysis - US Cohort
# Author: Claude AI
# Date: 2025-08-11
# Purpose: Analyze ORF changes at specific timepoints using factor(onset_timeline_combined)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)
library(gridExtra)
library(ggrepel)

cat("=== CATEGORICAL TIMELINE ANALYSIS - US COHORT ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA AND SET UP CATEGORICAL MODEL
# ============================================================================

cat("1. Loading US cohort data and setting up categorical timeline model...\n")

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

cat("Loaded abundance data with", nrow(abundance_data), "ORFs and", ncol(abundance_data), "samples\n")
cat("Loaded metadata with", nrow(metadata), "samples\n")

# Match samples
matching_samples <- intersect(colnames(abundance_data), rownames(metadata))
abundance_filtered <- abundance_data[, matching_samples]
metadata_filtered <- metadata[matching_samples, ]

cat("Using", length(matching_samples), "samples for analysis\n")

# Set reference levels
metadata_filtered$Dx.Status <- factor(metadata_filtered$Dx.Status, levels = c("CELIAC", "CONTROL"))

# CRITICAL: Set t0 as reference level for onset_timeline_combined
cat("Setting t0 as reference level for onset_timeline_combined\n")
metadata_filtered$onset_timeline_combined <- factor(metadata_filtered$onset_timeline_combined, 
                                                   levels = c("t0", "t0-6", "t0-12", "t0-18", "t0-24", "t0-30", "t0-36", "t0-over42"))

# Check factor levels
cat("Timeline factor levels:", levels(metadata_filtered$onset_timeline_combined), "\n")
cat("Timeline distribution:\n")
print(table(metadata_filtered$onset_timeline_combined))

# Create design matrix with CATEGORICAL timeline
linear.model.design <- model.matrix(~ Dx.Status * onset_timeline_combined + Sex + 
                                   Age.at.Gluten.Introduction..months. + HLA.Category + 
                                   feeding_first_year + Delivery.Mode, 
                                   data = metadata_filtered)

cat("Design matrix dimensions:", nrow(linear.model.design), "x", ncol(linear.model.design), "\n")
cat("Design matrix column names:\n")
print(colnames(linear.model.design))

# Convert abundance data to matrix for limma
abundance_matrix <- as.matrix(abundance_filtered)

# ============================================================================
# PART 2: RUN LIMMA ANALYSIS WITH CATEGORICAL TIMELINE
# ============================================================================

cat("\n2. Running limma analysis with categorical timeline...\n")

# Fit limma model with voomLmFit and patient blocking
linear.model.fit <- voomLmFit(abundance_matrix, linear.model.design, 
                              block = metadata_filtered$patientID, plot = FALSE)
linear.model.fit <- eBayes(linear.model.fit)

cat("Model fitted successfully. Extracting results for each timepoint...\n")

# Get all coefficients related to onset_timeline_combined
timeline_coeffs <- colnames(linear.model.design)[grep("onset_timeline_combined", colnames(linear.model.design))]
cat("Timeline coefficients found:", paste(timeline_coeffs, collapse = ", "), "\n")

# ============================================================================
# PART 3: EXTRACT RESULTS FOR EACH TIMEPOINT
# ============================================================================

cat("\n3. Extracting significant ORFs for each timepoint...\n")

# Store results for each timepoint
timepoint_results <- list()
timepoint_summary <- data.frame()

for(coeff in timeline_coeffs) {
  cat("Processing coefficient:", coeff, "\n")
  
  # Get results for this timepoint
  results <- topTable(linear.model.fit, coef = coeff, number = Inf, sort.by = "P")
  
  # Get significant ORFs
  sig_orfs <- results[results$adj.P.Val < 0.05, ]
  
  # Store results
  timepoint_results[[coeff]] <- results
  
  # Create summary row
  timepoint_name <- gsub("onset_timeline_combined", "", coeff)
  summary_row <- data.frame(
    Timepoint = timepoint_name,
    Total_ORFs = nrow(results),
    Significant_ORFs = nrow(sig_orfs),
    Percent_Significant = round(nrow(sig_orfs)/nrow(results)*100, 1),
    Positive_Changes = sum(sig_orfs$logFC > 0),
    Negative_Changes = sum(sig_orfs$logFC < 0),
    Max_logFC = ifelse(nrow(sig_orfs) > 0, round(max(sig_orfs$logFC), 2), 0),
    Min_logFC = ifelse(nrow(sig_orfs) > 0, round(min(sig_orfs$logFC), 2), 0),
    stringsAsFactors = FALSE
  )
  
  timepoint_summary <- rbind(timepoint_summary, summary_row)
  
  cat("  -", timepoint_name, ":", nrow(sig_orfs), "significant ORFs\n")
}

# Print summary
cat("\n=== TIMEPOINT SUMMARY ===\n")
print(timepoint_summary)

# Save detailed results
write.csv(timepoint_summary, "categorical_timeline_summary.csv", row.names = FALSE)

# ============================================================================
# PART 4: CREATE VISUALIZATION PLOTS
# ============================================================================

cat("\n4. Creating visualization plots...\n")

# Plot 1: Number of significant ORFs by timepoint
p1 <- ggplot(timepoint_summary, aes(x = factor(Timepoint, levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6")), 
                                   y = Significant_ORFs)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = Significant_ORFs), vjust = -0.3, size = 3) +
  labs(title = "US Cohort: Number of Significant ORFs by Timepoint",
       subtitle = "ORFs significantly different compared to T0 (diagnosis)",
       x = "Timepoint", y = "Number of Significant ORFs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Percentage of significant ORFs
p2 <- ggplot(timepoint_summary, aes(x = factor(Timepoint, levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6")), 
                                   y = Percent_Significant)) +
  geom_col(fill = "darkgreen", alpha = 0.7) +
  geom_text(aes(label = paste0(Percent_Significant, "%")), vjust = -0.3, size = 3) +
  labs(title = "US Cohort: Percentage of ORFs Significantly Changed",
       subtitle = "Compared to T0 (diagnosis time)",
       x = "Timepoint", y = "Percentage of ORFs (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 3: Direction of changes (positive vs negative logFC)
timepoint_summary_long <- reshape2::melt(timepoint_summary[, c("Timepoint", "Positive_Changes", "Negative_Changes")], 
                                        id.vars = "Timepoint")

p3 <- ggplot(timepoint_summary_long, aes(x = factor(Timepoint, levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6")), 
                                        y = value, fill = variable)) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_manual(values = c("Positive_Changes" = "#E31A1C", "Negative_Changes" = "#1F78B4"),
                    labels = c("Higher at timepoint vs T0", "Lower at timepoint vs T0")) +
  labs(title = "US Cohort: Direction of ORF Changes by Timepoint",
       subtitle = "Positive = Higher than T0, Negative = Lower than T0",
       x = "Timepoint", y = "Number of ORFs", fill = "Change Direction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Plot 4: Effect size ranges
p4 <- ggplot(timepoint_summary, aes(x = factor(Timepoint, levels = c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6")))) +
  geom_point(aes(y = Max_logFC), color = "red", size = 3, alpha = 0.7) +
  geom_point(aes(y = Min_logFC), color = "blue", size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "US Cohort: Effect Size Ranges by Timepoint",
       subtitle = "Red = Maximum logFC, Blue = Minimum logFC",
       x = "Timepoint", y = "Effect Size (logFC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)

# Save plots
ggsave("categorical_timeline_analysis_plots.png", combined_plot, width = 14, height = 10, dpi = 300)

# ============================================================================
# PART 5: DETAILED ANALYSIS OF EACH TIMEPOINT
# ============================================================================

cat("\n5. Creating detailed analysis for each significant timepoint...\n")

# For each timepoint with significant ORFs, create detailed output
for(coeff in timeline_coeffs) {
  timepoint_name <- gsub("onset_timeline_combined", "", coeff)
  results <- timepoint_results[[coeff]]
  sig_orfs <- results[results$adj.P.Val < 0.05, ]
  
  if(nrow(sig_orfs) > 0) {
    # Save detailed results for this timepoint
    filename <- paste0("significant_orfs_", timepoint_name, ".csv")
    write.csv(sig_orfs, filename, row.names = TRUE)
    
    cat("Saved", nrow(sig_orfs), "significant ORFs for", timepoint_name, "to", filename, "\n")
  }
}

# ============================================================================
# PART 6: VOLCANO PLOTS FOR KEY TIMEPOINTS
# ============================================================================

cat("\n6. Creating volcano plots for key timepoints...\n")

# Function to create volcano plot
create_volcano_plot <- function(results, timepoint_name, title_suffix = "") {
  # Add significance categories
  results$Significance <- "Not Significant"
  results$Significance[results$adj.P.Val < 0.05 & results$logFC > 0] <- "Higher vs T0"
  results$Significance[results$adj.P.Val < 0.05 & results$logFC < 0] <- "Lower vs T0"
  
  # Create plot
  p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Not Significant" = "gray60", 
                                  "Higher vs T0" = "#E31A1C", 
                                  "Lower vs T0" = "#1F78B4")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
    labs(title = paste("US Cohort:", timepoint_name, "vs T0", title_suffix),
         subtitle = paste("Significant ORFs:", sum(results$adj.P.Val < 0.05)),
         x = "Effect Size (logFC)", y = "-log10(Adjusted P-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Create volcano plots for timepoints with most significant ORFs
top_timepoints <- timepoint_summary[order(timepoint_summary$Significant_ORFs, decreasing = TRUE)[1:4], "Timepoint"]

volcano_plots <- list()
for(tp in top_timepoints) {
  coeff_name <- paste0("onset_timeline_combined", tp)
  if(coeff_name %in% names(timepoint_results)) {
    results <- timepoint_results[[coeff_name]]
    volcano_plots[[tp]] <- create_volcano_plot(results, tp)
  }
}

# Save volcano plots
if(length(volcano_plots) > 0) {
  combined_volcano <- do.call(grid.arrange, c(volcano_plots, ncol = 2))
  ggsave("categorical_timeline_volcano_plots.png", combined_volcano, width = 12, height = 10, dpi = 300)
}

# ============================================================================
# PART 7: FINAL SUMMARY
# ============================================================================

cat("\n=== CATEGORICAL TIMELINE ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- categorical_timeline_summary.csv: Summary of significant ORFs by timepoint\n")
cat("- categorical_timeline_analysis_plots.png: Overview visualization plots\n")
cat("- categorical_timeline_volcano_plots.png: Volcano plots for top timepoints\n")
cat("- significant_orfs_[timepoint].csv: Detailed results for each timepoint with significant ORFs\n")

cat("\nKey Findings Summary:\n")
most_changed_timepoint <- timepoint_summary[which.max(timepoint_summary$Significant_ORFs), ]
cat("- Timepoint with most changes:", most_changed_timepoint$Timepoint, 
    "(", most_changed_timepoint$Significant_ORFs, "significant ORFs)\n")
cat("- Highest percentage changed:", max(timepoint_summary$Percent_Significant), "% at", 
    timepoint_summary[which.max(timepoint_summary$Percent_Significant), "Timepoint"], "\n")

# Check for progressive pattern
if(nrow(timepoint_summary) > 2) {
  early_changes <- mean(timepoint_summary[timepoint_summary$Timepoint %in% c("t0-over42", "t0-36", "t0-30"), "Significant_ORFs"])
  late_changes <- mean(timepoint_summary[timepoint_summary$Timepoint %in% c("t0-12", "t0-6"), "Significant_ORFs"])
  
  if(late_changes > early_changes) {
    cat("- Pattern: Progressive increase in changes approaching diagnosis\n")
  } else {
    cat("- Pattern: More changes detected in earlier timepoints\n")
  }
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")