#!/usr/bin/env Rscript
# Create Volcano Plot with Taxonomic Coloring - Italy Cohort
# Author: Claude AI
# Date: 2025-08-06

library(dplyr)
library(ggplot2)
library(stringr)

cat("=== CREATING TAXONOMY VOLCANO PLOT - ITALY COHORT ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading data files...\n")

# Load limma results - Italy cohort
limma_results <- read.csv("Italy.limma_model_res.csv", row.names = 1)
cat("Loaded Italy cohort limma results with", nrow(limma_results), "ORFs\n")

# Load taxonomy data (same file as total cohort)
taxonomy_data <- read.table("../../mmseqs_taxonomy.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(taxonomy_data) <- c("orf_id_short", "taxid", "rank", "taxonomy", "score", "matches", "aligned", "coverage", "lineage")
cat("Loaded taxonomy data with", nrow(taxonomy_data), "entries\n")

# Filter for significant results only
significant_results <- limma_results[limma_results$adj.P.Val < 0.05, ]
cat("Filtered to", nrow(significant_results), "significant ORFs (adj.p < 0.05)\n")

# Check if we have any significant results
if(nrow(significant_results) == 0) {
  cat("WARNING: No significant ORFs found at adj.p < 0.05 threshold\n")
  cat("Proceeding with all ORFs for visualization purposes\n")
  significant_results <- limma_results
}

# ============================================================================
# PART 2: PROCESS ORF IDs
# ============================================================================

cat("\n2. Processing ORF IDs for matching...\n")

# Create modified ORF ID column in limma results by removing last part
significant_results$orf_id_original <- rownames(significant_results)
significant_results$orf_id_short <- gsub("_[0-9]+$", "", significant_results$orf_id_original)

cat("Example ORF ID transformations:\n")
examples <- head(significant_results[, c("orf_id_original", "orf_id_short")], 5)
for(i in 1:nrow(examples)) {
  cat(paste0("  ", examples$orf_id_original[i], " -> ", examples$orf_id_short[i], "\n"))
}

# Check matching between datasets
matching_ids <- intersect(significant_results$orf_id_short, taxonomy_data$orf_id_short)
cat("\nFound", length(matching_ids), "matching ORF IDs between datasets\n")
cat("ORFs in analysis:", nrow(significant_results), "ORFs\n")
cat("Taxonomy data:", nrow(taxonomy_data), "entries\n")
cat("Overlap:", length(matching_ids), "ORFs\n")

# ============================================================================
# PART 3: MERGE DATA
# ============================================================================

cat("\n3. Merging datasets...\n")

# Merge the datasets
merged_data <- merge(significant_results, taxonomy_data, by = "orf_id_short", all.x = TRUE)
cat("Merged dataset contains", nrow(merged_data), "ORFs\n")
cat("ORFs with taxonomy annotation:", sum(!is.na(merged_data$taxonomy)), "\n")
cat("ORFs without taxonomy annotation:", sum(is.na(merged_data$taxonomy)), "\n")

# Check taxonomy categories
if(sum(!is.na(merged_data$taxonomy)) > 0) {
  taxonomy_counts <- table(merged_data$taxonomy, useNA = "ifany")
  cat("\nTaxonomy categories found:\n")
  print(taxonomy_counts)
}

# ============================================================================
# PART 4: CREATE VOLCANO PLOT
# ============================================================================

cat("\n4. Creating volcano plot with taxonomic coloring...\n")

# Prepare data for plotting
plot_data <- merged_data
plot_data$logFC_flipped <- -plot_data$logFC  # Flip logFC to move significant dots to the left
plot_data$neg_log10_pval <- -log10(plot_data$adj.P.Val)

# Handle taxonomy categories - group rare ones
if(sum(!is.na(plot_data$taxonomy)) > 0) {
  taxonomy_counts <- table(plot_data$taxonomy, useNA = "ifany")
  
  # Keep top taxonomies, group others as "Other"
  top_taxonomies <- names(sort(taxonomy_counts, decreasing = TRUE))[1:min(8, length(taxonomy_counts))]
  
  plot_data$taxonomy_grouped <- ifelse(
    is.na(plot_data$taxonomy), 
    "Unknown",
    ifelse(plot_data$taxonomy %in% top_taxonomies, 
           plot_data$taxonomy, 
           "Other")
  )
} else {
  plot_data$taxonomy_grouped <- "Unknown"
}

# Count final categories
final_counts <- table(plot_data$taxonomy_grouped)
cat("\nFinal taxonomy groups for plotting:\n")
print(final_counts)

# Create annotation text with counts for each category
annotation_text <- paste0(names(final_counts), " (n=", final_counts, ")", collapse = "\n")
cat("\nAnnotation text to be added:\n")
cat(annotation_text)

# Create color palette
unique_taxonomies <- unique(plot_data$taxonomy_grouped)
n_taxonomies <- length(unique_taxonomies)

# Create distinct colors
if(n_taxonomies <= 8) {
  colors <- RColorBrewer::brewer.pal(max(3, n_taxonomies), "Set2")
} else {
  colors <- rainbow(n_taxonomies, alpha = 0.8)
}
names(colors) <- unique_taxonomies

# Make sure "Unknown" is grey if present
if("Unknown" %in% names(colors)) {
  colors["Unknown"] <- "grey60"
}

# Make "root" light gray so other categories are more obvious
if("root" %in% names(colors)) {
  colors["root"] <- "lightgray"
}

# Create annotation subtitle before creating the plot
annotation_subtitle <- paste(
  paste0(names(final_counts), ": ", final_counts), 
  collapse = " | "
)

# Determine title based on whether we have significant results
sig_count <- sum(significant_results$adj.P.Val < 0.05, na.rm = TRUE)
if(sig_count > 0) {
  plot_title <- "Volcano Plot: Italy Cohort ORFs Colored by Taxonomy"
  plot_subtitle <- paste("Showing", sig_count, "significant ORFs with taxonomic annotations\n", annotation_subtitle)
} else {
  plot_title <- "Volcano Plot: Italy Cohort ORFs Colored by Taxonomy"
  plot_subtitle <- paste("Showing all", nrow(plot_data), "ORFs (no significant results at adj.p < 0.05)\n", annotation_subtitle)
}

# Create the volcano plot
p <- ggplot(plot_data, aes(x = logFC_flipped, y = neg_log10_pval)) +
  geom_point(aes(color = taxonomy_grouped), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = colors) +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "Effect Size (-logFC)\nNegative = Higher in CONTROL, Positive = Higher in CELIAC",
    y = "-log10(adjusted P-value)",
    color = "Taxonomy"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey98", linewidth = 0.25)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# Add significance lines
p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7)
p <- p + geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred", alpha = 0.7)

# Save the plot - increase width to prevent text trimming
ggsave("taxonomy_volcano_plot_Italy.png", p, width = 14, height = 9, dpi = 300, bg = "white")
cat("Volcano plot saved as 'taxonomy_volcano_plot_Italy.png'\n")

# Also save as PDF
ggsave("taxonomy_volcano_plot_Italy.pdf", p, width = 14, height = 9, bg = "white")
cat("Volcano plot also saved as 'taxonomy_volcano_plot_Italy.pdf'\n")

# ============================================================================
# PART 5: SAVE MERGED DATA AND SUMMARY
# ============================================================================

cat("\n5. Saving merged data and summary statistics...\n")

# Save merged data
write.csv(merged_data, "limma_results_with_taxonomy_Italy.csv", row.names = FALSE)
cat("Merged data saved as 'limma_results_with_taxonomy_Italy.csv'\n")

# Create summary statistics
summary_stats <- data.frame(
  Category = c("Total ORFs in limma results", "Significant ORFs (adj.p < 0.05)", "ORFs with taxonomy annotation", 
               "ORFs without taxonomy annotation", "Unique taxonomies", 
               "Most common taxonomy", "Significant ORFs (adj.p < 0.01)"),
  Count = c(nrow(limma_results), sum(limma_results$adj.P.Val < 0.05, na.rm = TRUE),
            sum(!is.na(merged_data$taxonomy)),
            sum(is.na(merged_data$taxonomy)),
            length(unique(merged_data$taxonomy[!is.na(merged_data$taxonomy)])),
            if(sum(!is.na(merged_data$taxonomy)) > 0) names(sort(table(merged_data$taxonomy), decreasing = TRUE))[1] else "None",
            sum(merged_data$adj.P.Val < 0.01, na.rm = TRUE))
)

write.csv(summary_stats, "taxonomy_analysis_summary_Italy.csv", row.names = FALSE)
cat("Summary statistics saved as 'taxonomy_analysis_summary_Italy.csv'\n")

# Print summary
cat("\n=== TAXONOMY VOLCANO PLOT ANALYSIS COMPLETED - ITALY COHORT ===\n")
cat("Generated files:\n")
cat("- taxonomy_volcano_plot_Italy.png: Main volcano plot with taxonomic coloring\n")
cat("- taxonomy_volcano_plot_Italy.pdf: PDF version of the plot\n")
cat("- limma_results_with_taxonomy_Italy.csv: Merged dataset with taxonomy annotations\n")
cat("- taxonomy_analysis_summary_Italy.csv: Summary statistics\n")

cat("\nSummary Statistics:\n")
print(summary_stats)

cat("\nAnalysis completed at:", format(Sys.time()), "\n")