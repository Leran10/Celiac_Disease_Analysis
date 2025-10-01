#!/usr/bin/env Rscript
# Create Volcano Plot with Taxonomic Coloring
# Author: Claude AI
# Date: 2025-08-05

library(dplyr)
library(ggplot2)
library(stringr)

cat("=== CREATING TAXONOMY VOLCANO PLOT ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading data files...\n")

# Load limma results - use the SIGNIFICANT ORFs file for better volcano plot
limma_results <- read.csv("total_limma_model_Sig_res_final.csv", row.names = 1)
cat("Loaded significant limma results with", nrow(limma_results), "ORFs\n")

# Load taxonomy data
taxonomy_data <- read.table("../../mmseqs_taxonomy.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(taxonomy_data) <- c("orf_id_short", "taxid", "rank", "taxonomy", "score", "matches", "aligned", "coverage", "lineage")
cat("Loaded taxonomy data with", nrow(taxonomy_data), "entries\n")

# ============================================================================
# PART 2: PROCESS ORF IDs
# ============================================================================

cat("\n2. Processing ORF IDs for matching...\n")

# Create modified ORF ID column in limma results by removing last part
limma_results$orf_id_original <- rownames(limma_results)
limma_results$orf_id_short <- gsub("_[0-9]+$", "", limma_results$orf_id_original)

cat("Example ORF ID transformations:\n")
examples <- head(limma_results[, c("orf_id_original", "orf_id_short")], 5)
for(i in 1:nrow(examples)) {
  cat(paste0("  ", examples$orf_id_original[i], " -> ", examples$orf_id_short[i], "\n"))
}

# Check matching between datasets
matching_ids <- intersect(limma_results$orf_id_short, taxonomy_data$orf_id_short)
cat("\nFound", length(matching_ids), "matching ORF IDs between datasets\n")
cat("Limma results:", nrow(limma_results), "ORFs\n")
cat("Taxonomy data:", nrow(taxonomy_data), "entries\n")
cat("Overlap:", length(matching_ids), "ORFs\n")

# ============================================================================
# PART 3: MERGE DATA
# ============================================================================

cat("\n3. Merging datasets...\n")

# Merge the datasets
merged_data <- merge(limma_results, taxonomy_data, by = "orf_id_short", all.x = TRUE)
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

# Create the volcano plot
p <- ggplot(plot_data, aes(x = logFC_flipped, y = neg_log10_pval)) +
  geom_point(aes(color = taxonomy_grouped), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = colors) +
  labs(
    title = "Volcano Plot: Total Cohort ORFs Colored by Taxonomy",
    subtitle = paste("Showing", nrow(plot_data), "significant ORFs with taxonomic annotations\n", annotation_subtitle),
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
    panel.grid.major = element_line(color = "grey95", size = 0.5),
    panel.grid.minor = element_line(color = "grey98", size = 0.25)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))


# Add significance lines if there are significant results
if(any(plot_data$adj.P.Val < 0.05)) {
  p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7)
}
if(any(plot_data$adj.P.Val < 0.01)) {
  p <- p + geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred", alpha = 0.7)
}

# Save the plot - increase width to prevent text trimming
ggsave("taxonomy_volcano_plot.png", p, width = 14, height = 9, dpi = 300, bg = "white")
cat("Volcano plot saved as 'taxonomy_volcano_plot.png'\n")

# Also save as PDF
ggsave("taxonomy_volcano_plot.pdf", p, width = 14, height = 9, bg = "white")
cat("Volcano plot also saved as 'taxonomy_volcano_plot.pdf'\n")

# ============================================================================
# PART 5: SAVE MERGED DATA AND SUMMARY
# ============================================================================

cat("\n5. Saving merged data and summary statistics...\n")

# Save merged data
write.csv(merged_data, "limma_results_with_taxonomy.csv", row.names = FALSE)
cat("Merged data saved as 'limma_results_with_taxonomy.csv'\n")

# Create summary statistics
summary_stats <- data.frame(
  Category = c("Total ORFs in limma results", "ORFs with taxonomy annotation", 
               "ORFs without taxonomy annotation", "Unique taxonomies", 
               "Most common taxonomy", "Significant ORFs (adj.p < 0.05)",
               "Significant ORFs (adj.p < 0.01)"),
  Count = c(nrow(merged_data), 
            sum(!is.na(merged_data$taxonomy)),
            sum(is.na(merged_data$taxonomy)),
            length(unique(merged_data$taxonomy[!is.na(merged_data$taxonomy)])),
            if(sum(!is.na(merged_data$taxonomy)) > 0) names(sort(table(merged_data$taxonomy), decreasing = TRUE))[1] else "None",
            sum(merged_data$adj.P.Val < 0.05),
            sum(merged_data$adj.P.Val < 0.01))
)

write.csv(summary_stats, "taxonomy_analysis_summary.csv", row.names = FALSE)
cat("Summary statistics saved as 'taxonomy_analysis_summary.csv'\n")

# Print summary
cat("\n=== TAXONOMY VOLCANO PLOT ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- taxonomy_volcano_plot.png: Main volcano plot with taxonomic coloring\n")
cat("- taxonomy_volcano_plot.pdf: PDF version of the plot\n")
cat("- limma_results_with_taxonomy.csv: Merged dataset with taxonomy annotations\n")
cat("- taxonomy_analysis_summary.csv: Summary statistics\n")

cat("\nSummary Statistics:\n")
print(summary_stats)

cat("\nAnalysis completed at:", format(Sys.time()), "\n")