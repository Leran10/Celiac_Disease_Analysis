#!/usr/bin/env Rscript
# Create volcano plot for limma differential abundance results
# Highlights the 1,845 significant ORFs from trajectory analysis
# Author: Claude AI
# Date: 2025-07-21

library(ggplot2)
library(dplyr)
library(ggrepel)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/compositonal_analysis")

cat("=== LIMMA VOLCANO PLOT GENERATION ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading limma results...\n")

# Load limma results
limma_results <- read.csv("../total_limma_results.csv", row.names = 1)

cat("Limma results loaded:", nrow(limma_results), "ORFs\n")
cat("Columns:", paste(colnames(limma_results), collapse = ", "), "\n")

# Check for required columns
required_cols <- c("logFC", "adj.P.Val", "P.Value", "t", "AveExpr")
missing_cols <- setdiff(required_cols, colnames(limma_results))
if(length(missing_cols) > 0) {
  cat("Warning: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
}

# ============================================================================
# PART 2: PREPARE DATA FOR VOLCANO PLOT
# ============================================================================

cat("\n2. Preparing volcano plot data...\n")

# Add significance categories
limma_results$significant <- limma_results$adj.P.Val < 0.05
limma_results$highly_significant <- limma_results$adj.P.Val < 0.001
limma_results$neg_log10_padj <- -log10(limma_results$adj.P.Val)
limma_results$neg_log10_pval <- -log10(limma_results$P.Value)

# Add ORF IDs as a column
limma_results$orf_id <- rownames(limma_results)

# Create significance labels
limma_results$significance_label <- "Not Significant"
limma_results$significance_label[limma_results$significant] <- "Significant (adj.p < 0.05)"
limma_results$significance_label[limma_results$highly_significant] <- "Highly Significant (adj.p < 0.001)"

# Convert to factor for proper ordering in legend
limma_results$significance_label <- factor(limma_results$significance_label,
                                          levels = c("Not Significant", 
                                                    "Significant (adj.p < 0.05)",
                                                    "Highly Significant (adj.p < 0.001)"))

# Summary statistics
cat("Total ORFs:", nrow(limma_results), "\n")
cat("Significant ORFs (adj.p < 0.05):", sum(limma_results$significant), 
    "(", round(100 * sum(limma_results$significant) / nrow(limma_results), 1), "%)\n")
cat("Highly significant ORFs (adj.p < 0.001):", sum(limma_results$highly_significant),
    "(", round(100 * sum(limma_results$highly_significant) / nrow(limma_results), 1), "%)\n")

cat("LogFC range:", round(range(limma_results$logFC, na.rm = TRUE), 3), "\n")
cat("Adj.p range:", range(limma_results$adj.P.Val, na.rm = TRUE), "\n")

# ============================================================================
# PART 3: CREATE VOLCANO PLOT
# ============================================================================

cat("\n3. Creating volcano plot...\n")

# Define colors for significance levels
colors <- c("Not Significant" = "grey70",
           "Significant (adj.p < 0.05)" = "#E31A1C", 
           "Highly Significant (adj.p < 0.001)" = "#B2182B")

# Create the volcano plot
volcano_plot <- ggplot(limma_results, aes(x = logFC, y = neg_log10_padj)) +
  geom_point(aes(color = significance_label), 
             size = 0.8, alpha = 0.7) +
  scale_color_manual(values = colors, name = "Significance") +
  
  # Add significance threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "red", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", 
             color = "darkred", alpha = 0.7) +
  
  # Customize axes
  labs(
    title = "Viral ORF Differential Abundance: CELIAC vs CONTROL",
    subtitle = paste("Total:", nrow(limma_results), "ORFs |", 
                     "Significant:", sum(limma_results$significant),
                     "ORFs (adj.p < 0.05) |",
                     "Highly Significant:", sum(limma_results$highly_significant), 
                     "ORFs (adj.p < 0.001)"),
    x = "Log2 Fold Change (CELIAC vs CONTROL)",
    y = "-Log10(Adjusted P-value)",
    caption = "Red dashed lines: adj.p = 0.05 and 0.001 thresholds\nPositive logFC = Higher in CONTROL, Lower in CELIAC"
  ) +
  
  # Theme customization
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 9, color = "grey50"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  
  # Add some padding to the plot
  expand_limits(y = 0)

# ============================================================================
# PART 4: CREATE ENHANCED VERSION WITH TOP GENE LABELS
# ============================================================================

cat("4. Creating enhanced volcano plot with top gene labels...\n")

# Identify top genes to label (most significant)
top_genes <- limma_results %>%
  filter(significant) %>%
  arrange(adj.P.Val) %>%
  head(20)

cat("Top 20 most significant genes selected for labeling\n")

# Enhanced volcano plot with labels
enhanced_volcano_plot <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = orf_id),
    size = 3,
    color = "black",
    bg.color = "white",
    bg.r = 0.1,
    box.padding = 0.3,
    point.padding = 0.3,
    force = 2,
    max.overlaps = 20,
    seed = 42
  ) +
  labs(
    title = "Viral ORF Differential Abundance with Top 20 Significant Genes",
    subtitle = paste("Total:", nrow(limma_results), "ORFs |", 
                     "Significant:", sum(limma_results$significant),
                     "ORFs (adj.p < 0.05) |",
                     "Top 20 most significant genes labeled")
  )

# ============================================================================
# PART 5: SAVE PLOTS
# ============================================================================

cat("\n5. Saving volcano plots...\n")

# Save basic volcano plot
ggsave("limma_volcano_plot.png", volcano_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Save enhanced volcano plot
ggsave("limma_volcano_plot_with_labels.png", enhanced_volcano_plot,
       width = 14, height = 10, dpi = 300, bg = "white")

cat("Volcano plots saved:\n")
cat("- limma_volcano_plot.png: Basic volcano plot\n")
cat("- limma_volcano_plot_with_labels.png: Enhanced with gene labels\n")

# ============================================================================
# PART 6: GENERATE SUMMARY STATISTICS
# ============================================================================

cat("\n6. Generating summary statistics...\n")

# Create summary statistics table
summary_stats <- data.frame(
  Category = c("Total ORFs", "Significant ORFs (adj.p < 0.05)", 
               "Highly Significant (adj.p < 0.001)", "Non-significant ORFs"),
  Count = c(nrow(limma_results), 
           sum(limma_results$significant),
           sum(limma_results$highly_significant),
           sum(!limma_results$significant)),
  Percentage = c(100, 
                round(100 * sum(limma_results$significant) / nrow(limma_results), 1),
                round(100 * sum(limma_results$highly_significant) / nrow(limma_results), 1),
                round(100 * sum(!limma_results$significant) / nrow(limma_results), 1)),
  stringsAsFactors = FALSE
)

print(summary_stats)

# Save summary statistics
write.csv(summary_stats, "limma_volcano_plot_summary.csv", row.names = FALSE)

# Save top genes information
write.csv(top_genes[, c("orf_id", "logFC", "adj.P.Val", "P.Value", "t", "AveExpr")], 
         "limma_volcano_top20_genes.csv", row.names = FALSE)

cat("\nSummary files saved:\n")
cat("- limma_volcano_plot_summary.csv: Overall statistics\n") 
cat("- limma_volcano_top20_genes.csv: Top 20 most significant genes\n")

cat("\n=== VOLCANO PLOT ANALYSIS COMPLETED ===\n")
cat("Generated:", format(Sys.time()), "\n")