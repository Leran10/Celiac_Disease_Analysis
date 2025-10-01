#!/usr/bin/env Rscript
# Recreate Volcano Plot with Clean White Background (No Black Margins)
# Author: Claude AI
# Date: 2025-08-05

library(dplyr)
library(ggplot2)

cat("=== RECREATING CLEAN VOLCANO PLOT ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading data...\n")

# Load limma results
limma_results <- read.csv("total_limma_model_Sig_res_final.csv", row.names = 1)
cat("Loaded limma results with", nrow(limma_results), "significant ORFs\n")

# Load functional annotations
phold_data <- read.csv("phold_per_cds_predictions.tsv", sep = "\t")
cat("Loaded phold predictions with", nrow(phold_data), "entries\n")

# ============================================================================
# PART 2: MERGE DATA
# ============================================================================

cat("\n2. Merging data and preparing for visualization...\n")

# Merge limma results with phold data
merged_data <- merge(limma_results, phold_data, by.x = "row.names", by.y = "cds_id", all.x = TRUE)
rownames(merged_data) <- merged_data$Row.names
merged_data <- merged_data[, -1]  # Remove Row.names column

cat("Merged data has", nrow(merged_data), "rows\n")

# Check annotation coverage
annotated_count <- sum(!is.na(merged_data$phrog))
cat("Functional annotations available for", annotated_count, "out of", nrow(merged_data), 
    "significant ORFs\n")

# ============================================================================
# PART 3: CREATE CLEAN VOLCANO PLOT
# ============================================================================

cat("\n3. Creating clean volcano plot by function categories...\n")

# Create function categories - focus on KNOWN functions only
merged_data$function_category <- ifelse(is.na(merged_data$`function.`) | merged_data$`function.` == "", 
                                       NA,  # Set unknown to NA instead of labeling
                                       as.character(merged_data$`function.`))

# Filter to only known functions for this plot
known_functions_data <- merged_data[!is.na(merged_data$function_category), ]

cat("Using", nrow(known_functions_data), "ORFs with known functional annotations\n")

if(nrow(known_functions_data) > 0) {
  # Get top function categories
  function_counts <- table(known_functions_data$function_category)
  top_functions <- names(sort(function_counts, decreasing = TRUE))[1:20]
  
  # Group rare functions
  known_functions_data$function_category[!known_functions_data$function_category %in% top_functions] <- "other functions"
  
  # Create color palette
  unique_functions <- unique(known_functions_data$function_category)
  n_functions <- length(unique_functions)
  
  # Create distinct colors for functions
  function_colors <- rainbow(n_functions, alpha = 0.8)
  names(function_colors) <- unique_functions
  
  # Set specific color for "other functions"
  if("other functions" %in% names(function_colors)) {
    function_colors["other functions"] <- "grey60"
  }
  
  # Create the clean volcano plot
  p <- ggplot(known_functions_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = function_category), size = 2, alpha = 0.8) +
    scale_color_manual(values = function_colors) +
    labs(
      title = "Volcano Plot: Significant ORFs by Known Functional Categories",
      subtitle = paste("Showing", nrow(known_functions_data), "ORFs with functional annotations"),
      x = "Effect Size (logFC)\nPositive = Higher in CONTROL, Negative = Higher in CELIAC",
      y = "-log10(adjusted P-value)",
      color = "Functional Category"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid.major = element_line(color = "grey95", size = 0.5),
      panel.grid.minor = element_line(color = "grey98", size = 0.25),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.background = element_rect(fill = "white", color = "white"),
      legend.box.background = element_rect(fill = "white", color = "white"),
      plot.margin = margin(20, 20, 20, 20),  # Add margins for clean edges
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  
  # Save as PNG with white background and no margins
  ggsave("volcano_plot_functional_categories_clean.png", p, 
         width = 14, height = 10, dpi = 300, 
         bg = "white",  # Ensure white background
         device = "png")
  
  cat("Clean volcano plot saved as PNG\n")
  
  # Also save as PDF with white background
  ggsave("volcano_plot_functional_categories_clean.pdf", p, 
         width = 14, height = 10, 
         bg = "white",  # Ensure white background
         device = "pdf")
  
  cat("Clean volcano plot also saved as PDF\n")
  
} else {
  cat("No ORFs with known functional annotations found!\n")
}

# ============================================================================
# PART 4: SUMMARY
# ============================================================================

cat("\n=== CLEAN VOLCANO PLOT COMPLETED ===\n")
cat("Generated files:\n")
cat("- volcano_plot_functional_categories_clean.png: Clean PNG version (white background, no margins)\n")
cat("- volcano_plot_functional_categories_clean.pdf: Clean PDF version (white background, no margins)\n")

if(exists("known_functions_data")) {
  cat("- ORFs with known functions:", nrow(known_functions_data), "\n")
  cat("- Function categories displayed:", length(unique(known_functions_data$function_category)), "\n")
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")