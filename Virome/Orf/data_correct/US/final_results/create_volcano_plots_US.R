#!/usr/bin/env Rscript
# Volcano Plots with Functional Annotations - US Cohort Only
# Author: Claude AI
# Date: 2025-07-30

library(dplyr)
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/US/final_results")

cat("=== VOLCANO PLOTS ANALYSIS - US COHORT ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading data...\n")

# Load limma results
limma_results <- read.csv("../US.limma_model_res.csv")
cat("Loaded limma results with", nrow(limma_results), "ORFs\n")

# Load functional annotations (from total analysis directory)
phold_data <- read.csv("../../total/final_results/phold_per_cds_predictions.tsv", sep = "\t")
cat("Loaded phold predictions with", nrow(phold_data), "entries\n")

# Filter significant ORFs (adj.p < 0.01)
sig_data <- limma_results[limma_results$adj.P.Val < 0.01, ]
cat("Found", nrow(sig_data), "significant ORFs at adj.p < 0.01\n")

# ============================================================================
# PART 2: MERGE WITH FUNCTIONAL ANNOTATIONS
# ============================================================================

cat("\n2. Merging with functional annotations...\n")

# Merge limma results with phold data
# Match by cds_id (phold) to X (limma)
merged_data <- merge(sig_data, phold_data, by.x = "X", by.y = "cds_id", all.x = TRUE)
cat("Merged data has", nrow(merged_data), "rows\n")

# Check annotation coverage
annotated_count <- sum(!is.na(merged_data$phrog))
cat("Functional annotations available for", annotated_count, "out of", nrow(merged_data), 
    "significant ORFs (", round(annotated_count/nrow(merged_data)*100, 1), "%)\n")

# Note: Using direct column access for "function" column due to R reserved word

# ============================================================================
# PART 3: VOLCANO PLOT BY FUNCTION CATEGORIES
# ============================================================================

cat("\n3. Creating volcano plot by function categories...\n")

create_function_volcano <- function(data) {
  
  # Create function categories (note: column is "function." due to R reserved word handling)
  func_col <- data[["function."]]
  data$function_category <- ifelse(is.na(func_col) | func_col == "", 
                                   "unknown function", 
                                   as.character(func_col))
  
  # Get top function categories
  function_counts <- table(data$function_category)
  top_functions <- names(sort(function_counts, decreasing = TRUE))[1:20]
  
  # Group rare functions
  data$function_category[!data$function_category %in% top_functions] <- "other functions"
  
  # Create color palette
  unique_functions <- unique(data$function_category)
  n_functions <- length(unique_functions)
  
  # Identify gray categories
  gray_function_categories <- grep("unknown|other", unique_functions, value = TRUE, ignore.case = TRUE)
  
  # Create colors
  function_colors <- rainbow(n_functions, alpha = 0.8)
  names(function_colors) <- unique_functions
  
  # Set gray colors for unknown/other
  for(gray_cat in gray_function_categories) {
    if(gray_cat %in% names(function_colors)) {
      function_colors[gray_cat] <- "grey60"
    }
  }
  
  # Reorder data so gray function points are plotted first
  data_ordered <- data[order(data$function_category %in% gray_function_categories, decreasing = TRUE), ]
  
  # Create volcano plot
  p <- ggplot(data_ordered, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = function_category, 
                   size = ifelse(function_category %in% gray_function_categories, 1.2, 1.8),
                   alpha = ifelse(function_category %in% gray_function_categories, 0.4, 0.9))) +
    scale_color_manual(values = function_colors) +
    scale_size_identity() +
    scale_alpha_identity() +
    labs(
      title = "Volcano Plot: Significant ORFs by Functional Categories (US Cohort)",
      subtitle = paste("Showing", nrow(data), "ORFs with adj.p < 0.01"),
      x = "Effect Size (logFC)\nPositive = Higher in CONTROL, Negative = Higher in CELIAC",
      y = "-log10(adjusted P-value)",
      color = "Functional Category"
    ) +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  return(p)
}

function_plot <- create_function_volcano(merged_data)
ggsave("volcano_plot_by_function_US.png", function_plot, width = 14, height = 10, dpi = 300, bg = "white")
cat("Function volcano plot saved\n")

# ============================================================================
# PART 4: VOLCANO PLOT BY PHROG CATEGORIES
# ============================================================================

cat("\n4. Creating volcano plot by phrog categories...\n")

create_phrog_volcano <- function(data) {
  
  # Create phrog categories
  data$phrog_category <- ifelse(is.na(data$phrog) | data$phrog == "", 
                                "unknown phrog", 
                                paste0("phrog_", data$phrog))
  
  # Get top phrog categories
  phrog_counts <- table(data$phrog_category)
  top_phrogs <- names(sort(phrog_counts, decreasing = TRUE))[1:20]
  
  # Group rare phrogs
  data$phrog_category[!data$phrog_category %in% top_phrogs] <- "other phrogs"
  
  # Create color palette
  unique_phrogs <- unique(data$phrog_category)
  n_phrogs <- length(unique_phrogs)
  
  # Identify gray categories
  gray_phrog_categories <- grep("unknown|other", unique_phrogs, value = TRUE, ignore.case = TRUE)
  
  # Create colors
  phrog_colors <- rainbow(n_phrogs, alpha = 0.8)
  names(phrog_colors) <- unique_phrogs
  
  # Set gray colors for unknown/other
  for(gray_cat in gray_phrog_categories) {
    if(gray_cat %in% names(phrog_colors)) {
      phrog_colors[gray_cat] <- "grey60"
    }
  }
  
  # Reorder data so gray points are plotted first
  data_ordered <- data[order(data$phrog_category %in% gray_phrog_categories, decreasing = TRUE), ]
  
  # Create volcano plot
  p <- ggplot(data_ordered, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = phrog_category,
                   size = ifelse(phrog_category %in% gray_phrog_categories, 1.2, 1.8),
                   alpha = ifelse(phrog_category %in% gray_phrog_categories, 0.4, 0.9))) +
    scale_color_manual(values = phrog_colors) +
    scale_size_identity() +
    scale_alpha_identity() +
    labs(
      title = "Volcano Plot: Significant ORFs by Phrog Categories (US Cohort)", 
      subtitle = paste("Showing", nrow(data), "ORFs with adj.p < 0.01"),
      x = "Effect Size (logFC)\nPositive = Higher in CONTROL, Negative = Higher in CELIAC",
      y = "-log10(adjusted P-value)",
      color = "Phrog Category"
    ) +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  return(p)
}

phrog_plot <- create_phrog_volcano(merged_data)
ggsave("volcano_plot_by_phrog_US.png", phrog_plot, width = 14, height = 10, dpi = 300, bg = "white")
cat("Phrog volcano plot saved\n")

# ============================================================================
# PART 5: VOLCANO PLOT BY PRODUCT CATEGORIES
# ============================================================================

cat("\n5. Creating volcano plot by product categories...\n")

create_product_volcano <- function(data) {
  
  # Create product categories
  data$product_category <- ifelse(is.na(data$product) | data$product == "", 
                                  "unknown product", 
                                  as.character(data$product))
  
  # Get top product categories
  product_counts <- table(data$product_category)
  top_products <- names(sort(product_counts, decreasing = TRUE))[1:20]
  
  # Group rare products
  data$product_category[!data$product_category %in% top_products] <- "other products"
  
  # Create color palette
  unique_products <- unique(data$product_category)
  n_products <- length(unique_products)
  
  # Identify gray categories
  gray_product_categories <- grep("unknown|other", unique_products, value = TRUE, ignore.case = TRUE)
  
  # Create colors
  product_colors <- rainbow(n_products, alpha = 0.8)
  names(product_colors) <- unique_products
  
  # Set gray colors for unknown/other
  for(gray_cat in gray_product_categories) {
    if(gray_cat %in% names(product_colors)) {
      product_colors[gray_cat] <- "grey60"
    }
  }
  
  # Reorder data so gray points are plotted first
  data_ordered <- data[order(data$product_category %in% gray_product_categories, decreasing = TRUE), ]
  
  # Create volcano plot
  p <- ggplot(data_ordered, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = product_category,
                   size = ifelse(product_category %in% gray_product_categories, 1.2, 1.8),
                   alpha = ifelse(product_category %in% gray_product_categories, 0.4, 0.9))) +
    scale_color_manual(values = product_colors) +
    scale_size_identity() +
    scale_alpha_identity() +
    labs(
      title = "Volcano Plot: Significant ORFs by Product Categories (US Cohort)",
      subtitle = paste("Showing", nrow(data), "ORFs with adj.p < 0.01"),
      x = "Effect Size (logFC)\nPositive = Higher in CONTROL, Negative = Higher in CELIAC",
      y = "-log10(adjusted P-value)",
      color = "Product Category"
    ) +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  return(p)
}

product_plot <- create_product_volcano(merged_data)
ggsave("volcano_plot_by_product_US.png", product_plot, width = 14, height = 10, dpi = 300, bg = "white")
cat("Product volcano plot saved\n")

# ============================================================================
# PART 6: SAVE SUMMARY DATA
# ============================================================================

cat("\n6. Saving summary data...\n")

# Save merged data with functional annotations
write.csv(merged_data, "functional_annotations_summary_US.csv", row.names = FALSE)

# Create summary statistics
annotation_summary <- data.frame(
  category = c("Total significant ORFs", "With phrog annotation", "With function annotation", "With product annotation"),
  count = c(nrow(merged_data),
            sum(!is.na(merged_data$phrog)),
            sum(!is.na(merged_data$func_category)),
            sum(!is.na(merged_data$product))),
  percentage = c(100,
                 round(sum(!is.na(merged_data$phrog))/nrow(merged_data)*100, 1),
                 round(sum(!is.na(merged_data$func_category))/nrow(merged_data)*100, 1),
                 round(sum(!is.na(merged_data$product))/nrow(merged_data)*100, 1))
)

write.csv(annotation_summary, "annotation_coverage_summary_US.csv", row.names = FALSE)

cat("Summary data saved\n")

# ============================================================================
# PART 7: FINAL SUMMARY
# ============================================================================

cat("\n=== VOLCANO PLOTS ANALYSIS COMPLETED (US COHORT) ===\n")
cat("Generated files:\n")
cat("- volcano_plot_by_function_US.png: Volcano plot colored by functional categories\n")
cat("- volcano_plot_by_phrog_US.png: Volcano plot colored by phrog categories\n")
cat("- volcano_plot_by_product_US.png: Volcano plot colored by product categories\n")
cat("- functional_annotations_summary_US.csv: Merged data with annotations\n")
cat("- annotation_coverage_summary_US.csv: Annotation coverage statistics\n")

cat("\nAnnotation coverage (US Cohort):\n")
print(annotation_summary)

cat("\nAnalysis completed at:", format(Sys.time()), "\n")