#!/usr/bin/env Rscript
# Create Volcano Plots with Functional Annotations
# Colors significant ORFs by phrog, function, and product categories
# Author: Claude AI
# Date: 2025-07-25

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data_correct/total/final_results")

cat("=== VOLCANO PLOT ANALYSIS WITH FUNCTIONAL ANNOTATIONS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading data...\n")

# Load limma results
limma_results <- read.csv("total_limma_model_res.csv", row.names = 1)
limma_results$orf_id <- rownames(limma_results)

# Load phold functional predictions
phold_data <- read.delim("phold_per_cds_predictions.tsv", stringsAsFactors = FALSE)
# Rename function column to avoid R reserved word conflict (note the dot in column name)
names(phold_data)[names(phold_data) == "function."] <- "func_category"

cat("Data loaded:\n")
cat("  Limma results:", nrow(limma_results), "ORFs\n")
cat("  Phold predictions:", nrow(phold_data), "ORFs\n")

# ============================================================================
# PART 2: MERGE DATA AND PREPARE FOR PLOTTING
# ============================================================================

cat("\n2. Merging data and preparing for plotting...\n")

# Merge limma results with functional annotations
merged_data <- merge(limma_results, phold_data, 
                     by.x = "orf_id", by.y = "cds_id", 
                     all.x = TRUE)

cat("Merged data:", nrow(merged_data), "ORFs\n")
cat("ORFs with functional annotations:", sum(!is.na(merged_data$func_category)), "\n")

# Clean up function categories - handle missing values and empty strings
merged_data$func_category <- ifelse(is.na(merged_data$func_category) | merged_data$func_category == "", 
                                   "no annotation", merged_data$func_category)
merged_data$phrog <- ifelse(is.na(merged_data$phrog) | merged_data$phrog == "", 
                           "No_PHROG", merged_data$phrog)
merged_data$product <- ifelse(is.na(merged_data$product) | merged_data$product == "", 
                             "no product annotation", merged_data$product)

# Add significance categories (using original logFC for thresholds)
merged_data$significance <- case_when(
  merged_data$adj.P.Val < 0.01 & abs(merged_data$logFC) >= 0.5 ~ "Significant (adj.p < 0.01, |logFC| >= 0.5)",
  merged_data$adj.P.Val < 0.01 ~ "Significant (adj.p < 0.01)",
  merged_data$P.Value < 0.05 ~ "Nominally significant (p < 0.05)",
  TRUE ~ "Not significant"
)

# Create -log10(P.Value) for plotting
merged_data$neg_log10_pval <- -log10(merged_data$P.Value)
merged_data$neg_log10_adj_pval <- -log10(merged_data$adj.P.Val)

# Flip logFC so CONTROL is on left (negative) and CELIAC is on right (positive)
merged_data$logFC_flipped <- -merged_data$logFC

# Cap extreme values for better visualization
merged_data$neg_log10_pval[merged_data$neg_log10_pval > 20] <- 20
merged_data$neg_log10_adj_pval[merged_data$neg_log10_adj_pval > 20] <- 20

cat("Significance summary:\n")
cat(table(merged_data$significance), "\n")

# ============================================================================
# PART 3: FUNCTION-BASED VOLCANO PLOT
# ============================================================================

cat("\n3. Creating function-based volcano plot...\n")

# Filter to only show significant ORFs (adj.p < 0.01)
sig_data <- merged_data[merged_data$adj.P.Val < 0.01, ]
cat("Plotting", nrow(sig_data), "significant ORFs (adj.p < 0.01)\n")

# Get top functions for coloring (exclude unknown function to focus on known functions)
function_counts <- table(sig_data$func_category)
top_functions <- names(sort(function_counts, decreasing = TRUE))[1:8]

# Create function categories for coloring
sig_data$function_category <- ifelse(sig_data$func_category %in% top_functions, 
                                    sig_data$func_category, "Other/Unknown")

# Define colors for functions
function_colors <- c(
  "head and packaging" = "#E31A1C",
  "DNA, RNA and nucleotide metabolism" = "#1F78B4", 
  "tail" = "#33A02C",
  "transcription regulation" = "#FF7F00",
  "lysis" = "#6A3D9A",
  "connector" = "#A6CEE3",
  "integration and excision" = "#FB9A99",
  "other" = "#FDBF6F",
  "Other/Unknown" = "grey70",
  "unknown function" = "grey60"
)

# Create volcano plot colored by function (only significant ORFs)
# Identify all gray-colored categories
gray_function_categories <- names(function_colors)[grepl("grey", function_colors)]
# Reorder data so gray function points are plotted first (will appear behind)
# Use decreasing=TRUE to put gray points (TRUE) first, colored points (FALSE) last
sig_data_ordered <- sig_data[order(sig_data$function_category %in% gray_function_categories, decreasing = TRUE), ]

p1 <- ggplot(sig_data_ordered, aes(x = logFC_flipped, y = neg_log10_adj_pval)) +
  geom_point(aes(color = function_category, 
                 size = ifelse(function_category %in% gray_function_categories, 1.2, 1.8),
                 alpha = ifelse(function_category %in% gray_function_categories, 0.4, 0.9))) +
  scale_color_manual(values = function_colors, name = "Function Category") +
  scale_size_identity() +
  scale_alpha_identity() +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.7) +
  labs(
    title = paste("Volcano Plot:", nrow(sig_data), "Significant Viral ORFs (adj.p < 0.01)"),
    subtitle = "Colored by Function Category\nLeft = Higher in CONTROL, Right = Higher in CELIAC",
    x = "log2 Fold Change (CONTROL ← → CELIAC)",
    y = "-log10(adjusted P-value)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white")
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 2))
  )

ggsave("volcano_plot_by_function.png", p1, width = 14, height = 10, dpi = 300)

# ============================================================================
# PART 4: PHROG-BASED VOLCANO PLOT  
# ============================================================================

cat("\n4. Creating PHROG-based volcano plot...\n")

# For PHROG, we'll group by broad categories since there are many individual PHROGs
# Create broader PHROG categories based on the function mapping (using sig_data)
sig_data$phrog_category <- case_when(
  sig_data$phrog == "No_PHROG" ~ "No PHROG annotation",
  sig_data$func_category == "head and packaging" ~ "Head and packaging PHROGs",
  sig_data$func_category == "DNA, RNA and nucleotide metabolism" ~ "DNA/RNA metabolism PHROGs", 
  sig_data$func_category == "tail" ~ "Tail PHROGs",
  sig_data$func_category == "transcription regulation" ~ "Transcription PHROGs",
  sig_data$func_category == "lysis" ~ "Lysis PHROGs",
  sig_data$func_category == "connector" ~ "Connector PHROGs",
  sig_data$func_category == "integration and excision" ~ "Integration/excision PHROGs",
  sig_data$func_category == "other" ~ "Other PHROGs",
  TRUE ~ "Unknown function PHROGs"
)

# Define colors for PHROG categories
phrog_colors <- c(
  "Head and packaging PHROGs" = "#E31A1C",
  "DNA/RNA metabolism PHROGs" = "#1F78B4",
  "Tail PHROGs" = "#33A02C", 
  "Transcription PHROGs" = "#FF7F00",
  "Lysis PHROGs" = "#6A3D9A",
  "Connector PHROGs" = "#A6CEE3",
  "Integration/excision PHROGs" = "#FB9A99",
  "Other PHROGs" = "#FDBF6F",
  "Unknown function PHROGs" = "grey60",
  "No PHROG annotation" = "grey80"
)

# Identify all gray-colored PHROG categories
gray_phrog_categories <- names(phrog_colors)[grepl("grey", phrog_colors)]
# Reorder data so gray PHROG points are plotted first (will appear behind)
# Use decreasing=TRUE to put gray points (TRUE) first, colored points (FALSE) last
sig_data_ordered_phrog <- sig_data[order(sig_data$phrog_category %in% gray_phrog_categories, decreasing = TRUE), ]

p2 <- ggplot(sig_data_ordered_phrog, aes(x = logFC_flipped, y = neg_log10_adj_pval)) +
  geom_point(aes(color = phrog_category, 
                 size = ifelse(phrog_category %in% gray_phrog_categories, 1.2, 1.8),
                 alpha = ifelse(phrog_category %in% gray_phrog_categories, 0.4, 0.9))) +
  scale_color_manual(values = phrog_colors, name = "PHROG Category") +
  scale_size_identity() +
  scale_alpha_identity() +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.7) +
  labs(
    title = paste("Volcano Plot:", nrow(sig_data), "Significant Viral ORFs (adj.p < 0.01)"),
    subtitle = "Colored by PHROG Category\nLeft = Higher in CONTROL, Right = Higher in CELIAC",
    x = "log2 Fold Change (CONTROL ← → CELIAC)",
    y = "-log10(adjusted P-value)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white")
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 2))
  )

ggsave("volcano_plot_by_phrog.png", p2, width = 14, height = 10, dpi = 300)

# ============================================================================
# PART 5: PRODUCT-BASED VOLCANO PLOT
# ============================================================================

cat("\n5. Creating product-based volcano plot...\n")

# Get top products for significant ORFs (focusing on informative product names)
product_counts <- table(sig_data$product)
# Filter out generic terms
informative_products <- names(product_counts)[!grepl("hypothetical protein|no product annotation", names(product_counts))]
top_products <- names(sort(product_counts[informative_products], decreasing = TRUE))[1:10]

# Create product categories
sig_data$product_category <- case_when(
  sig_data$product %in% top_products ~ sig_data$product,
  grepl("hypothetical protein|no product annotation", sig_data$product) ~ "Hypothetical/Unknown",
  TRUE ~ "Other annotated products"
)

# Define colors for products (using a diverse palette)
n_product_cats <- length(unique(sig_data$product_category))
product_colors <- c(
  RColorBrewer::brewer.pal(min(n_product_cats-2, 11), "Spectral"),
  "Hypothetical/Unknown" = "grey70",
  "Other annotated products" = "grey50"
)
names(product_colors)[1:(n_product_cats-2)] <- setdiff(unique(sig_data$product_category), 
                                                       c("Hypothetical/Unknown", "Other annotated products"))

# Identify all gray-colored product categories
gray_product_categories <- names(product_colors)[grepl("grey", product_colors)]
# Reorder data so gray product points are plotted first (will appear behind)
# Use decreasing=TRUE to put gray points (TRUE) first, colored points (FALSE) last
sig_data_ordered_product <- sig_data[order(sig_data$product_category %in% gray_product_categories, decreasing = TRUE), ]

p3 <- ggplot(sig_data_ordered_product, aes(x = logFC_flipped, y = neg_log10_adj_pval)) +
  geom_point(aes(color = product_category, 
                 size = ifelse(product_category %in% gray_product_categories, 1.2, 1.8),
                 alpha = ifelse(product_category %in% gray_product_categories, 0.4, 0.9))) +
  scale_color_manual(values = product_colors, name = "Product Category") +
  scale_size_identity() +
  scale_alpha_identity() +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.7) +
  labs(
    title = paste("Volcano Plot:", nrow(sig_data), "Significant Viral ORFs (adj.p < 0.01)"),
    subtitle = "Colored by Product Category\nLeft = Higher in CONTROL, Right = Higher in CELIAC",
    x = "log2 Fold Change (CONTROL ← → CELIAC)",
    y = "-log10(adjusted P-value)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white")
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 1)
  )

ggsave("volcano_plot_by_product.png", p3, width = 16, height = 10, dpi = 300)

# ============================================================================
# PART 6: SUMMARY STATISTICS
# ============================================================================

cat("\n6. Generating summary statistics...\n")

# Function distribution among significant ORFs
sig_function_summary <- sig_data %>%
  group_by(func_category) %>%
  summarise(
    count = n(),
    mean_logFC = mean(logFC),
    mean_adj_pval = mean(adj.P.Val),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

write.csv(sig_function_summary, "significant_orfs_function_summary.csv", row.names = FALSE)

# PHROG category distribution
sig_phrog_summary <- sig_data %>%
  group_by(phrog_category) %>%
  summarise(
    count = n(),
    mean_logFC = mean(logFC),
    mean_adj_pval = mean(adj.P.Val),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

write.csv(sig_phrog_summary, "significant_orfs_phrog_summary.csv", row.names = FALSE)

# Product distribution among significant ORFs (top 20)
sig_product_summary <- sig_data %>%
  group_by(product) %>%
  summarise(
    count = n(),
    mean_logFC = mean(logFC),
    mean_adj_pval = mean(adj.P.Val),
    .groups = "drop"
  ) %>%
  arrange(desc(count)) %>%
  head(20)

write.csv(sig_product_summary, "significant_orfs_product_summary.csv", row.names = FALSE)

# Overall summary
overall_summary <- data.frame(
  metric = c("Total ORFs", "ORFs with functional annotation", "Significant ORFs (adj.p<0.01) - plotted",
             "Highly significant ORFs (adj.p<0.01 & |logFC|>=0.5)", "Unknown function (significant)",
             "Known function (significant)"),
  count = c(
    nrow(merged_data),
    sum(merged_data$func_category != "no annotation"),
    nrow(sig_data),
    sum(sig_data$adj.P.Val < 0.01 & abs(sig_data$logFC) >= 0.5),
    sum(sig_data$func_category == "unknown function"),
    sum(sig_data$func_category != "unknown function")
  )
)

write.csv(overall_summary, "volcano_plot_summary.csv", row.names = FALSE)

cat("\n=== VOLCANO PLOT ANALYSIS COMPLETED ===\n")
cat("Generated files:\n")
cat("- volcano_plot_by_function.png: Volcano plot colored by function categories\n")
cat("- volcano_plot_by_phrog.png: Volcano plot colored by PHROG categories\n") 
cat("- volcano_plot_by_product.png: Volcano plot colored by product annotations\n")
cat("- significant_orfs_function_summary.csv: Function distribution among significant ORFs\n")
cat("- significant_orfs_phrog_summary.csv: PHROG distribution among significant ORFs\n")
cat("- significant_orfs_product_summary.csv: Product distribution among significant ORFs\n")
cat("- volcano_plot_summary.csv: Overall summary statistics\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")