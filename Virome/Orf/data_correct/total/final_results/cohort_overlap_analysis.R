#!/usr/bin/env Rscript
# Cohort Overlap Analysis: Significant ORFs across Total, US, and Italy Cohorts
# Author: Claude AI
# Date: 2025-08-06

library(dplyr)
library(ggplot2)
library(VennDiagram)
library(gridExtra)

cat("=== COHORT OVERLAP ANALYSIS ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

cat("1. Loading limma results from all three cohorts...\n")

# Load Total cohort results (adj.p < 0.01)
total_results <- read.csv("total_limma_model_Sig_res_final.csv", row.names = 1)
cat("Total cohort: Loaded", nrow(total_results), "significant ORFs (adj.p < 0.01)\n")

# Load US cohort results (adj.p < 0.05)  
us_results <- read.csv("../../US/final_results/US_limma_model_Sig_res_final.csv", row.names = 1)
cat("US cohort: Loaded", nrow(us_results), "significant ORFs (adj.p < 0.05)\n")

# Load Italy cohort results (adj.p < 0.05)
italy_results <- read.csv("../../Italy/final_results/Italy_limma_model_Sig_res_final.csv", row.names = 1)
cat("Italy cohort: Loaded", nrow(italy_results), "significant ORFs (adj.p < 0.05)\n")

# ============================================================================
# PART 2: EXTRACT ORF LISTS
# ============================================================================

cat("\n2. Extracting ORF identifiers from each cohort...\n")

# Extract ORF names
total_orfs <- rownames(total_results)
us_orfs <- rownames(us_results) 
italy_orfs <- rownames(italy_results)

cat("ORF lists extracted:\n")
cat("- Total cohort:", length(total_orfs), "ORFs\n")
cat("- US cohort:", length(us_orfs), "ORFs\n") 
cat("- Italy cohort:", length(italy_orfs), "ORFs\n")

# ============================================================================
# PART 3: OVERLAP ANALYSIS
# ============================================================================

cat("\n3. Performing overlap analysis...\n")

# Two-way overlaps
total_us_overlap <- intersect(total_orfs, us_orfs)
total_italy_overlap <- intersect(total_orfs, italy_orfs)
us_italy_overlap <- intersect(us_orfs, italy_orfs)

# Three-way overlap
all_three_overlap <- intersect(intersect(total_orfs, us_orfs), italy_orfs)

cat("\nOverlap Results:\n")
cat("Total vs US:", length(total_us_overlap), "ORFs\n")
cat("Total vs Italy:", length(total_italy_overlap), "ORFs\n")
cat("US vs Italy:", length(us_italy_overlap), "ORFs\n")
cat("All three cohorts:", length(all_three_overlap), "ORFs\n")

# Calculate percentages
cat("\nOverlap Percentages:\n")
cat("Total-US overlap as % of Total:", round(100 * length(total_us_overlap) / length(total_orfs), 1), "%\n")
cat("Total-US overlap as % of US:", round(100 * length(total_us_overlap) / length(us_orfs), 1), "%\n")
cat("Total-Italy overlap as % of Total:", round(100 * length(total_italy_overlap) / length(total_orfs), 1), "%\n")
cat("Total-Italy overlap as % of Italy:", round(100 * length(total_italy_overlap) / length(italy_orfs), 1), "%\n")

# ============================================================================
# PART 4: DETAILED ANALYSIS OF TOTAL-US OVERLAP
# ============================================================================

cat("\n4. Detailed analysis of Total-US overlap...\n")

# Get detailed information for overlapping ORFs
total_us_details <- data.frame(
  orf_id = total_us_overlap,
  total_logFC = total_results[total_us_overlap, "logFC"],
  total_adj_p = total_results[total_us_overlap, "adj.P.Val"],
  us_logFC = us_results[total_us_overlap, "logFC"],
  us_adj_p = us_results[total_us_overlap, "adj.P.Val"],
  stringsAsFactors = FALSE
)

# Calculate effect size correlation
logfc_correlation <- cor(total_us_details$total_logFC, total_us_details$us_logFC)
cat("LogFC correlation between Total and US cohorts:", round(logfc_correlation, 3), "\n")

# Direction consistency
same_direction <- sum(sign(total_us_details$total_logFC) == sign(total_us_details$us_logFC))
cat("ORFs with same direction (sign) in both cohorts:", same_direction, "out of", length(total_us_overlap), 
    "(", round(100 * same_direction / length(total_us_overlap), 1), "%)\n")

# ============================================================================
# PART 5: CREATE VISUALIZATIONS
# ============================================================================

cat("\n5. Creating visualizations...\n")

# Create summary data frame for plotting
overlap_summary <- data.frame(
  Comparison = c("Total vs US", "Total vs Italy", "US vs Italy", "All Three"),
  Count = c(length(total_us_overlap), length(total_italy_overlap), 
            length(us_italy_overlap), length(all_three_overlap)),
  Percentage_of_smaller = c(
    round(100 * length(total_us_overlap) / min(length(total_orfs), length(us_orfs)), 1),
    round(100 * length(total_italy_overlap) / min(length(total_orfs), length(italy_orfs)), 1),
    round(100 * length(us_italy_overlap) / min(length(us_orfs), length(italy_orfs)), 1),
    round(100 * length(all_three_overlap) / min(length(total_orfs), length(us_orfs), length(italy_orfs)), 1)
  )
)

# Bar plot of overlaps
p1 <- ggplot(overlap_summary, aes(x = Comparison, y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = Count), vjust = -0.3, size = 4, fontface = "bold") +
  labs(title = "Significant ORF Overlaps Between Cohorts",
       subtitle = "Number of overlapping significant ORFs",
       x = "Cohort Comparison",
       y = "Number of Overlapping ORFs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)) +
  scale_fill_brewer(palette = "Set2")

# Effect size correlation plot (Total vs US)
p2 <- ggplot(total_us_details, aes(x = total_logFC, y = us_logFC)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "Effect Size Correlation: Total vs US Cohort",
       subtitle = paste("Correlation:", round(logfc_correlation, 3), "| n =", length(total_us_overlap), "overlapping ORFs"),
       x = "Total Cohort logFC",
       y = "US Cohort logFC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

# Combined plot
combined_plot <- grid.arrange(p1, p2, nrow = 1, ncol = 2)

# Save plots
ggsave("cohort_overlap_analysis.png", combined_plot, width = 16, height = 8, dpi = 300, bg = "white")
cat("Combined overlap analysis plot saved as 'cohort_overlap_analysis.png'\n")

# ============================================================================
# PART 6: SAVE RESULTS
# ============================================================================

cat("\n6. Saving detailed results...\n")

# Save overlap lists
write.csv(total_us_details, "total_us_overlapping_orfs.csv", row.names = FALSE)
write.csv(data.frame(orf_id = total_italy_overlap), "total_italy_overlapping_orfs.csv", row.names = FALSE)
write.csv(data.frame(orf_id = all_three_overlap), "all_cohorts_overlapping_orfs.csv", row.names = FALSE)

# Save summary statistics
summary_stats <- data.frame(
  Cohort = c("Total", "US", "Italy"),
  Significant_ORFs = c(length(total_orfs), length(us_orfs), length(italy_orfs)),
  Significance_Threshold = c("adj.p < 0.01", "adj.p < 0.05", "adj.p < 0.05"),
  Overlap_with_Total = c(length(total_orfs), length(total_us_overlap), length(total_italy_overlap)),
  Overlap_Percentage = c(100, 
                        round(100 * length(total_us_overlap) / length(us_orfs), 1),
                        round(100 * length(total_italy_overlap) / length(italy_orfs), 1))
)

write.csv(summary_stats, "cohort_overlap_summary.csv", row.names = FALSE)
write.csv(overlap_summary, "overlap_comparison_details.csv", row.names = FALSE)

cat("Detailed results saved:\n")
cat("- total_us_overlapping_orfs.csv: Detailed overlap between Total and US\n")
cat("- total_italy_overlapping_orfs.csv: ORFs overlapping between Total and Italy\n") 
cat("- all_cohorts_overlapping_orfs.csv: ORFs significant in all three cohorts\n")
cat("- cohort_overlap_summary.csv: Summary statistics\n")
cat("- overlap_comparison_details.csv: Comparison details\n")

# ============================================================================
# PART 7: FINAL SUMMARY
# ============================================================================

cat("\n=== COHORT OVERLAP ANALYSIS COMPLETED ===\n")
cat("Key Findings:\n")
cat("1. Total-US Overlap:", length(total_us_overlap), "ORFs (", 
    round(100 * length(total_us_overlap) / length(total_orfs), 1), "% of Total,", 
    round(100 * length(total_us_overlap) / length(us_orfs), 1), "% of US)\n")
cat("2. Effect Size Correlation (Total-US):", round(logfc_correlation, 3), "\n")
cat("3. Direction Consistency (Total-US):", round(100 * same_direction / length(total_us_overlap), 1), "%\n")
cat("4. Total-Italy Overlap:", length(total_italy_overlap), "ORFs (", 
    round(100 * length(total_italy_overlap) / length(italy_orfs), 1), "% of Italy)\n")
cat("5. All-cohort Core ORFs:", length(all_three_overlap), "ORFs\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")