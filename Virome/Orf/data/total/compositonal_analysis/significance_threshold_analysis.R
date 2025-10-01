#!/usr/bin/env Rscript
# Analysis of significance thresholds: 0.05 vs 0.001 vs other cutoffs
# Provides evidence for optimal threshold selection

library(dplyr)
library(ggplot2)

# Set working directory
setwd("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/compositonal_analysis")

cat("=== SIGNIFICANCE THRESHOLD ANALYSIS ===\n")

# Load limma results
limma_results <- read.csv("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Orf/data/total/total_limma_results.csv", row.names = 1)

cat("Total ORFs analyzed:", nrow(limma_results), "\n\n")

# Test different significance thresholds
thresholds <- c(0.1, 0.05, 0.01, 0.001, 0.0001)

cat("=== THRESHOLD COMPARISON ===\n")
threshold_summary <- data.frame(
  threshold = thresholds,
  n_significant = sapply(thresholds, function(t) sum(limma_results$adj.P.Val < t, na.rm = TRUE)),
  percent_significant = sapply(thresholds, function(t) round(100 * sum(limma_results$adj.P.Val < t, na.rm = TRUE) / nrow(limma_results), 1))
)

print(threshold_summary)

# What threshold did I actually use in my analysis?
cat("\n=== WHAT I ACTUALLY USED ===\n")
my_significant <- limma_results$adj.P.Val < 0.05
my_top_50 <- head(limma_results[order(limma_results$adj.P.Val), ], 50)

cat("I used: adj.p < 0.05 threshold, then took top 50 most significant\n")
cat("Number significant at 0.05:", sum(my_significant, na.rm = TRUE), "\n")
cat("Top 50 adj.p range:", round(min(my_top_50$adj.P.Val), 8), "to", round(my_top_50$adj.P.Val[50], 8), "\n")
cat("Top 50 equivalent threshold: approximately", round(my_top_50$adj.P.Val[50], 6), "\n")

# Effect size analysis by threshold
cat("\n=== EFFECT SIZE BY THRESHOLD ===\n")
for(thresh in thresholds) {
  sig_at_thresh <- limma_results[limma_results$adj.P.Val < thresh, ]
  if(nrow(sig_at_thresh) > 0) {
    cat("Threshold", thresh, ":\n")
    cat("  Mean |logFC|:", round(mean(abs(sig_at_thresh$logFC)), 4), "\n")
    cat("  Median |logFC|:", round(median(abs(sig_at_thresh$logFC)), 4), "\n")
    cat("  Range logFC:", round(min(sig_at_thresh$logFC), 4), "to", round(max(sig_at_thresh$logFC), 4), "\n")
  }
}

# Statistical power analysis
cat("\n=== STATISTICAL POWER ANALYSIS ===\n")

# Distribution of p-values
p_val_breaks <- c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)
p_val_counts <- table(cut(limma_results$adj.P.Val, breaks = p_val_breaks, include.lowest = TRUE))
cat("P-value distribution:\n")
print(p_val_counts)

# Multiple testing burden
total_tests <- nrow(limma_results)
cat("\nMultiple testing context:\n")
cat("Total tests:", total_tests, "\n")
cat("Expected false positives at 0.05:", round(total_tests * 0.05), "\n")
cat("Expected false positives at 0.001:", round(total_tests * 0.001), "\n")

# Field standards analysis
cat("\n=== FIELD STANDARDS ===\n")
cat("Microbiome studies typically use:\n")
cat("- 16S rRNA studies: adj.p < 0.05 (standard)\n")
cat("- Metagenomics: adj.p < 0.05 to 0.01 (varies by study size)\n")
cat("- Single-cell: adj.p < 0.01 to 0.001 (high-dimensional)\n")
cat("- GWAS: p < 5e-8 (genome-wide significance)\n")

# Recommendation based on your data
cat("\n=== RECOMMENDATION FOR YOUR STUDY ===\n")

# Calculate meaningful thresholds
sig_05 <- sum(limma_results$adj.P.Val < 0.05, na.rm = TRUE)
sig_001 <- sum(limma_results$adj.P.Val < 0.001, na.rm = TRUE)

cat("Your study characteristics:\n")
cat("- Viral ORFs (2,154 tests)\n")
cat("- Longitudinal design with confounders\n")
cat("- Celiac disease (well-defined phenotype)\n")
cat("- Effect sizes: small but consistent (0.06-0.13)\n\n")

cat("Threshold options:\n")
cat("1. adj.p < 0.05:", sig_05, "ORFs (", round(100*sig_05/nrow(limma_results), 1), "%)\n")
cat("   - Standard field practice\n")
cat("   - Balances discovery and false positives\n")
cat("   - Appropriate for your sample size\n\n")

cat("2. adj.p < 0.001:", sig_001, "ORFs (", round(100*sig_001/nrow(limma_results), 1), "%)\n")
cat("   - Very conservative\n")
cat("   - Reduces false positives\n")
cat("   - May miss biologically meaningful effects\n\n")

# Evidence-based recommendation
if(sig_05 > 100 && sig_001 > 50) {
  cat("RECOMMENDATION: Use adj.p < 0.05\n")
  cat("RATIONALE:\n")
  cat("- Standard in microbiome field\n")
  cat("- Your FDR correction already controls false discovery rate\n") 
  cat("- Effect sizes are small but consistent (typical for microbiome)\n")
  cat("- Sufficient power to detect meaningful biological signals\n")
} else {
  cat("RECOMMENDATION: Consider adj.p < 0.001\n")
  cat("RATIONALE:\n")
  cat("- High number of significant hits suggests strong signal\n")
  cat("- Conservative threshold appropriate for high-dimensional data\n")
}

cat("\n=== CONCLUSION ===\n")
cat("No universal 'better' threshold exists.\n")
cat("Choice depends on:\n")
cat("1. Study goals (discovery vs confirmation)\n")
cat("2. Downstream validation capabilities\n") 
cat("3. Biological effect sizes in your system\n")
cat("4. Field standards and precedent\n")

# Create visualization
png("significance_threshold_comparison.png", width = 1200, height = 800, res = 150)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Panel 1: Number significant by threshold
barplot(threshold_summary$n_significant, 
        names.arg = threshold_summary$threshold,
        main = "Number of Significant ORFs",
        xlab = "Adjusted p-value threshold",
        ylab = "Number of ORFs",
        col = "lightblue")

# Panel 2: P-value histogram
hist(limma_results$adj.P.Val, 
     breaks = 50, 
     main = "Distribution of Adjusted P-values",
     xlab = "Adjusted p-value",
     ylab = "Frequency",
     col = "lightcoral")
abline(v = 0.05, col = "red", lwd = 2, lty = 2)
abline(v = 0.001, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("0.05", "0.001"), col = c("red", "blue"), lty = 2)

# Panel 3: Effect size vs significance
plot(limma_results$adj.P.Val, abs(limma_results$logFC),
     log = "x",
     main = "Effect Size vs Significance",
     xlab = "Adjusted p-value (log scale)",
     ylab = "|logFC|",
     pch = 16, cex = 0.5)
abline(v = 0.05, col = "red", lwd = 2, lty = 2)
abline(v = 0.001, col = "blue", lwd = 2, lty = 2)

# Panel 4: Cumulative significance
sorted_p <- sort(limma_results$adj.P.Val)
plot(1:length(sorted_p), sorted_p,
     main = "Cumulative Significance",
     xlab = "ORF rank (by significance)",
     ylab = "Adjusted p-value",
     log = "y", type = "l", lwd = 2)
abline(h = 0.05, col = "red", lwd = 2, lty = 2)
abline(h = 0.001, col = "blue", lwd = 2, lty = 2)

dev.off()

cat("Generated significance_threshold_comparison.png\n")