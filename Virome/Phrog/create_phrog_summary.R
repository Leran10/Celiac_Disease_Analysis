# Create PHROG Analysis Summary and Figures
library(dplyr)
library(ggplot2)
library(readr)

# Read the completed results
total_logit <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_logit_results.csv")
total_cloglog <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_cloglog_results.csv")

# Summary statistics
cat("PHROG Analysis Summary\n")
cat("=====================\n")
cat("Total cohort - glmmTMB logit results:", nrow(total_logit), "rows\n")
cat("Total cohort - glmmTMB cloglog results:", nrow(total_cloglog), "rows\n")
cat("Unique PHROGs in logit results:", length(unique(total_logit$phrog_id)), "\n")
cat("Unique PHROGs in cloglog results:", length(unique(total_cloglog$phrog_id)), "\n")

# Count significant results (p < 0.05)
logit_sig <- total_logit %>% filter(p.value < 0.05 & !is.na(p.value))
cloglog_sig <- total_cloglog %>% filter(p.value < 0.05 & !is.na(p.value))

cat("Significant results (p < 0.05):\n")
cat("- Logit model:", nrow(logit_sig), "coefficients\n")
cat("- Cloglog model:", nrow(cloglog_sig), "coefficients\n")

# Focus on Dx.Status effects
dx_effects_logit <- total_logit %>% 
  filter(grepl("Dx.StatusCELIAC", term)) %>%
  filter(!is.na(p.value))

dx_effects_cloglog <- total_cloglog %>% 
  filter(grepl("Dx.StatusCELIAC", term)) %>%
  filter(!is.na(p.value))

cat("Dx.Status effects:\n")
cat("- Logit model:", nrow(dx_effects_logit), "PHROGs with valid results\n")
cat("- Cloglog model:", nrow(dx_effects_cloglog), "PHROGs with valid results\n")

# Count significant Dx.Status effects
dx_sig_logit <- dx_effects_logit %>% filter(p.value < 0.05)
dx_sig_cloglog <- dx_effects_cloglog %>% filter(p.value < 0.05)

cat("Significant Dx.Status effects (p < 0.05):\n")
cat("- Logit model:", nrow(dx_sig_logit), "PHROGs\n")
cat("- Cloglog model:", nrow(dx_sig_cloglog), "PHROGs\n")

# Create effect size plot for Dx.Status main effects
if (nrow(dx_effects_logit) > 0) {
  p1 <- ggplot(dx_effects_logit, aes(x = estimate)) +
    geom_histogram(bins = 50, fill = "lightblue", alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    labs(title = "PHROG PA Analysis - Dx.Status Effect Sizes (Logit Model)",
         subtitle = paste("Total cohort -", nrow(dx_effects_logit), "PHROGs with valid results,", nrow(dx_sig_logit), "significant"),
         x = "Log Odds Ratio (Dx.StatusCELIAC)",
         y = "Count") +
    theme_minimal()
  
  ggsave("../Orf_Contig_Phrog_compositional/report_figures/phrog_PA_dx_effects_logit.pdf", p1, width = 10, height = 6)
}

# Create p-value distribution plot
if (nrow(dx_effects_logit) > 0) {
  p2 <- ggplot(dx_effects_logit, aes(x = p.value)) +
    geom_histogram(bins = 50, fill = "lightgreen", alpha = 0.7) +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
    labs(title = "PHROG PA Analysis - P-value Distribution (Dx.Status)",
         subtitle = paste("Total cohort - Logit model -", nrow(dx_sig_logit), "significant out of", nrow(dx_effects_logit)),
         x = "P-value",
         y = "Count") +
    theme_minimal()
  
  ggsave("../Orf_Contig_Phrog_compositional/report_figures/phrog_PA_pvalue_dist_logit.pdf", p2, width = 10, height = 6)
}

# Volcano plot
if (nrow(dx_effects_logit) > 0) {
  dx_effects_logit <- dx_effects_logit %>%
    mutate(
      log10_p = -log10(p.value),
      significant = p.value < 0.05
    )
  
  p3 <- ggplot(dx_effects_logit, aes(x = estimate, y = log10_p, color = significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
    labs(title = "PHROG PA Analysis - Volcano Plot (Dx.Status Effect)",
         subtitle = "Total cohort - Logit model",
         x = "Log Odds Ratio (Dx.StatusCELIAC)",
         y = "-log10(P-value)",
         color = "Significant") +
    theme_minimal()
  
  ggsave("../Orf_Contig_Phrog_compositional/report_figures/phrog_PA_volcano_logit.pdf", p3, width = 10, height = 6)
}

# Summary table of top significant PHROGs
if (nrow(dx_sig_logit) > 0) {
  top_phrogs <- dx_sig_logit %>%
    arrange(p.value) %>%
    select(phrog_id, estimate, std.error, p.value) %>%
    head(20)
  
  write.csv(top_phrogs, "../Orf_Contig_Phrog_compositional/phrog_results/top_significant_phrogs_dx_status.csv", row.names = FALSE)
  cat("Top 20 significant PHROGs saved to top_significant_phrogs_dx_status.csv\n")
}

cat("Summary figures saved to report_figures directory\n")