# Create Comprehensive PHROG Analysis Figures
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(viridis)
library(pheatmap)

# Read the PHROG analysis results
total_logit <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_logit_results.csv")
total_cloglog <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_cloglog_results.csv")
total_abundance <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_limma_results.csv")

cat("Creating comprehensive PHROG analysis figures...\n")

# Figure 1: Effect Size Distribution for Dx.Status (PA Analysis)
dx_effects_logit <- total_logit %>% 
  filter(grepl("Dx.StatusCELIAC", term)) %>%
  filter(!is.na(p.value))

p1 <- ggplot(dx_effects_logit, aes(x = estimate)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "PHROG Presence/Absence Analysis: Celiac vs Control Effect Sizes",
       subtitle = paste("Total cohort - glmmTMB logit model -", nrow(dx_effects_logit), "PHROGs analyzed"),
       x = "Log Odds Ratio (Celiac vs Control)",
       y = "Number of PHROGs",
       caption = "Red line indicates no effect (log OR = 0)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_effect_sizes_dx_status.pdf", 
       p1, width = 12, height = 8)

# Figure 2: P-value Distribution with Significance Highlighting
dx_effects_logit$significant <- dx_effects_logit$p.value < 0.05

p2 <- ggplot(dx_effects_logit, aes(x = p.value, fill = significant)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "white") +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "red"), 
                    name = "Significant", labels = c("No", "Yes (p < 0.05)")) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "PHROG Presence/Absence Analysis: P-value Distribution",
       subtitle = paste("Total cohort -", sum(dx_effects_logit$significant), "significant out of", nrow(dx_effects_logit), "PHROGs"),
       x = "P-value",
       y = "Number of PHROGs",
       caption = "Red line indicates significance threshold (p = 0.05)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom")

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_pvalue_distribution.pdf", 
       p2, width = 12, height = 8)

# Figure 3: Volcano Plot (Effect Size vs Significance)
dx_effects_logit$log10_p <- -log10(dx_effects_logit$p.value)
dx_effects_logit$effect_category <- case_when(
  dx_effects_logit$p.value >= 0.05 ~ "Not Significant",
  dx_effects_logit$estimate > 0 ~ "Higher in Celiac",
  dx_effects_logit$estimate < 0 ~ "Lower in Celiac"
)

p3 <- ggplot(dx_effects_logit, aes(x = estimate, y = log10_p, color = effect_category)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Not Significant" = "gray70", 
                               "Higher in Celiac" = "red", 
                               "Lower in Celiac" = "blue"),
                    name = "PHROG Category") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  labs(title = "PHROG Presence/Absence Analysis: Volcano Plot",
       subtitle = "Effect Size vs Statistical Significance (Celiac vs Control)",
       x = "Log Odds Ratio (Celiac vs Control)",
       y = "-log10(P-value)",
       caption = "Dashed lines indicate significance threshold (p = 0.05) and no effect") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom")

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_volcano_plot.pdf", 
       p3, width = 12, height = 8)

# Figure 4: Abundance Analysis - Dx.Status Effects
dx_abundance_effects <- total_abundance %>% 
  filter(coefficient == "Dx.StatusCELIAC") %>%
  filter(!is.na(P.Value))

p4 <- ggplot(dx_abundance_effects, aes(x = logFC)) +
  geom_histogram(bins = 50, fill = "forestgreen", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "PHROG Abundance Analysis: Celiac vs Control Effect Sizes",
       subtitle = paste("Total cohort - limma-voom -", nrow(dx_abundance_effects), "PHROGs analyzed"),
       x = "Log Fold Change (Celiac vs Control)",
       y = "Number of PHROGs",
       caption = "Red line indicates no change (log FC = 0)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_abundance_effect_sizes_dx_status.pdf", 
       p4, width = 12, height = 8)

# Figure 5: Abundance Analysis Volcano Plot
dx_abundance_effects$log10_p <- -log10(dx_abundance_effects$P.Value)
dx_abundance_effects$significant <- dx_abundance_effects$P.Value < 0.05
dx_abundance_effects$effect_category <- case_when(
  dx_abundance_effects$P.Value >= 0.05 ~ "Not Significant",
  dx_abundance_effects$logFC > 0 ~ "Higher in Celiac",
  dx_abundance_effects$logFC < 0 ~ "Lower in Celiac"
)

p5 <- ggplot(dx_abundance_effects, aes(x = logFC, y = log10_p, color = effect_category)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Not Significant" = "gray70", 
                               "Higher in Celiac" = "red", 
                               "Lower in Celiac" = "blue"),
                    name = "PHROG Category") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  labs(title = "PHROG Abundance Analysis: Volcano Plot",
       subtitle = "Log Fold Change vs Statistical Significance (Celiac vs Control)",
       x = "Log Fold Change (Celiac vs Control)",
       y = "-log10(P-value)",
       caption = "Dashed lines indicate significance threshold (p = 0.05) and no change") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom")

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_abundance_volcano_plot.pdf", 
       p5, width = 12, height = 8)

# Figure 6: Comparison of PA vs Abundance Results
# Get top significant PHROGs from both analyses
top_pa_phrogs <- dx_effects_logit %>% 
  filter(significant) %>% 
  arrange(p.value) %>% 
  head(20) %>%
  pull(phrog_id)

top_abundance_phrogs <- dx_abundance_effects %>% 
  filter(significant) %>% 
  arrange(P.Value) %>% 
  head(20) %>%
  pull(phrog_id)

# Compare overlap
overlap_phrogs <- intersect(top_pa_phrogs, top_abundance_phrogs)
pa_only <- setdiff(top_pa_phrogs, top_abundance_phrogs)
abundance_only <- setdiff(top_abundance_phrogs, top_pa_phrogs)

comparison_data <- data.frame(
  Category = c("PA Only", "Abundance Only", "Both Analyses", "Total Unique"),
  Count = c(length(pa_only), length(abundance_only), length(overlap_phrogs), 
            length(unique(c(top_pa_phrogs, top_abundance_phrogs))))
)

p6 <- ggplot(comparison_data, aes(x = Category, y = Count, fill = Category)) +
  geom_col(alpha = 0.8, color = "white", linewidth = 1) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "PHROG Analysis: Top 20 Significant PHROGs Comparison",
       subtitle = "Overlap between Presence/Absence and Abundance Analyses",
       x = "Analysis Type",
       y = "Number of PHROGs",
       caption = "Comparison of top 20 significant PHROGs from each analysis") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  ylim(0, max(comparison_data$Count) * 1.1)

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_analysis_comparison.pdf", 
       p6, width = 10, height = 8)

# Figure 7: Timeline Effects (across all timepoints)
timeline_effects <- total_logit %>%
  filter(grepl("onset_timeline_combined", term)) %>%
  filter(!is.na(p.value)) %>%
  mutate(
    timepoint = gsub("onset_timeline_combined", "", term),
    significant = p.value < 0.05
  )

timeline_summary <- timeline_effects %>%
  group_by(timepoint) %>%
  summarise(
    total_effects = n(),
    significant_effects = sum(significant, na.rm = TRUE),
    prop_significant = significant_effects / total_effects
  )

p7 <- ggplot(timeline_summary, aes(x = timepoint, y = significant_effects)) +
  geom_col(fill = "orange", alpha = 0.8, color = "white") +
  geom_text(aes(label = significant_effects), vjust = -0.5, fontface = "bold") +
  labs(title = "PHROG Presence/Absence: Temporal Effects Across Disease Timeline",
       subtitle = "Number of significant PHROG effects at each timepoint vs diagnosis (t0)",
       x = "Timepoint",
       y = "Number of Significant PHROG Effects",
       caption = "All timepoints compared to t0 (diagnosis)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, max(timeline_summary$significant_effects) * 1.1)

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_temporal_effects.pdf", 
       p7, width = 12, height = 8)

# Figure 8: Summary Statistics
summary_stats <- data.frame(
  Metric = c("Total PHROGs Analyzed", "PA - Significant Dx.Status Effects", 
             "Abundance - Significant Dx.Status Effects", "Total Significant PHROGs",
             "PA - Total Coefficients", "Abundance - Total Coefficients"),
  Value = c(1111, 
            sum(dx_effects_logit$significant),
            sum(dx_abundance_effects$significant),
            length(unique(c(dx_effects_logit$phrog_id[dx_effects_logit$significant],
                           dx_abundance_effects$phrog_id[dx_abundance_effects$significant]))),
            nrow(total_logit),
            nrow(total_abundance))
)

p8 <- ggplot(summary_stats, aes(x = reorder(Metric, Value), y = Value)) +
  geom_col(fill = "purple", alpha = 0.7, color = "white") +
  geom_text(aes(label = Value), hjust = -0.1, fontface = "bold") +
  coord_flip() +
  labs(title = "PHROG Analysis: Summary Statistics",
       subtitle = "Comprehensive analysis results for total cohort",
       x = "Metric",
       y = "Count",
       caption = "PA = Presence/Absence analysis; Dx.Status = Celiac vs Control comparison") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  expand_limits(y = max(summary_stats$Value) * 1.1)

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_summary_statistics.pdf", 
       p8, width = 12, height = 8)

cat("Created 8 comprehensive PHROG analysis figures:\n")
cat("1. PHROG_PA_effect_sizes_dx_status.pdf\n")
cat("2. PHROG_PA_pvalue_distribution.pdf\n")
cat("3. PHROG_PA_volcano_plot.pdf\n")
cat("4. PHROG_abundance_effect_sizes_dx_status.pdf\n")
cat("5. PHROG_abundance_volcano_plot.pdf\n")
cat("6. PHROG_analysis_comparison.pdf\n")
cat("7. PHROG_temporal_effects.pdf\n")
cat("8. PHROG_summary_statistics.pdf\n")
cat("All figures saved to ../Orf_Contig_Phrog_compositional/phrog_figures/\n")