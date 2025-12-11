# Create comprehensive PHROG figures showing all completed models
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(viridis)
library(pheatmap)

# Read all available PHROG analysis results
cat("Reading PHROG analysis results...\n")

# PA Models
total_logit <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_logit_results.csv")
total_cloglog <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_cloglog_results.csv")
total_glmer <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmer_results.csv")

# Abundance Models
total_abundance_limma <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_limma_results.csv")
total_abundance_negbinom <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_abundance_negbinom_subset_results.csv")

cat("PA Models loaded:\n")
cat("- glmmTMB logit:", nrow(total_logit), "rows\n")
cat("- glmmTMB cloglog:", nrow(total_cloglog), "rows\n") 
cat("- glmer:", nrow(total_glmer), "rows\n")

cat("Abundance Models loaded:\n")
cat("- limma-voom:", nrow(total_abundance_limma), "rows\n")
cat("- Negative Binomial:", nrow(total_abundance_negbinom), "rows\n")

# Figure 1: Model Comparison for PA Analysis - All Models
pa_models_comparison <- data.frame(
  Model = c("glmmTMB (logit)", "glmmTMB (cloglog)", "glmer"),
  Total_Results = c(nrow(total_logit), nrow(total_cloglog), nrow(total_glmer)),
  Dx_Status_Results = c(
    nrow(filter(total_logit, grepl("Dx.StatusCELIAC", term) & !is.na(p.value))),
    nrow(filter(total_cloglog, grepl("Dx.StatusCELIAC", term) & !is.na(p.value))),
    nrow(filter(total_glmer, grepl("Dx.StatusCELIAC", term) & !is.na(p.value)))
  ),
  Significant_Dx_Status = c(
    nrow(filter(total_logit, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05)),
    nrow(filter(total_cloglog, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05)),
    nrow(filter(total_glmer, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05))
  )
)

# Reshape for plotting
pa_comparison_long <- pa_models_comparison %>%
  select(Model, Dx_Status_Results, Significant_Dx_Status) %>%
  pivot_longer(cols = c(Dx_Status_Results, Significant_Dx_Status), 
               names_to = "Result_Type", values_to = "Count") %>%
  mutate(Result_Type = case_when(
    Result_Type == "Dx_Status_Results" ~ "Total Dx.Status Results",
    Result_Type == "Significant_Dx_Status" ~ "Significant Dx.Status (p<0.05)"
  ))

p1 <- ggplot(pa_comparison_long, aes(x = Model, y = Count, fill = Result_Type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_manual(values = c("Total Dx.Status Results" = "lightblue", 
                              "Significant Dx.Status (p<0.05)" = "darkblue")) +
  labs(title = "PHROG PA Analysis: Model Comparison",
       subtitle = "Comparison across different mixed-effects approaches",
       x = "Statistical Model",
       y = "Number of Results",
       fill = "Result Type") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_model_comparison.pdf", 
       p1, width = 12, height = 8)

# Figure 2: Abundance Model Comparison
abundance_models_comparison <- data.frame(
  Model = c("limma-voom", "Negative Binomial"),
  Total_Coefficients = c(nrow(total_abundance_limma), nrow(total_abundance_negbinom)),
  Dx_Status_Results = c(
    nrow(filter(total_abundance_limma, coefficient == "Dx.StatusCELIAC" & !is.na(P.Value))),
    nrow(filter(total_abundance_negbinom, grepl("Dx.StatusCELIAC", term) & !is.na(p.value)))
  ),
  Significant_Dx_Status = c(
    nrow(filter(total_abundance_limma, coefficient == "Dx.StatusCELIAC" & !is.na(P.Value) & P.Value < 0.05)),
    nrow(filter(total_abundance_negbinom, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05))
  )
)

abundance_comparison_long <- abundance_models_comparison %>%
  select(Model, Dx_Status_Results, Significant_Dx_Status) %>%
  pivot_longer(cols = c(Dx_Status_Results, Significant_Dx_Status), 
               names_to = "Result_Type", values_to = "Count") %>%
  mutate(Result_Type = case_when(
    Result_Type == "Dx_Status_Results" ~ "Total Dx.Status Results",
    Result_Type == "Significant_Dx_Status" ~ "Significant Dx.Status (p<0.05)"
  ))

p2 <- ggplot(abundance_comparison_long, aes(x = Model, y = Count, fill = Result_Type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_manual(values = c("Total Dx.Status Results" = "lightgreen", 
                              "Significant Dx.Status (p<0.05)" = "darkgreen")) +
  labs(title = "PHROG Abundance Analysis: Model Comparison", 
       subtitle = "Comparison across different abundance modeling approaches",
       x = "Statistical Model",
       y = "Number of Results",
       fill = "Result Type") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_abundance_model_comparison.pdf", 
       p2, width = 12, height = 8)

# Figure 3: Heatmap for PA Model 2 (cloglog)
dx_effects_cloglog <- total_cloglog %>% 
  filter(grepl("Dx.StatusCELIAC", term)) %>%
  filter(p.value < 0.05) %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  head(50)

if(nrow(dx_effects_cloglog) > 0) {
  p3 <- ggplot(dx_effects_cloglog, aes(x = 1, y = reorder(phrog_id, -estimate), fill = estimate)) +
    geom_tile(color = "white", width = 0.8) +
    geom_text(aes(label = paste0("p=", round(p.value, 4))), size = 2.5, color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Log Odds\nRatio") +
    labs(title = "PHROG PA Model 2 (cloglog): Top 50 Celiac vs Control Effects",
         subtitle = "Complementary log-log model results ranked by significance",
         x = "",
         y = "PHROG ID",
         caption = "Blue = Lower in Celiac, Red = Higher in Celiac") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_model2_cloglog_heatmap.pdf", 
         p3, width = 8, height = 16)
}

# Figure 4: Heatmap for PA Model 3 (glmer)
dx_effects_glmer <- total_glmer %>% 
  filter(grepl("Dx.StatusCELIAC", term)) %>%
  filter(p.value < 0.05) %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  head(50)

if(nrow(dx_effects_glmer) > 0) {
  p4 <- ggplot(dx_effects_glmer, aes(x = 1, y = reorder(phrog_id, -estimate), fill = estimate)) +
    geom_tile(color = "white", width = 0.8) +
    geom_text(aes(label = paste0("p=", round(p.value, 4))), size = 2.5, color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Log Odds\nRatio") +
    labs(title = "PHROG PA Model 3 (glmer): Top 50 Celiac vs Control Effects",
         subtitle = "lme4 mixed-effects model results ranked by significance",
         x = "",
         y = "PHROG ID",
         caption = "Blue = Lower in Celiac, Red = Higher in Celiac") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_model3_glmer_heatmap.pdf", 
         p4, width = 8, height = 16)
}

# Figure 5: Abundance Model 2 (Negative Binomial) Heatmap
dx_effects_negbinom <- total_abundance_negbinom %>% 
  filter(grepl("Dx.StatusCELIAC", term)) %>%
  filter(p.value < 0.05) %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  head(min(50, nrow(.)))

if(nrow(dx_effects_negbinom) > 0) {
  p5 <- ggplot(dx_effects_negbinom, aes(x = 1, y = reorder(phrog_id, -estimate), fill = estimate)) +
    geom_tile(color = "white", width = 0.8) +
    geom_text(aes(label = paste0("p=", round(p.value, 4))), size = 2.5, color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Log Fold\nChange") +
    labs(title = "PHROG Abundance Model 2 (Negative Binomial): Top Celiac vs Control Effects",
         subtitle = "glmmTMB negative binomial model results (subset)",
         x = "",
         y = "PHROG ID",
         caption = "Blue = Lower in Celiac, Red = Higher in Celiac") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_abundance_model2_negbinom_heatmap.pdf", 
         p5, width = 8, height = 12)
}

# Figure 6: Complete Model Overview
all_models_summary <- data.frame(
  Analysis_Type = c(rep("Presence/Absence", 3), rep("Abundance", 2)),
  Model = c("glmmTMB (logit)", "glmmTMB (cloglog)", "glmer", "limma-voom", "Negative Binomial"),
  Status = c("✅ Complete", "✅ Complete", "✅ Complete", "✅ Complete", "✅ Subset Complete"),
  Total_Results = c(nrow(total_logit), nrow(total_cloglog), nrow(total_glmer), 
                   nrow(total_abundance_limma), nrow(total_abundance_negbinom)),
  Significant_Dx_Status = c(
    nrow(filter(total_logit, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05)),
    nrow(filter(total_cloglog, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05)),
    nrow(filter(total_glmer, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05)),
    nrow(filter(total_abundance_limma, coefficient == "Dx.StatusCELIAC" & !is.na(P.Value) & P.Value < 0.05)),
    nrow(filter(total_abundance_negbinom, grepl("Dx.StatusCELIAC", term) & !is.na(p.value) & p.value < 0.05))
  )
)

p6 <- ggplot(all_models_summary, aes(x = Model, y = Significant_Dx_Status, fill = Analysis_Type)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = Significant_Dx_Status), vjust = -0.5, fontface = "bold") +
  scale_fill_manual(values = c("Presence/Absence" = "lightblue", "Abundance" = "lightgreen")) +
  facet_wrap(~Analysis_Type, scales = "free_x") +
  labs(title = "PHROG Analysis: Complete Model Results Summary",
       subtitle = "Significant Dx.Status effects across all implemented models",
       x = "Statistical Model",
       y = "Number of Significant PHROGs (p < 0.05)",
       fill = "Analysis Type") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"))

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_complete_model_overview.pdf", 
       p6, width = 14, height = 8)

cat("Created comprehensive PHROG model comparison figures:\n")
cat("11. PHROG_PA_model_comparison.pdf\n")
cat("12. PHROG_abundance_model_comparison.pdf\n")
cat("13. PHROG_PA_model2_cloglog_heatmap.pdf\n")
cat("14. PHROG_PA_model3_glmer_heatmap.pdf\n")
cat("15. PHROG_abundance_model2_negbinom_heatmap.pdf\n")
cat("16. PHROG_complete_model_overview.pdf\n")
cat("All figures demonstrate multiple model implementation!\n")