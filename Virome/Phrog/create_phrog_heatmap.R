# Create PHROG Analysis Heatmap
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(viridis)

# Read the PHROG analysis results
total_logit <- read.csv("../Orf_Contig_Phrog_compositional/phrog_results/total_PA_glmmTMB_logit_results.csv")

cat("Creating PHROG analysis heatmap...\n")

# Focus on significant results for main effects
sig_results <- total_logit %>%
  filter(p.value < 0.05) %>%
  filter(!is.na(p.value)) %>%
  filter(grepl("Dx.StatusCELIAC|onset_timeline", term)) %>%
  select(phrog_id, term, estimate, p.value) %>%
  # Get top 30 most significant PHROGs by minimum p-value per PHROG
  group_by(phrog_id) %>%
  mutate(min_p = min(p.value)) %>%
  ungroup() %>%
  filter(phrog_id %in% (filter(., min_p <= sort(unique(min_p))[30]) %>% pull(phrog_id) %>% unique()))

cat("Selected", length(unique(sig_results$phrog_id)), "top significant PHROGs\n")
cat("Covering", length(unique(sig_results$term)), "model terms\n")

# Clean up term names for better display
sig_results$term_clean <- case_when(
  sig_results$term == "Dx.StatusCELIAC" ~ "Celiac vs Control",
  grepl("onset_timeline_combinedt0-6", sig_results$term) ~ "t0-6 vs t0",
  grepl("onset_timeline_combinedt0-12", sig_results$term) ~ "t0-12 vs t0", 
  grepl("onset_timeline_combinedt0-18", sig_results$term) ~ "t0-18 vs t0",
  grepl("onset_timeline_combinedt0-24", sig_results$term) ~ "t0-24 vs t0",
  grepl("onset_timeline_combinedt0-30", sig_results$term) ~ "t0-30 vs t0",
  grepl("onset_timeline_combinedt0-36", sig_results$term) ~ "t0-36 vs t0",
  grepl("onset_timeline_combinedt0-over42", sig_results$term) ~ "t0-over42 vs t0",
  TRUE ~ sig_results$term
)

# Create significance indicator
sig_results$significance <- case_when(
  sig_results$p.value < 0.001 ~ "***",
  sig_results$p.value < 0.01 ~ "**", 
  sig_results$p.value < 0.05 ~ "*",
  TRUE ~ ""
)

# Create the heatmap
p_heatmap <- ggplot(sig_results, aes(x = term_clean, y = reorder(phrog_id, min_p), fill = estimate)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = significance), color = "white", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                      name = "Log Odds\nRatio") +
  labs(title = "PHROG Presence/Absence Analysis: Top Significant Effects Heatmap",
       subtitle = paste("Top", length(unique(sig_results$phrog_id)), "PHROGs with most significant effects (p < 0.05)"),
       x = "Model Terms",
       y = "PHROG ID",
       caption = "Color indicates effect size; symbols indicate significance: *** p<0.001, ** p<0.01, * p<0.05") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank()
  )

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_heatmap_top_significant.pdf", 
       p_heatmap, width = 12, height = 16)

# Create a focused heatmap just for Dx.Status effects
dx_effects <- total_logit %>%
  filter(term == "Dx.StatusCELIAC") %>%
  filter(p.value < 0.05) %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  head(50)  # Top 50 most significant

p_dx_heatmap <- ggplot(dx_effects, aes(x = 1, y = reorder(phrog_id, -estimate), fill = estimate)) +
  geom_tile(color = "white", width = 0.8) +
  geom_text(aes(label = paste0("p=", round(p.value, 4))), size = 2.5, color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                      name = "Log Odds\nRatio") +
  labs(title = "PHROG Presence/Absence: Top 50 Celiac vs Control Effects",
       subtitle = "Ranked by statistical significance (all p < 0.05)",
       x = "",
       y = "PHROG ID",
       caption = "Blue = Lower in Celiac, Red = Higher in Celiac") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank()
  )

ggsave("../Orf_Contig_Phrog_compositional/phrog_figures/PHROG_PA_dx_status_top50_heatmap.pdf", 
       p_dx_heatmap, width = 8, height = 16)

cat("Created additional heatmap figures:\n")
cat("9. PHROG_PA_heatmap_top_significant.pdf\n")
cat("10. PHROG_PA_dx_status_top50_heatmap.pdf\n")
cat("All figures saved to ../Orf_Contig_Phrog_compositional/phrog_figures/\n")