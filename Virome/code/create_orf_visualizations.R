# Create ORF Analysis Visualizations
library(dplyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(tibble)
library(tidyr)

cat("Loading existing PA model results...\n")
total_pa_results <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/total_PA_model1_glmmTMB_logit.csv")
US_pa_results <- read.csv("../Orf_Contig_Phrog_compositional/orf_results/US_PA_model1_glmmTMB_logit.csv")

cat("Total PA results:", nrow(total_pa_results), "rows\n")
cat("US PA results:", nrow(US_pa_results), "rows\n")

# Extract timepoint-specific results function
extract_timepoint_results <- function(model_results, model_name, gene_col = "orf_id") {
  if(is.null(model_results) || nrow(model_results) == 0) {
    return(data.frame())
  }
  
  # Extract interaction terms for timepoints
  timepoint_results <- model_results %>%
    filter(grepl("Dx.StatusCELIAC:onset_timeline_combined", term)) %>%
    mutate(
      timepoint = gsub("Dx.StatusCELIAC:onset_timeline_combined", "", term),
      model = model_name,
      gene = get(gene_col)
    ) %>%
    select(gene, timepoint, estimate, std.error, p.value, conf.low, conf.high, model)
  
  return(timepoint_results)
}

# Extract timepoint results
cat("Extracting timepoint-specific results...\n")
total_pa_timepoints <- extract_timepoint_results(total_pa_results, "PA_glmmTMB_logit")
US_pa_timepoints <- extract_timepoint_results(US_pa_results, "PA_glmmTMB_logit")

# Save timepoint results
write.csv(total_pa_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/total_PA_model1_glmmTMB_logit_timepoint_specific_results.csv", row.names = FALSE)
write.csv(US_pa_timepoints, "../Orf_Contig_Phrog_compositional/orf_results/US_PA_model1_glmmTMB_logit_timepoint_specific_results.csv", row.names = FALSE)

cat("Total timepoint results:", nrow(total_pa_timepoints), "entries\n")
cat("US timepoint results:", nrow(US_pa_timepoints), "entries\n")

# Create summary statistics
total_summary <- total_pa_timepoints %>%
  filter(!is.na(p.value)) %>%
  group_by(timepoint) %>%
  summarise(
    n_genes = n(),
    significant = sum(p.value < 0.05, na.rm = TRUE),
    mean_effect = mean(estimate, na.rm = TRUE),
    .groups = "drop"
  )

US_summary <- US_pa_timepoints %>%
  filter(!is.na(p.value)) %>%
  group_by(timepoint) %>%
  summarise(
    n_genes = n(),
    significant = sum(p.value < 0.05, na.rm = TRUE),
    mean_effect = mean(estimate, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(total_summary, "../Orf_Contig_Phrog_compositional/orf_results/total_PA_timepoint_summary.csv", row.names = FALSE)
write.csv(US_summary, "../Orf_Contig_Phrog_compositional/orf_results/US_PA_timepoint_summary.csv", row.names = FALSE)

cat("Summary statistics saved\n")
print("Total cohort summary:")
print(total_summary)
print("US cohort summary:")
print(US_summary)

# Create heatmap for significant results
create_orf_heatmap <- function(timepoint_data, cohort_name, output_dir = "../Orf_Contig_Phrog_compositional/orf_figures/") {
  # Get significant results
  sig_data <- timepoint_data %>%
    filter(!is.na(p.value) & p.value < 0.05)
  
  if(nrow(sig_data) == 0) {
    cat("No significant results for", cohort_name, "heatmap\n")
    return(NULL)
  }
  
  # Prepare data for heatmap
  timepoint_order <- c("t0-over42", "t0-36", "t0-30", "t0-24", "t0-18", "t0-12", "t0-6", "t0")
  
  heatmap_data <- sig_data %>%
    select(gene, timepoint, estimate) %>%
    mutate(timepoint = factor(timepoint, levels = timepoint_order)) %>%
    filter(timepoint %in% timepoint_order) %>%
    tidyr::pivot_wider(names_from = timepoint, values_from = estimate, values_fill = 0) %>%
    column_to_rownames("gene")
  
  if(nrow(heatmap_data) == 0) {
    cat("No valid timepoint data for", cohort_name, "heatmap\n")
    return(NULL)
  }
  
  # Create heatmap
  pdf(paste0(output_dir, cohort_name, "_PA_model1_orf_heatmap.pdf"), width = 10, height = 8)
  
  pheatmap(
    as.matrix(heatmap_data),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",
    cluster_cols = FALSE,  # Keep chronological order
    cluster_rows = TRUE,
    fontsize = 8,
    main = paste(cohort_name, "Cohort - ORF PA Analysis\nLog Odds Ratios (CELIAC vs CONTROL)"),
    na_col = "grey90"
  )
  
  dev.off()
  
  cat("Heatmap saved for", cohort_name, ":", nrow(heatmap_data), "genes\n")
  return(heatmap_data)
}

# Create heatmaps
dir.create("../Orf_Contig_Phrog_compositional/orf_figures", showWarnings = FALSE, recursive = TRUE)

total_heatmap <- create_orf_heatmap(total_pa_timepoints, "total")
US_heatmap <- create_orf_heatmap(US_pa_timepoints, "US")

# Create effect size plots
create_orf_effect_plot <- function(results_data, cohort_name, output_dir = "../Orf_Contig_Phrog_compositional/orf_figures/") {
  
  # Filter main effects (not interaction terms)
  main_effects <- results_data %>%
    filter(term == "Dx.StatusCELIAC") %>%
    filter(!is.na(estimate))
  
  if(nrow(main_effects) == 0) {
    cat("No main effects data for", cohort_name, "effect plot\n")
    return(NULL)
  }
  
  # Create effect size plot
  p <- ggplot(main_effects, aes(x = estimate, y = -log10(p.value + 1e-10))) +
    geom_point(aes(color = p.value < 0.05), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), 
                       name = "Significant\n(p < 0.05)") +
    labs(
      title = paste(cohort_name, "Cohort - ORF PA Effect Sizes"),
      subtitle = "Main effect: CELIAC vs CONTROL",
      x = "Log Odds Ratio (Effect Size)",
      y = "-log10(p-value)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  
  # Save plot
  ggsave(paste0(output_dir, cohort_name, "_PA_model1_orf_effect_sizes.pdf"), 
         p, width = 10, height = 8, device = "pdf")
  
  cat("Effect size plot saved for", cohort_name, ":", nrow(main_effects), "effects\n")
  
  return(p)
}

# Create effect size plots
total_effect_plot <- create_orf_effect_plot(total_pa_results, "total")
US_effect_plot <- create_orf_effect_plot(US_pa_results, "US")

# Create analysis summary plot
create_summary_barplot <- function(total_summary, US_summary, output_dir = "../Orf_Contig_Phrog_compositional/orf_figures/") {
  
  # Combine summaries
  combined_summary <- bind_rows(
    total_summary %>% mutate(cohort = "Total"),
    US_summary %>% mutate(cohort = "US")
  )
  
  if(nrow(combined_summary) == 0) {
    cat("No summary data for barplot\n")
    return(NULL)
  }
  
  # Create barplot
  p <- ggplot(combined_summary, aes(x = timepoint, y = significant, fill = cohort)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Total" = "#2E8B57", "US" = "#4682B4")) +
    labs(
      title = "ORF-based PA Analysis - Significant Results by Timepoint",
      subtitle = "Number of ORFs with significant effects (p < 0.05)",
      x = "Timepoint (months before diagnosis)",
      y = "Number of Significant ORFs",
      fill = "Cohort"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save plot
  ggsave(paste0(output_dir, "orf_analysis_summary_barplot.pdf"), 
         p, width = 12, height = 8, device = "pdf")
  
  cat("Summary barplot saved\n")
  return(p)
}

# Create summary plot
if(nrow(total_summary) > 0 || nrow(US_summary) > 0) {
  summary_plot <- create_summary_barplot(total_summary, US_summary)
}

cat("\n=== ORF VISUALIZATION CREATION COMPLETED ===\n")
cat("Files created in: ../Orf_Contig_Phrog_compositional/orf_figures/\n")
cat("Results files in: ../Orf_Contig_Phrog_compositional/orf_results/\n")