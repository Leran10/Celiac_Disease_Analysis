# Basic visualizations for Celiac Phage Analysis
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
})

cat("Creating basic visualizations...\n")

# Create output directory
dir.create("Orf_Contig_Phrog_compositional/figures", showWarnings = FALSE)

# Load diversity data
if(file.exists("Orf_Contig_Phrog_compositional/results/total_diversity_metrics.csv")) {
  diversity_data <- read.csv("Orf_Contig_Phrog_compositional/results/total_diversity_metrics.csv", row.names = 1)
  
  # Shannon diversity by diagnosis status
  p1 <- ggplot(diversity_data, aes(x = Dx.Status, y = shannon, fill = Dx.Status)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_viridis_d() +
    labs(title = "Shannon Diversity by Diagnosis Status",
         x = "Diagnosis Status", y = "Shannon Diversity") +
    theme_minimal()
  
  ggsave("Orf_Contig_Phrog_compositional/figures/shannon_diversity_by_diagnosis.pdf", p1, width = 8, height = 6)
  
  # Shannon diversity by timepoint
  p2 <- ggplot(diversity_data, aes(x = onset_timeline_combined, y = shannon, fill = Dx.Status)) +
    geom_boxplot() +
    scale_fill_viridis_d() +
    labs(title = "Shannon Diversity by Timepoint and Diagnosis",
         x = "Timepoint", y = "Shannon Diversity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("Orf_Contig_Phrog_compositional/figures/shannon_diversity_by_timepoint.pdf", p2, width = 10, height = 6)
  
  cat("Diversity plots created\n")
}

# Load and visualize timepoint-specific results
timepoint_files <- list.files("Orf_Contig_Phrog_compositional/results/", 
                             pattern = "timepoint_specific_results.csv", 
                             full.names = TRUE)

for(file in timepoint_files) {
  data <- read.csv(file)
  cohort_model <- gsub("_timepoint_specific_results.csv", "", basename(file))
  
  if(nrow(data) > 0) {
    # Volcano plot equivalent - effect size vs significance
    data$neg_log10_p <- -log10(data$p.value)
    
    p3 <- ggplot(data, aes(x = estimate, y = neg_log10_p, color = timepoint)) +
      geom_point(alpha = 0.6) +
      scale_color_viridis_d() +
      labs(title = paste("Effect Sizes by Timepoint -", cohort_model),
           x = "Log Odds Ratio (Estimate)", 
           y = "-log10(p-value)") +
      theme_minimal()
    
    output_file <- paste0("Orf_Contig_Phrog_compositional/figures/", cohort_model, "_effect_sizes.pdf")
    ggsave(output_file, p3, width = 10, height = 6)
  }
}

# Load analysis summary
if(file.exists("Orf_Contig_Phrog_compositional/results/analysis_summary.csv")) {
  summary_data <- read.csv("Orf_Contig_Phrog_compositional/results/analysis_summary.csv")
  
  # Bar plot of number of results by cohort and model
  p4 <- ggplot(summary_data, aes(x = cohort, y = n_results, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis_d() +
    labs(title = "Number of Model Results by Cohort and Model Type",
         x = "Cohort", y = "Number of Results") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("Orf_Contig_Phrog_compositional/figures/analysis_summary_barplot.pdf", p4, width = 10, height = 6)
  
  cat("Summary plot created\n")
}

cat("Basic visualizations completed!\n")
cat("Plots saved in Orf_Contig_Phrog_compositional/figures/\n")