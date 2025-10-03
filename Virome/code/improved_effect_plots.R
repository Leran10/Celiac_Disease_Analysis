# Create improved effect size visualizations that show the full picture
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(viridis)
})

cat("Creating improved effect size visualizations...\n")

# Function to create comprehensive effect size plots
create_improved_effect_plots <- function(timepoint_file, output_name, title_prefix) {
  if(!file.exists(timepoint_file)) {
    cat("File not found:", timepoint_file, "\n")
    return(NULL)
  }
  
  # Load timepoint results
  data <- read.csv(timepoint_file)
  
  if(nrow(data) == 0) {
    cat("No data in file:", timepoint_file, "\n")
    return(NULL)
  }
  
  cat("Processing", nrow(data), "results from", basename(timepoint_file), "\n")
  
  # Prepare data for plotting
  data$has_pvalue <- !is.na(data$p.value)
  data$neg_log10_p <- ifelse(data$has_pvalue, -log10(data$p.value), NA)
  data$is_significant <- data$has_pvalue & data$p.value < 0.05
  
  # Summary statistics
  n_total <- nrow(data)
  n_with_pvalue <- sum(data$has_pvalue)
  n_significant <- sum(data$is_significant, na.rm = TRUE)
  
  cat("Data summary for", output_name, ":\n")
  cat("  Total results:", n_total, "\n")
  cat("  With valid p-values:", n_with_pvalue, "(", round(100*n_with_pvalue/n_total, 1), "%)\n")
  cat("  Significant (p<0.05):", n_significant, "\n")
  cat("  Effect size range:", round(min(data$estimate, na.rm=TRUE), 2), "to", round(max(data$estimate, na.rm=TRUE), 2), "\n")
  
  # Create multiple plots to show the full picture
  
  # Plot 1: All effect sizes (including NAs)
  p1 <- ggplot(data, aes(x = estimate, y = timepoint)) +
    geom_point(aes(color = has_pvalue, shape = is_significant), 
               alpha = 0.7, size = 2) +
    scale_color_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#F24236"),
                       labels = c("TRUE" = "Valid p-value", "FALSE" = "No p-value (NA)"),
                       name = "Statistical Status") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
                       name = "Significance",
                       na.value = 4) +
    labs(title = paste(title_prefix, "- All Effect Sizes"),
         subtitle = paste0("Total: ", n_total, " | Valid p-values: ", n_with_pvalue, 
                          " (", round(100*n_with_pvalue/n_total, 1), "%)"),
         x = "Effect Size (Log Odds/Fold Change)",
         y = "Timepoint") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9, color = "gray60")) +
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3)))
  
  # Plot 2: Distribution of effect sizes
  p2 <- ggplot(data, aes(x = estimate)) +
    geom_histogram(aes(fill = has_pvalue), bins = 30, alpha = 0.7, position = "stack") +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#F24236"),
                      labels = c("TRUE" = "Valid p-value", "FALSE" = "No p-value (NA)"),
                      name = "Statistical Status") +
    labs(title = "Distribution of Effect Sizes",
         subtitle = "Red bars show large effects with no p-values (model convergence issues)",
         x = "Effect Size (Log Odds/Fold Change)",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9, color = "gray60"))
  
  # Plot 3: Traditional volcano plot (only valid p-values)
  data_valid <- data[data$has_pvalue, ]
  
  if(nrow(data_valid) > 0) {
    p3 <- ggplot(data_valid, aes(x = estimate, y = neg_log10_p)) +
      geom_point(aes(color = is_significant, shape = timepoint), 
                 alpha = 0.7, size = 2) +
      scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                         labels = c("TRUE" = "Significant (p<0.05)", "FALSE" = "Not significant"),
                         name = "Significance") +
      labs(title = "Traditional Volcano Plot (Valid P-values Only)",
           subtitle = paste0("Showing ", nrow(data_valid), " out of ", n_total, " results"),
           x = "Effect Size (Log Odds/Fold Change)",
           y = "-log10(p-value)") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 11),
            plot.subtitle = element_text(size = 9, color = "gray60"))
  } else {
    # Create empty plot if no valid p-values
    p3 <- ggplot() + 
      annotate("text", x = 0, y = 0, 
               label = "No valid p-values found\n(All models failed to converge)", 
               size = 5, color = "red") +
      labs(title = "Traditional Volcano Plot (Valid P-values Only)",
           subtitle = "No convergent models found",
           x = "Effect Size", y = "-log10(p-value)") +
      theme_void()
  }
  
  # Combine plots
  combined_plot <- plot_grid(p1, p2, p3, 
                            ncol = 1, 
                            rel_heights = c(1, 0.8, 1),
                            align = "v")
  
  # Save improved effect size plot
  output_file <- paste0("../Orf_Contig_Phrog_compositional/figures/", output_name, "_improved_effects.pdf")
  ggsave(output_file, combined_plot, width = 10, height = 12, device = "pdf")
  
  cat("Improved effect size plot saved:", basename(output_file), "\n")
  
  return(list(
    total_results = n_total,
    valid_pvalues = n_with_pvalue,
    significant = n_significant,
    effect_range = c(min(data$estimate, na.rm=TRUE), max(data$estimate, na.rm=TRUE))
  ))
}

# Create improved plots for all timepoint-specific results
timepoint_files <- list.files("../Orf_Contig_Phrog_compositional/results/", 
                             pattern = "timepoint_specific_results.csv", 
                             full.names = TRUE)

summary_stats <- list()

for(file in timepoint_files) {
  base_name <- gsub("_timepoint_specific_results.csv", "", basename(file))
  
  # Extract cohort and model info for title
  parts <- strsplit(base_name, "_")[[1]]
  cohort <- parts[1]
  model_type <- ifelse(grepl("PA", base_name), "PA", "Abundance")
  
  title_prefix <- paste(toupper(cohort), model_type, "Model")
  
  stats <- create_improved_effect_plots(file, base_name, title_prefix)
  if(!is.null(stats)) {
    summary_stats[[base_name]] <- stats
  }
}

# Create summary table
if(length(summary_stats) > 0) {
  cat("\n=== EFFECT SIZE ANALYSIS SUMMARY ===\n")
  for(name in names(summary_stats)) {
    stats <- summary_stats[[name]]
    cat(sprintf("%-40s: %3d total | %3d valid p-values (%5.1f%%) | %2d significant | Effects: %6.2f to %6.2f\n",
                name,
                stats$total_results,
                stats$valid_pvalues,
                100 * stats$valid_pvalues / stats$total_results,
                stats$significant,
                stats$effect_range[1],
                stats$effect_range[2]))
  }
}

cat("\n✅ Improved effect size visualizations completed!\n")
cat("New files saved with '_improved_effects' suffix\n")
cat("\nKey insights:\n")
cat("- Red points/bars show large effects that couldn't get p-values (sparse data)\n")
cat("- Blue points show effects with valid statistical inference\n")
cat("- The 'clustering near zero' was due to hiding the large effects!\n")