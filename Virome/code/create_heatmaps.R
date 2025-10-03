# Create heatmaps for timepoint-specific differential abundance results
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(viridis)
  library(RColorBrewer)
  library(reshape2)
})

cat("Creating heatmaps for timepoint-specific differential abundance...\n")

# Function to create heatmap from timepoint results
create_timepoint_heatmap <- function(timepoint_file, output_name, title_prefix) {
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
  
  # Filter for significant results (p < 0.05) and create matrix
  sig_data <- data %>%
    filter(p.value < 0.05) %>%
    select(gene, timepoint, estimate, p.value)
  
  if(nrow(sig_data) == 0) {
    cat("No significant results found for", output_name, "\n")
    return(NULL)
  }
  
  # Create matrix for heatmap (genes x timepoints)
  estimate_matrix <- sig_data %>%
    select(gene, timepoint, estimate) %>%
    dcast(gene ~ timepoint, value.var = "estimate", fill = 0)
  
  # Convert to matrix
  rownames(estimate_matrix) <- estimate_matrix$gene
  estimate_matrix <- estimate_matrix[, -1, drop = FALSE]
  estimate_matrix <- as.matrix(estimate_matrix)
  
  # Remove genes with all zeros
  estimate_matrix <- estimate_matrix[rowSums(abs(estimate_matrix)) > 0, , drop = FALSE]
  
  if(nrow(estimate_matrix) == 0) {
    cat("No genes with non-zero estimates for", output_name, "\n")
    return(NULL)
  }
  
  cat("Creating heatmap with", nrow(estimate_matrix), "genes and", ncol(estimate_matrix), "timepoints\n")
  
  # Create heatmap
  pdf(paste0("../Orf_Contig_Phrog_compositional/figures/", output_name, "_heatmap.pdf"), width = 10, height = 12)
  
  pheatmap(estimate_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "none",
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-max(abs(estimate_matrix)), max(abs(estimate_matrix)), length.out = 101),
           main = paste(title_prefix, "- Log Odds Ratios (CELIAC vs CONTROL)"),
           fontsize = 8,
           fontsize_row = 6,
           fontsize_col = 8,
           show_rownames = ifelse(nrow(estimate_matrix) <= 50, TRUE, FALSE),
           border_color = NA)
  
  dev.off()
  
  # Also create a top genes heatmap (top 30 by max absolute effect size)
  if(nrow(estimate_matrix) > 30) {
    max_effects <- apply(abs(estimate_matrix), 1, max)
    top_genes <- names(sort(max_effects, decreasing = TRUE))[1:30]
    top_matrix <- estimate_matrix[top_genes, , drop = FALSE]
    
    pdf(paste0("../Orf_Contig_Phrog_compositional/figures/", output_name, "_top30_heatmap.pdf"), width = 10, height = 8)
    
    pheatmap(top_matrix,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             scale = "none",
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             breaks = seq(-max(abs(top_matrix)), max(abs(top_matrix)), length.out = 101),
             main = paste(title_prefix, "- Top 30 Genes by Effect Size"),
             fontsize = 10,
             fontsize_row = 8,
             fontsize_col = 10,
             show_rownames = TRUE,
             border_color = NA)
    
    dev.off()
    cat("Created top 30 genes heatmap\n")
  }
  
  return(estimate_matrix)
}

# Create heatmaps for all timepoint-specific results
timepoint_files <- list.files("../Orf_Contig_Phrog_compositional/results/", 
                             pattern = "timepoint_specific_results.csv", 
                             full.names = TRUE)

for(file in timepoint_files) {
  base_name <- gsub("_timepoint_specific_results.csv", "", basename(file))
  
  # Extract cohort and model info for title
  parts <- strsplit(base_name, "_")[[1]]
  cohort <- parts[1]
  model_type <- ifelse(grepl("PA", base_name), "PA", "Abundance")
  
  title_prefix <- paste(toupper(cohort), model_type, "Model")
  
  create_timepoint_heatmap(file, base_name, title_prefix)
}

# Create combined heatmap across all cohorts for PA models
cat("Creating combined PA heatmap across all cohorts...\n")

# Load all PA timepoint results
pa_files <- timepoint_files[grepl("PA_model1", timepoint_files)]
combined_pa_data <- list()

for(file in pa_files) {
  cohort <- gsub("_.*", "", basename(file))
  data <- read.csv(file)
  if(nrow(data) > 0) {
    data$cohort <- cohort
    combined_pa_data[[cohort]] <- data
  }
}

if(length(combined_pa_data) > 0) {
  all_pa_data <- do.call(rbind, combined_pa_data)
  
  # Filter significant results and create matrix
  sig_pa_data <- all_pa_data %>%
    filter(p.value < 0.05) %>%
    select(gene, timepoint, cohort, estimate) %>%
    mutate(gene_cohort = paste(gene, cohort, sep = "_"))
  
  if(nrow(sig_pa_data) > 0) {
    pa_matrix <- sig_pa_data %>%
      dcast(gene_cohort ~ timepoint, value.var = "estimate", fill = 0)
    
    rownames(pa_matrix) <- pa_matrix$gene_cohort
    pa_matrix <- pa_matrix[, -1, drop = FALSE]
    pa_matrix <- as.matrix(pa_matrix)
    
    # Remove genes with all zeros
    pa_matrix <- pa_matrix[rowSums(abs(pa_matrix)) > 0, , drop = FALSE]
    
    if(nrow(pa_matrix) > 0) {
      pdf("../Orf_Contig_Phrog_compositional/figures/combined_PA_all_cohorts_heatmap.pdf", width = 12, height = 16)
      
      pheatmap(pa_matrix,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               breaks = seq(-max(abs(pa_matrix)), max(abs(pa_matrix)), length.out = 101),
               main = "PA Models - All Cohorts Combined\nLog Odds Ratios (CELIAC vs CONTROL)",
               fontsize = 8,
               fontsize_row = 4,
               fontsize_col = 10,
               show_rownames = ifelse(nrow(pa_matrix) <= 100, TRUE, FALSE),
               border_color = NA)
      
      dev.off()
      cat("Combined PA heatmap created with", nrow(pa_matrix), "genes\n")
    }
  }
}

# Create summary heatmap of significant gene counts by timepoint and cohort
cat("Creating summary heatmap of significant gene counts...\n")

summary_data <- data.frame()
for(file in timepoint_files) {
  base_name <- gsub("_timepoint_specific_results.csv", "", basename(file))
  data <- read.csv(file)
  
  if(nrow(data) > 0) {
    sig_counts <- data %>%
      filter(p.value < 0.05) %>%
      group_by(timepoint) %>%
      summarise(n_significant = n(), .groups = 'drop') %>%
      mutate(cohort_model = base_name)
    
    summary_data <- rbind(summary_data, sig_counts)
  }
}

if(nrow(summary_data) > 0) {
  summary_matrix <- summary_data %>%
    dcast(cohort_model ~ timepoint, value.var = "n_significant", fill = 0)
  
  rownames(summary_matrix) <- summary_matrix$cohort_model
  summary_matrix <- summary_matrix[, -1, drop = FALSE]
  summary_matrix <- as.matrix(summary_matrix)
  
  pdf("../Orf_Contig_Phrog_compositional/figures/significant_genes_summary_heatmap.pdf", width = 10, height = 8)
  
  pheatmap(summary_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "none",
           color = colorRampPalette(c("white", "red"))(100),
           main = "Number of Significant Genes by Timepoint and Model",
           fontsize = 10,
           fontsize_row = 9,
           fontsize_col = 10,
           show_rownames = TRUE,
           border_color = "grey60",
           display_numbers = TRUE,
           number_format = "%.0f")
  
  dev.off()
  cat("Summary heatmap created\n")
}

cat("All heatmaps completed!\n")
cat("Heatmaps saved in Orf_Contig_Phrog_compositional/figures/\n")