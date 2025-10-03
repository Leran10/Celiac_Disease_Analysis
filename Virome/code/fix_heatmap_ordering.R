# Fix heatmap timepoint ordering
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(viridis)
  library(RColorBrewer)
  library(reshape2)
})

cat("Creating heatmaps with proper timepoint ordering...\n")

# Define correct timepoint order
timepoint_order <- c("t0", "t0-6", "t0-12", "t0-18", "t0-24", "t0-30", "t0-36", "t0-over42")

# Function to create heatmap with proper column ordering
create_ordered_heatmap <- function(timepoint_file, output_name, title_prefix) {
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
  
  # IMPORTANT: Reorder columns to match timepoint order
  available_timepoints <- intersect(timepoint_order, colnames(estimate_matrix))
  if(length(available_timepoints) > 0) {
    estimate_matrix <- estimate_matrix[, available_timepoints, drop = FALSE]
  }
  
  cat("Creating ordered heatmap with", nrow(estimate_matrix), "genes and", ncol(estimate_matrix), "timepoints\n")
  cat("Timepoint order:", colnames(estimate_matrix), "\n")
  
  # Calculate dimensions for compact layout
  n_genes <- nrow(estimate_matrix)
  n_timepoints <- ncol(estimate_matrix)
  
  # Compact sizing
  cell_height <- 0.3
  cell_width <- 0.8
  
  fig_height <- max(4, min(12, n_genes * cell_height + 2))
  fig_width <- max(6, n_timepoints * cell_width + 3)
  
  # Create ordered heatmap
  pdf(paste0("../Orf_Contig_Phrog_compositional/figures/", output_name, "_ordered_heatmap.pdf"), 
      width = fig_width, height = fig_height)
  
  pheatmap(estimate_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,  # KEEP FALSE to maintain timepoint order
           scale = "none",
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-max(abs(estimate_matrix)), max(abs(estimate_matrix)), length.out = 101),
           main = paste(title_prefix, "- Log Odds/Fold Change (CELIAC vs CONTROL)\nTimepoints in Chronological Order"),
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 10,
           show_rownames = TRUE,
           border_color = "white",
           cellwidth = 25,
           cellheight = 12,
           treeheight_row = 30,
           treeheight_col = 0,  # Remove column dendrogram since we want chronological order
           angle_col = 45)      # Angle column labels for better readability
  
  dev.off()
  
  return(estimate_matrix)
}

# Create ordered heatmaps for all timepoint-specific results
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
  
  create_ordered_heatmap(file, base_name, title_prefix)
}

# Create ordered combined heatmap across all cohorts for PA models
cat("Creating ordered combined PA heatmap across all cohorts...\n")

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
      # IMPORTANT: Reorder columns to match timepoint order
      available_timepoints <- intersect(timepoint_order, colnames(pa_matrix))
      if(length(available_timepoints) > 0) {
        pa_matrix <- pa_matrix[, available_timepoints, drop = FALSE]
      }
      
      cat("Combined PA timepoint order:", colnames(pa_matrix), "\n")
      
      # Compact dimensions for combined plot
      n_genes <- nrow(pa_matrix)
      n_timepoints <- ncol(pa_matrix)
      fig_height <- max(5, min(10, n_genes * 0.25 + 2))
      fig_width <- max(7, n_timepoints * 0.8 + 3)
      
      pdf("../Orf_Contig_Phrog_compositional/figures/combined_PA_all_cohorts_ordered_heatmap.pdf", 
          width = fig_width, height = fig_height)
      
      pheatmap(pa_matrix,
               cluster_rows = TRUE,
               cluster_cols = FALSE,  # KEEP FALSE to maintain timepoint order
               scale = "none",
               color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               breaks = seq(-max(abs(pa_matrix)), max(abs(pa_matrix)), length.out = 101),
               main = "PA Models - All Cohorts Combined\nLog Odds Ratios (CELIAC vs CONTROL)\nTimepoints in Chronological Order",
               fontsize = 10,
               fontsize_row = 8,
               fontsize_col = 10,
               show_rownames = TRUE,
               border_color = "white",
               cellwidth = 25,
               cellheight = 12,
               treeheight_row = 30,
               treeheight_col = 0,  # Remove column dendrogram
               angle_col = 45)
      
      dev.off()
      cat("Ordered combined PA heatmap created with", nrow(pa_matrix), "genes\n")
    }
  }
}

# Create ordered summary heatmap of significant gene counts
cat("Creating ordered summary heatmap of significant gene counts...\n")

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
  
  # IMPORTANT: Reorder columns to match timepoint order
  available_timepoints <- intersect(timepoint_order, colnames(summary_matrix))
  if(length(available_timepoints) > 0) {
    summary_matrix <- summary_matrix[, available_timepoints, drop = FALSE]
  }
  
  cat("Summary timepoint order:", colnames(summary_matrix), "\n")
  
  # Ordered summary plot
  pdf("../Orf_Contig_Phrog_compositional/figures/significant_genes_summary_ordered_heatmap.pdf", 
      width = 10, height = 6)
  
  pheatmap(summary_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,  # KEEP FALSE to maintain timepoint order
           scale = "none",
           color = colorRampPalette(c("white", "red"))(100),
           main = "Number of Significant Genes by Timepoint and Model\nTimepoints in Chronological Order",
           fontsize = 11,
           fontsize_row = 10,
           fontsize_col = 11,
           show_rownames = TRUE,
           border_color = "grey70",
           display_numbers = TRUE,
           number_format = "%.0f",
           cellwidth = 40,
           cellheight = 20,
           angle_col = 45)
  
  dev.off()
  cat("Ordered summary heatmap created\n")
}

cat("All ordered heatmaps completed!\n")
cat("New files saved with '_ordered' suffix in Orf_Contig_Phrog_compositional/figures/\n")