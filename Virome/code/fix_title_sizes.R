# Fix heatmap title sizing and formatting
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(viridis)
  library(RColorBrewer)
  library(reshape2)
})

cat("Creating heatmaps with properly sized titles...\n")

# Define correct timepoint order
timepoint_order <- c("t0", "t0-6", "t0-12", "t0-18", "t0-24", "t0-30", "t0-36", "t0-over42")

# Function to create heatmap with proper title sizing
create_fixed_title_heatmap <- function(timepoint_file, output_name, title_prefix) {
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
  
  # Reorder columns to match timepoint order
  available_timepoints <- intersect(timepoint_order, colnames(estimate_matrix))
  if(length(available_timepoints) > 0) {
    estimate_matrix <- estimate_matrix[, available_timepoints, drop = FALSE]
  }
  
  cat("Creating heatmap with", nrow(estimate_matrix), "genes and", ncol(estimate_matrix), "timepoints\n")
  
  # Calculate dimensions
  n_genes <- nrow(estimate_matrix)
  n_timepoints <- ncol(estimate_matrix)
  
  # Adjust figure dimensions to accommodate title
  cell_height <- 0.3
  cell_width <- 0.8
  
  fig_height <- max(4, min(12, n_genes * cell_height + 3))  # Extra space for title
  fig_width <- max(6, n_timepoints * cell_width + 3)
  
  # Create shortened, multi-line title
  model_type <- ifelse(grepl("PA", output_name), "PA Model", "Abundance Model")
  cohort <- toupper(gsub("_.*", "", output_name))
  
  # Shorter title with line breaks
  main_title <- paste0(cohort, " ", model_type, "\nCELIAC vs CONTROL")
  
  # Create heatmap with smaller title font
  pdf(paste0("../Orf_Contig_Phrog_compositional/figures/", output_name, "_fixed_title.pdf"), 
      width = fig_width, height = fig_height)
  
  pheatmap(estimate_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "none",
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-max(abs(estimate_matrix)), max(abs(estimate_matrix)), length.out = 101),
           main = main_title,
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 10,
           fontsize_main = 12,     # Smaller main title font
           show_rownames = TRUE,
           border_color = "white",
           cellwidth = 25,
           cellheight = 12,
           treeheight_row = 30,
           treeheight_col = 0,
           angle_col = 45)
  
  dev.off()
  
  return(estimate_matrix)
}

# Create fixed title heatmaps for all timepoint-specific results
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
  
  create_fixed_title_heatmap(file, base_name, title_prefix)
}

# Create fixed title combined heatmap
cat("Creating combined PA heatmap with fixed title...\n")

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
      # Reorder columns
      available_timepoints <- intersect(timepoint_order, colnames(pa_matrix))
      if(length(available_timepoints) > 0) {
        pa_matrix <- pa_matrix[, available_timepoints, drop = FALSE]
      }
      
      # Dimensions with extra space for title
      n_genes <- nrow(pa_matrix)
      n_timepoints <- ncol(pa_matrix)
      fig_height <- max(5, min(10, n_genes * 0.25 + 3))  # Extra space
      fig_width <- max(7, n_timepoints * 0.8 + 3)
      
      # Shorter title
      combined_title <- "Combined PA Models\nAll Cohorts - CELIAC vs CONTROL"
      
      pdf("../Orf_Contig_Phrog_compositional/figures/combined_PA_all_cohorts_fixed_title.pdf", 
          width = fig_width, height = fig_height)
      
      pheatmap(pa_matrix,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               breaks = seq(-max(abs(pa_matrix)), max(abs(pa_matrix)), length.out = 101),
               main = combined_title,
               fontsize = 10,
               fontsize_row = 8,
               fontsize_col = 10,
               fontsize_main = 12,     # Smaller main title font
               show_rownames = TRUE,
               border_color = "white",
               cellwidth = 25,
               cellheight = 12,
               treeheight_row = 30,
               treeheight_col = 0,
               angle_col = 45)
      
      dev.off()
      cat("Combined PA heatmap with fixed title created\n")
    }
  }
}

# Create fixed title summary heatmap
cat("Creating summary heatmap with fixed title...\n")

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
  
  # Reorder columns
  available_timepoints <- intersect(timepoint_order, colnames(summary_matrix))
  if(length(available_timepoints) > 0) {
    summary_matrix <- summary_matrix[, available_timepoints, drop = FALSE]
  }
  
  # Shorter title
  summary_title <- "Significant Gene Counts\nBy Timepoint and Model"
  
  # Fixed summary plot with extra height for title
  pdf("../Orf_Contig_Phrog_compositional/figures/significant_genes_summary_fixed_title.pdf", 
      width = 10, height = 7)  # Increased height
  
  pheatmap(summary_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "none",
           color = colorRampPalette(c("white", "red"))(100),
           main = summary_title,
           fontsize = 11,
           fontsize_row = 10,
           fontsize_col = 11,
           fontsize_main = 14,     # Smaller main title font
           show_rownames = TRUE,
           border_color = "grey70",
           display_numbers = TRUE,
           number_format = "%.0f",
           cellwidth = 40,
           cellheight = 20,
           angle_col = 45)
  
  dev.off()
  cat("Summary heatmap with fixed title created\n")
}

cat("All heatmaps with fixed titles completed!\n")
cat("New files saved with '_fixed_title' suffix\n")