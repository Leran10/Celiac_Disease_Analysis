# Timepoint-specific analysis and compositional analysis
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(pheatmap)
  library(vegan)
})

cat("Starting timepoint-specific and compositional analysis...\n")

# Function to extract timepoint-specific results
extract_timepoint_results <- function(results_file, cohort_name) {
  cat("Extracting timepoint results for", cohort_name, "cohort...\n")
  
  if(!file.exists(results_file)) {
    cat("File not found:", results_file, "\n")
    return(NULL)
  }
  
  results <- read.csv(results_file)
  
  # Extract interaction terms for each timepoint
  timepoints <- c("t0-6", "t0-12", "t0-18", "t0-24", "t0-30", "t0-36", "t0-over42")
  
  timepoint_results <- list()
  
  for(tp in timepoints) {
    interaction_term <- paste0("Dx.StatusCELIAC:onset_timeline_combined", tp)
    tp_results <- results[results$term == interaction_term, ]
    
    if(nrow(tp_results) > 0) {
      tp_results$timepoint <- tp
      tp_results$comparison <- "CELIAC_vs_CONTROL"
      timepoint_results[[tp]] <- tp_results
    }
  }
  
  if(length(timepoint_results) > 0) {
    combined_results <- do.call(rbind, timepoint_results)
    output_file <- paste0("Orf_Contig_Phrog_compositional/results/", cohort_name, "_timepoint_specific_results.csv")
    write.csv(combined_results, output_file, row.names = FALSE)
    cat("Timepoint results saved:", output_file, "\n")
    return(combined_results)
  } else {
    cat("No timepoint results found for", cohort_name, "\n")
    return(NULL)
  }
}

# Extract timepoint results for all cohorts and models
cohorts <- c("total", "US", "Italy")
models <- c("PA_model1_glmmTMB_logit", "abundance_model1_nbinom")

for(cohort in cohorts) {
  for(model in models) {
    file_path <- paste0("Orf_Contig_Phrog_compositional/results/", cohort, "_", model, ".csv")
    extract_timepoint_results(file_path, paste0(cohort, "_", model))
  }
}

# Load cohort data for compositional analysis
load_cohort_data <- function(cohort_name) {
  abundance_file <- paste0("../Contig/data/", cohort_name, "/", cohort_name, ".contig.abundance.clean_0.75_0.03.csv")
  metadata_file <- paste0("../Contig/data/", cohort_name, "/", cohort_name, ".contig.metadata.clean_0.75_0.03.csv")
  
  abundance_data <- read.csv(abundance_file) %>% 
    tibble::column_to_rownames("X")
  colnames(abundance_data) <- gsub("X","",colnames(abundance_data))
  
  metadata_data <- read.csv(metadata_file) %>% 
    tibble::column_to_rownames("X")
  
  # Process metadata
  metadata_data$Dx.Status <- factor(metadata_data$Dx.Status, levels = c("CONTROL","CELIAC"))
  metadata_data$onset_timeline_combined <- factor(metadata_data$onset_timeline_combined, 
                                                 levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))
  
  return(list(abundance = abundance_data, metadata = metadata_data))
}

# Compositional analysis for total cohort
cat("Running compositional analysis for total cohort...\n")
total_data <- load_cohort_data("total")

# Calculate diversity metrics
calculate_diversity <- function(abundance_data) {
  # Convert to relative abundance
  rel_abundance <- sweep(abundance_data, 2, colSums(abundance_data), "/")
  
  # Calculate diversity metrics
  shannon <- diversity(t(rel_abundance), index = "shannon")
  simpson <- diversity(t(rel_abundance), index = "simpson")
  richness <- apply(abundance_data > 0, 2, sum)
  evenness <- shannon / log(richness)
  
  return(data.frame(
    shannon = shannon,
    simpson = simpson,
    richness = richness,
    evenness = evenness
  ))
}

# Calculate diversity for total cohort
diversity_metrics <- calculate_diversity(total_data$abundance)
diversity_data <- data.frame(diversity_metrics, total_data$metadata)

# Save diversity results
write.csv(diversity_data, "Orf_Contig_Phrog_compositional/results/total_diversity_metrics.csv", row.names = TRUE)

# Individual trajectory analysis
calculate_individual_trajectories <- function(abundance_data, metadata) {
  patients <- unique(metadata$patientID)
  trajectory_results <- list()
  
  for(patient in patients) {
    patient_samples <- rownames(metadata)[metadata$patientID == patient]
    if(length(patient_samples) >= 3) {  # At least 3 timepoints
      patient_data <- abundance_data[, patient_samples]
      patient_meta <- metadata[patient_samples, ]
      
      # Calculate slopes for each gene
      slopes <- apply(patient_data, 1, function(gene_counts) {
        timepoints <- as.numeric(gsub("t0-?", "", patient_meta$onset_timeline_combined))
        timepoints[patient_meta$onset_timeline_combined == "t0"] <- 0
        timepoints[patient_meta$onset_timeline_combined == "t0-over42"] <- 48
        
        if(length(unique(timepoints)) >= 2) {
          lm_result <- lm(gene_counts ~ timepoints)
          return(coef(lm_result)[2])  # slope
        } else {
          return(NA)
        }
      })
      
      # Calculate stability (coefficient of variation)
      stability <- apply(patient_data, 1, function(x) sd(x) / mean(x))
      
      trajectory_results[[patient]] <- data.frame(
        patientID = patient,
        gene = names(slopes),
        slope = slopes,
        stability = stability,
        dx_status = unique(patient_meta$Dx.Status),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if(length(trajectory_results) > 0) {
    combined_trajectories <- do.call(rbind, trajectory_results)
    return(combined_trajectories)
  } else {
    return(NULL)
  }
}

# Calculate trajectories for total cohort
trajectories <- calculate_individual_trajectories(total_data$abundance, total_data$metadata)
if(!is.null(trajectories)) {
  write.csv(trajectories, "Orf_Contig_Phrog_compositional/results/total_individual_trajectories.csv", row.names = FALSE)
  cat("Individual trajectories calculated and saved\n")
}

# Create summary statistics
cat("Creating summary statistics...\n")

# Count successful models
summary_stats <- data.frame(
  cohort = character(),
  model = character(),
  n_genes = integer(),
  n_results = integer(),
  stringsAsFactors = FALSE
)

result_files <- list.files("Orf_Contig_Phrog_compositional/results/", pattern = "*.csv", full.names = TRUE)
for(file in result_files) {
  if(grepl("_PA_model|_abundance_model", basename(file))) {
    data <- read.csv(file)
    cohort <- gsub("_.*", "", basename(file))
    model <- gsub(".*_model", "model", gsub(".csv", "", basename(file)))
    n_genes <- length(unique(data$gene))
    n_results <- nrow(data)
    
    summary_stats <- rbind(summary_stats, data.frame(
      cohort = cohort,
      model = model,
      n_genes = n_genes,
      n_results = n_results
    ))
  }
}

write.csv(summary_stats, "Orf_Contig_Phrog_compositional/results/analysis_summary.csv", row.names = FALSE)

cat("Timepoint and compositional analysis completed!\n")
cat("Results saved in Orf_Contig_Phrog_compositional/results/\n")