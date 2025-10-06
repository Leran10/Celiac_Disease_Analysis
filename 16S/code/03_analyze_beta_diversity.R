#!/usr/bin/env Rscript

# Beta Diversity Analysis with Mixed.MDMR
# Using R² as effect size proxy with bar chart visualization

# Section: Main function - Beta diversity analysis... ----
main <- function() {
  cat("Loading phyloseq objects...\n")
  ps_usa <- readRDS(p_f_ps_usa)
  ps_italy <- readRDS(p_f_ps_italy)
  cat("Preparing metadata...\n")
  df_metadata_usa <- prepare_metadata(ps_usa)
  df_metadata_italy <- prepare_metadata(ps_italy)
  cat("Calculating aitchison distance...\n")
  m_dist_usa <- calculate_aitchison_dist(ps_usa)
  m_dist_italy <- calculate_aitchison_dist(ps_italy)
  cat("Fitting mdmr models...\n")
  df_results_usa <- fit_model(m_dist_usa, df_metadata_usa, "USA")
  df_results_italy <- fit_model(m_dist_italy, df_metadata_italy, "Italy")
  cat("Saving results...\n")
  df_results <- rbind(df_results_usa, df_results_italy)
  
  # Combine results


  # # Section: Create & save plots... ----
  pdf(p_f_figs, width = 10, height = 8)
  
  cat("Creating bar plot...\n")
  plt_barplot <- create_barplot(df_results)
  print(plt_barplot)

  cat("Creating combined PCoA plot...\n")
  l_res_pcoa_usa <- perform_pcoa(m_dist_usa, df_metadata_usa)
  l_res_pcoa_italy <- perform_pcoa(m_dist_italy, df_metadata_italy)
  plt_pcoa_combined <- create_pcoa_plots(l_res_pcoa_usa, l_res_pcoa_italy, df_results_usa, df_results_italy)
  print(plt_pcoa_combined)

  cat("Creating combined trajectory plots...\n")
  plt_traject_combined <- create_trajectory_plots(
    l_res_pcoa_usa, 
    l_res_pcoa_italy,
    df_results_usa,
    df_results_italy
  )
  print(plt_traject_combined$pc1_combined)
  print(plt_traject_combined$pc2_combined)
  
  # Calculate and display correlations
  cat("Calculating PC correlations...\n")
  correlations <- calculate_pc_correlations(l_res_pcoa_usa, l_res_pcoa_italy, df_metadata_usa, df_metadata_italy)
  
  cat("\n=== PC1 and PC2 Correlations with Covariates ===\n")
  print(correlations)
  
  # Save correlations to file
  write.csv(correlations, file.path(p_dir_res, "03_pc_correlations.csv"), row.names = FALSE)
  cat("Correlations saved to:", file.path(p_dir_res, "03_pc_correlations.csv"), "\n")
  
  dev.off()
}


# Section: Define helper functions for this file
prepare_metadata <- function(ps) {

  # Get metadata
  df_metadata <- data.frame(sample_data(ps))
  df_metadata$time_to_onset <-
    as.numeric(as.character(df_metadata$time_to_onset))
  df_metadata$age_at_gluten_introduction <-
    as.numeric(as.character(df_metadata$age_at_gluten_introduction))
  df_metadata$age_at_solid_introduction <-
    as.numeric(as.character(df_metadata$age_at_solid_introduction))
  df_metadata$age <-
    as.numeric(as.character(df_metadata$age))
  df_metadata$disease_status <-
    relevel(factor(df_metadata$disease_status), ref = "CONTROL")
  df_metadata$sex <-
    relevel(factor(df_metadata$sex), ref = "Male")
  df_metadata$hla_risk_category <-
    relevel(factor(df_metadata$hla_risk_category), ref = "Standard Risk")
  df_metadata$delivery_mode <-
    relevel(factor(df_metadata$delivery_mode), ref = "Vaginal")
  df_metadata$feeding_type_first_year <-
    relevel(factor(df_metadata$feeding_type_first_year), ref = "Breast_fed")

  return(df_metadata)
}
calculate_aitchison_dist <- function(ps) {
  cat("  Calculating Aitchison distance...\n")
  # Extract OTU table
  m_otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    m_otu <- t(m_otu)
  }
  cat("   Applying CZM zero replacement...\n")
  m_otu_czm <- cmultRepl(
    m_otu,
    method = "CZM", output = "p-counts",
    z.delete = FALSE, z.warning = FALSE
  )
  cat("   Applying CLR transformation...\n")
  m_clr <- clr(m_otu_czm)
  cat("   Calculating Euclidean distance (Aitchison)...\n")
  m_dist <- dist(m_clr, method = "euclidean")
  return(m_dist)
}


fit_model <- function(m_dist, df_metadata, cohort_name) {
  str_formula <- paste0(
    " ~ disease_status * time_to_onset + sex + ",
    "age_at_gluten_introduction + age_at_solid_introduction + ",
    "hla_risk_category + delivery_mode + feeding_type_first_year + ",
    "(1|patient_id)"
  )

  model <- mixed.mdmr(
    fmla = as.formula(str_formula),
    data = df_metadata,
    D = m_dist,
    ncores = 1
  )

  # Extract coefficients and p-values
  # Mixed.MDMR returns test statistics (Pseudo F-statistic) and p-values
  # It doesn't return estimates of coefficients, standard errors, etc
  v_f_ratio <- model$stat
  v_p_value <- model$pv
  # Remove Omnibus and Intercept terms
  v_predictor <-
    names(v_f_ratio)[!names(v_f_ratio) %in% c("Omnibus", "(Intercept)")]

  # Create dataframe with results
  df_results <- data.frame(
    cohort = cohort_name,
    metric = "aitchison_distance",
    predictor = v_predictor,
    F_ratio = as.numeric(v_f_ratio[v_predictor]),
    p_value = as.numeric(v_p_value[v_predictor]),
    stringsAsFactors = FALSE
  )

  # FDR correction across predictors within this model
  df_results$p_value_adj <- p.adjust(df_results$p_value, method = "fdr")

  # Add significance asterisks
  df_results$sig_adj <- sapply(df_results$p_value_adj, add_significance)

  # Set rownames as concatenation of metric and predictor
  rownames(df_results) <-
    paste(df_results$cohort, df_results$metric, df_results$predictor, sep = "_")

  return(df_results)
}


# Function to create bar plot for beta diversity results
create_barplot <- function(df_results) {
  # Define manual order of predictors (customize this as needed)
  predictor_order <- c(
    "age_at_solid_introduction", "age_at_gluten_introduction",
    "feeding_type_first_year", "delivery_mode", "hla_risk_category",
    "sex", "disease_status:time_to_onset", "disease_status", "time_to_onset"
  )

  # Filter and order predictors
  df_results_plt <- df_results %>%
    filter(predictor %in% predictor_order) %>%
    mutate(predictor = factor(predictor, levels = predictor_order)) %>%
    mutate(predictor_label = PREDICTOR_LABELS[as.character(predictor)])

  # Create horizontal bar plot
  plt <- ggplot(df_results_plt, aes(x = F_ratio, y = predictor, fill = cohort)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.6) +
    scale_fill_manual(
      values = COLORS_COHORT,
      breaks = c("USA", "Italy")
    ) + # USA first in legend
    scale_y_discrete(labels = df_results_plt$predictor_label) +
    labs(
      # title = "F-ratios for Aitchison Distance",
      x = "F-ratio",
      y = "Predictor"
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      legend.position = LEGEND_POSITION,
      plot.title = element_text(size = TITLE_SIZE, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5),
      axis.text = element_text(size = AXIS_TEXT_SIZE),
      axis.title = element_text(size = AXIS_TITLE_SIZE),
      legend.title = element_blank(),
      legend.text = element_text(size = LEGEND_TEXT_SIZE)
    )

  # Add significance asterisks based on p-values
  if (any(df_results_plt$sig_adj != "")) {
    sig_data <- df_results_plt[df_results_plt$sig_adj != "", ]

    # Calculate y positions for dodged bars - center in the middle of bars
    sig_data$y_pos <- as.numeric(sig_data$predictor)
    sig_data$y_pos <- sig_data$y_pos + ifelse(sig_data$cohort == "USA", 0.2, -0.2)

    plt <- plt + geom_text(
      data = sig_data,
      aes(
        x = F_ratio + 0.05,
        y = y_pos,
        label = sig_adj,
        color = cohort # Color asterisks to match bar colors
      ),
      size = BETA_ASTERISK_SIZE,
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5
    ) +
      scale_color_manual(
        values = COLORS_COHORT,
        breaks = c("USA", "Italy"), # Match legend order
        guide = "none"
      ) # Hide color legend since we have fill legend
  }

  return(plt)
}

# Function to perform PCoA
perform_pcoa <- function(m_dist, df_metadata) {
  l_res_pcoa <-
    cmdscale(m_dist, k = min(nrow(as.matrix(m_dist)) - 1, 10), eig = TRUE)

  # Calculate variance explained
  v_var_expl <- l_res_pcoa$eig / sum(abs(l_res_pcoa$eig)) * 100

  # Create PCoA data frame
  df_pcoa <- data.frame(
    PC1 = l_res_pcoa$points[, 1],
    PC2 = l_res_pcoa$points[, 2]
  )
  df_pcoa_metadata <- df_pcoa %>%
    rownames_to_column(var = "sample_id") %>%
    left_join(
      df_metadata %>%
        rownames_to_column(var = "sample_id"),
      by = "sample_id"
    )

  return(list(pcoa = df_pcoa_metadata, var_expl = v_var_expl))
}


# Function to create combined PCoA plot for USA and Italy cohorts side by side
create_pcoa_plots <- function(l_res_pcoa_usa, l_res_pcoa_italy, df_results_usa, df_results_italy) {
  # Helper function to create a single PCoA plot
  create_one_pcoa_plot <- function(l_res_pcoa, cohort_name) {
    plt <- ggplot(l_res_pcoa$pcoa, aes(x = PC1, y = PC2, color = disease_status)) +
      geom_point(size = BETA_POINT_SIZE, alpha = BETA_POINT_ALPHA) +
      stat_ellipse(aes(group = disease_status), type = "norm", level = 0.95, size = BETA_ELLIPSE_SIZE) +
      scale_color_manual(values = COLORS_DISEASE) +
      labs(
        title = cohort_name,
        x = paste0("PC1 (", round(l_res_pcoa$var_expl[1], 1), "%)"),
        y = paste0("PC2 (", round(l_res_pcoa$var_expl[2], 1), "%)"),
        color = "Disease Status"
      ) +
      theme_bw(base_size = BASE_SIZE) +
      theme(
        legend.position = LEGEND_POSITION,
        plot.title = element_text(size = TITLE_SIZE, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5)
      )
    
    return(plt)
  }

  plt_usa <- create_one_pcoa_plot(l_res_pcoa_usa, "USA")
  plt_italy <- create_one_pcoa_plot(l_res_pcoa_italy, "Italy")

  # Combine plots side by side using patchwork
  plt_combined <- (plt_usa | plt_italy) +
    plot_layout(guides = "collect") &
    theme(legend.position = LEGEND_POSITION)

  return(plt_combined)
}


# Function to create combined trajectory plots for USA and Italy cohorts side by side
create_trajectory_plots <- function(l_res_pcoa_usa, l_res_pcoa_italy, df_results_usa, df_results_italy) {
  # helper function to create a single trajectory plot
  create_one_trajectory_plot <- function(l_res_pcoa, pc, cohort_name) {
    # Format p-values
    format_p <- function(p) {
      if (p < 0.001) return("< 0.001")
      if (p < 0.01) return(sprintf("%.3f", p))
      return(sprintf("%.3f", p))
    }
    
    if (pc == "PC1") {
      plt <- ggplot(l_res_pcoa$pcoa, aes(x = time_to_onset, y = PC1, color = disease_status))
      y_label <- paste0("PC1 (", round(l_res_pcoa$var_expl[1], 1), "%)")
    } else if (pc == "PC2") {
      plt <- ggplot(l_res_pcoa$pcoa, aes(x = time_to_onset, y = PC2, color = disease_status))
      y_label <- paste0("PC2 (", round(l_res_pcoa$var_expl[2], 1), "%)")
    } else {
      stop("Invalid PC number")
    }
    
    plt <- plt + 
      geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
      geom_smooth(
        method = TRAJECTORY_SMOOTH_METHOD, 
        formula = y ~ x,
        se = TRAJECTORY_CONFIDENCE_INTERVAL,
        size = TRAJECTORY_LINE_SIZE
      ) +
      scale_color_manual(values = COLORS_DISEASE) +
      labs(
        title = cohort_name,
        x = "Time to Onset (months)",
        y = y_label
      ) +
      scale_x_continuous(breaks = sort(unique(l_res_pcoa$pcoa$time_to_onset))) +
      theme_bw(base_size = BASE_SIZE) +
      theme(
        legend.position = LEGEND_POSITION,
        plot.title = element_text(size = TITLE_SIZE, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5)
      )
    
    return(plt)
  }
  # Create PC1 trajectory plots
  plt_pc1_usa <- create_one_trajectory_plot(l_res_pcoa_usa, "PC1", "USA")
  plt_pc1_italy <- create_one_trajectory_plot(l_res_pcoa_italy, "PC1", "Italy")

  # Create PC2 trajectory plots
  plt_pc2_usa <- create_one_trajectory_plot(l_res_pcoa_usa, "PC2", "USA")
  plt_pc2_italy <- create_one_trajectory_plot(l_res_pcoa_italy, "PC2", "Italy")
  
  # Combine PC1 plots side by side
  plt_pc1_combined <- (plt_pc1_usa | plt_pc1_italy) +
    plot_layout(guides = "collect") &
    theme(legend.position = LEGEND_POSITION)

  # Combine PC2 plots side by side
  plt_pc2_combined <- (plt_pc2_usa | plt_pc2_italy) +
    plot_layout(guides = "collect") &
    theme(legend.position = LEGEND_POSITION)

  return(list(pc1_combined = plt_pc1_combined, pc2_combined = plt_pc2_combined))
}

# Function to calculate correlations between covariates and PC1/PC2
calculate_pc_correlations <- function(l_res_pcoa_usa, l_res_pcoa_italy, df_metadata_usa, df_metadata_italy) {
  # Function to calculate correlations for a single cohort
  calculate_cohort_correlations <- function(l_res_pcoa, df_metadata, cohort_name) {
    # Get the merged data
    df_data <- l_res_pcoa$pcoa
    
    # Select numeric covariates for correlation
    numeric_vars <- c("time_to_onset", "age_at_gluten_introduction", 
                     "age_at_solid_introduction", "age")
    
    # Check which variables exist in the data
    available_vars <- numeric_vars[numeric_vars %in% names(df_data)]
    
    # Calculate correlations with PC1 and PC2
    correlations <- data.frame(
      cohort = cohort_name,
      variable = character(),
      PC1_correlation = numeric(),
      PC1_p_value = numeric(),
      PC2_correlation = numeric(),
      PC2_p_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (var in available_vars) {
      # Remove missing values
      complete_data <- df_data[complete.cases(df_data[c("PC1", "PC2", var)]), ]
      
      if (nrow(complete_data) > 3) {  # Need at least 4 observations for correlation
        # PC1 correlation
        pc1_cor <- cor.test(complete_data$PC1, complete_data[[var]])
        # PC2 correlation  
        pc2_cor <- cor.test(complete_data$PC2, complete_data[[var]])
        
        # Add to results
        correlations <- rbind(correlations, data.frame(
          cohort = cohort_name,
          variable = var,
          PC1_correlation = pc1_cor$estimate,
          PC1_p_value = pc1_cor$p.value,
          PC2_correlation = pc2_cor$estimate,
          PC2_p_value = pc2_cor$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    return(correlations)
  }
  
  # Calculate correlations for both cohorts
  cor_usa <- calculate_cohort_correlations(l_res_pcoa_usa, df_metadata_usa, "USA")
  cor_italy <- calculate_cohort_correlations(l_res_pcoa_italy, df_metadata_italy, "Italy")
  
  # Combine results
  all_correlations <- rbind(cor_usa, cor_italy)
  
  # Add significance indicators
  all_correlations$PC1_significant <- all_correlations$PC1_p_value < 0.05
  all_correlations$PC2_significant <- all_correlations$PC2_p_value < 0.05
  
  # Round values for display
  all_correlations$PC1_correlation <- round(all_correlations$PC1_correlation, 3)
  all_correlations$PC1_p_value <- round(all_correlations$PC1_p_value, 4)
  all_correlations$PC2_correlation <- round(all_correlations$PC2_correlation, 3)
  all_correlations$PC2_p_value <- round(all_correlations$PC2_p_value, 4)
  
  return(all_correlations)
}


# Section: Setup function... ----
setup <- function() {
  # Load required libraries
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(compositions)
  library(zCompositions) # For CZM zero replacement
  library(patchwork)
  library(MDMR)
  library(gridExtra)

  # Set random seed for reproducibility
  set.seed(42)

  # Set up paths
  p_dir_wd <<- Sys.getenv(
    "PROJ_CELIAC_PATH",
    "/Users/dabid/Documents/proj_celiac"
  )
  p_dir_out <<- file.path(p_dir_wd, "output")
  p_dir_fig <<- file.path(p_dir_wd, "figures")
  p_dir_res <<- file.path(p_dir_wd, "results")

  # Define input files
  p_f_ps_usa <<- file.path(p_dir_out, "ps_gg2_cleaned_usa.RDS")
  p_f_ps_italy <<- file.path(p_dir_out, "ps_gg2_cleaned_italy.RDS")

  # Define output files
  p_f_figs <<- file.path(p_dir_fig, "03_beta_diversity.pdf")
  p_f_excel <<- file.path(p_dir_res, "03_beta_diversity.xlsx")

  # Source utility functions
  source(file.path(p_dir_wd, "code", "utils", "config.R"))
  source(file.path(p_dir_wd, "code", "utils", "utils.R"))
}


# Section: Execute... ----
cat("
╔══════════════════════════════════════════════════════════════╗
║                    BETA DIVERSITY ANALYSIS                  ║
║                        STARTING NOW                          ║
╚══════════════════════════════════════════════════════════════╝
")
setup()
main()
