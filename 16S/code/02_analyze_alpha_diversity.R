#!/usr/bin/env Rscript

# Section: Main function - Alpha diversity analysis... ----
main <- function() {
  # Section: Load data & calculate alpha diversity... ----
  cat("Loading phyloseq objects...\n")
  ps_usa <- readRDS(p_f_ps_usa)
  ps_italy <- readRDS(p_f_ps_italy)
  cat("Rarefying phyloseq objects...\n")
  ps_usa_rare <- rarefy_phyloseq(ps_usa)
  ps_italy_rare <- rarefy_phyloseq(ps_italy)
  cat("Calculating alpha diversity...\n")
  # first prepare_metadata to set numeric variables as numeric 
  # and categorical variables as factors with proper reference levels
  ps_usa_rare_prepared <- prepare_metadata(ps_usa_rare)
  ps_italy_rare_prepared <- prepare_metadata(ps_italy_rare)
  # then calculate alpha diversity
  df_alpha_usa <- calculate_alpha_diversity(ps_usa_rare_prepared)
  df_alpha_italy <- calculate_alpha_diversity(ps_italy_rare_prepared)

  # Section: Fit mixed effects models... ----
  v_metrics <- c("shannon", "simpson", "observed", "chao1")
  df_results <- data.frame()
  cat("Fitting mixed effects models...\n")
  # loop over all metrics to fit models for both cohorts
  for (metric in v_metrics) {
    cat("  Processing", metric, "...\n")
    cat("    USA cohort...\n")
    df_results_usa <- fit_model(df_alpha_usa, metric, cohort_name = "USA")
    cat("    Italy cohort...\n")
    df_results_italy <- fit_model(df_alpha_italy, metric, cohort_name = "Italy")
    df_results <- rbind(df_results, df_results_usa, df_results_italy)
  }
  # Section: Save results to excel and pdf files ... ----
  cat("Save results to excel file...\n")
  save_results_to_excel_file(df_results, v_metrics)

  cat("Creating plots...\n")
  pdf(p_f_figs, width = 10, height = 8)
  cat("  Forest plots...\n")
  plt_forest_combined <- create_forest_plots(df_results, v_metrics)
  print(plt_forest_combined)
  cat("  Box plots...\n")
  create_boxplots(df_results, v_metrics, df_alpha_usa, df_alpha_italy)
  cat("  Trajectory plots...\n")
  create_trajectory_plots(df_results, v_metrics, df_alpha_usa, df_alpha_italy)
  dev.off()
  # End of analysis
  cat("End of alpha diversity analysis...\n")
}


# Section: Define helper functions for this file... ----
# rarefy_phyloseq: rarefy the phyloseq object to the minimum sample depth
rarefy_phyloseq <- function(ps, seed = 42) {
  # Rarefication is a statistic technitque used in ecology and microbiome
  # analysis to standardize the sample sizes for fair comparison of samples.
  set.seed(seed)
  min_depth <- min(sample_sums(ps))
  cat("Rarefying to minimum depth without replacement:", min_depth, "\n")
  ps_rare <- rarefy_even_depth(ps,
    sample.size = min_depth,
    rngseed = seed, replace = FALSE
  )
  return(ps_rare)
}

calculate_alpha_diversity <- function(ps) {
  # calculate_df_alphaersity: calculate the alpha diversity v_metrics
  # Shannon, Simpson, Observed, Chao1 are calculated. Then, diversity
  # v_metrics are merged with metadf_data.
  # Extract OTU table
  otu_df_data <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu_df_data <- t(otu_df_data)
  }
  # Calculate diversity v_metrics: Shannon, Simpson, Observed, Chao1
  df_alpha <- data.frame(
    sample_id = rownames(otu_df_data),
    shannon = vegan::diversity(otu_df_data, index = "shannon"),
    simpson = vegan::diversity(otu_df_data, index = "simpson"),
    observed = vegan::specnumber(otu_df_data),
    chao1 = apply(otu_df_data, 1, function(x) vegan::estimateR(x)["S.chao1"]),
    row.names = rownames(otu_df_data)
  )
  # Get metadf_data
  df_metadata <- as.data.frame(sample_data(ps))
  # Merge alpha diversity with sample df_data
  df_alpha_full <- cbind(df_alpha, df_metadata)
  df_alpha_full
}

prepare_metadata <- function(ps) {
  # prepare_metadata: prepare the metadata for the analysis
  # set numeric variables as numeric and categorical variables 
  # as factors with proper reference levels
  # Get metadata
  df_metadata <- as.data.frame(sample_data(ps))
  df_metadata$time_to_onset <-
    as.numeric(as.character(df_metadata$time_to_onset))
  df_metadata$age_at_gluten_introduction <-
    as.numeric(as.character(df_metadata$age_at_gluten_introduction))
  df_metadata$age_at_solid_introduction <-
    as.numeric(as.character(df_metadata$age_at_solid_introduction))
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

  sample_data(ps) <- sample_data(df_metadata)
  return(ps)
}

# Function to fit mixed effects model for each metric
fit_model <- function(df_data, metric, cohort_name){
  formula_str <- paste0(
    metric, " ~ disease_status * time_to_onset + sex + ",
    "age_at_gluten_introduction + age_at_solid_introduction + ",
    "hla_risk_category + delivery_mode + feeding_type_first_year + ",
    "(1|patient_id)"
  )
  model <- lmerTest::lmer(as.formula(formula_str),
    data = df_data, REML = FALSE
  )
  # Extract coefficients and p-values
  coef_summary <- summary(model)$coefficients
  p_values <- coef_summary[, "Pr(>|t|)"]
  # Create dataframe with results
  df_results <- data.frame(
    cohort = cohort_name,
    metric = metric,
    predictor = rownames(coef_summary),
    estimate = coef_summary[, "Estimate"],
    std_error = coef_summary[, "Std. Error"],
    t_value = coef_summary[, "t value"],
    p_value = p_values,
    stringsAsFactors = FALSE
  )

  # Rename predictors to cleaner names
  df_results$predictor <-
    gsub("disease_statusCELIAC", "celiac", df_results$predictor)
  df_results$predictor <-
    gsub("sexFemale", "female", df_results$predictor)
  df_results$predictor <-
    gsub("hla_risk_categoryHigh Risk", "hla_high_risk", df_results$predictor)
  df_results$predictor <-
    gsub("hla_risk_categoryLow/No Risk", "hla_low_risk", df_results$predictor)
  df_results$predictor <-
    gsub("delivery_modeC-Section", "c_section", df_results$predictor)
  df_results$predictor <-
    gsub("feeding_type_first_yearFormula", "formula", df_results$predictor)
  df_results$predictor <-
    gsub(
      "feeding_type_first_yearBreastmilk_and_formula",
      "breastmilk_and_formula", df_results$predictor
    )
  # FDR correction across predictors within this model
  df_results$p_value_adj <- p.adjust(df_results$p_value, method = "fdr")
  # Helper function to add significance stars
  
  # Add significance asterisks
  df_results$sig_adj <- sapply(df_results$p_value_adj, add_significance)
  # Set rownames as concatenation of metric and predictor
  rownames(df_results) <-
    paste(df_results$cohort, df_results$metric, df_results$predictor, sep = "_")
  df_results
}

save_results_to_excel_file <- function(df_results, v_metrics) {
  # Create workbook
  wb <- openxlsx::createWorkbook()
  # Add a sheet per metric with formatted headers and auto-width columns
  for (metric_name in v_metrics) {
    df_results_metric <- df_results[df_results$metric == metric_name, ]
    openxlsx::addWorksheet(wb, sheetName = metric_name)
    openxlsx::writeData(wb, sheet = metric_name, x = df_results_metric)
    openxlsx::addStyle(
      wb,
      sheet = metric_name,
      style = openxlsx::createStyle(textDecoration = "bold"),
      rows = 1,
      cols = seq_len(ncol(df_results_metric))
    )
    openxlsx::setColWidths(
      wb,
      sheet = metric_name,
      cols = seq_len(ncol(df_results_metric)),
      widths = "auto"
    )
  }
  openxlsx::saveWorkbook(wb, p_f_excel, overwrite = TRUE)
  cat("Excel file saved to:", p_f_excel, "\n")
}
# Function to create forest plots (both combined and individual)
create_forest_plots <- function(df_results, v_metrics, p_dir_fig) {
  # Helper function to create single forest plot
  create_one_forest_plot <- function(df_results, metric_name) {
    # Map metric to a human-friendly panel title
    metric_titles <- c(
      shannon = "Shannon",
      simpson = "Simpson",
      observed = "Observed",
      chao1 = "Chao1"
    )
    panel_title <- metric_titles[[metric_name]]
    # Define manual order of predictors (customize this as needed)
    predictor_order <- c(
      "age_at_solid_introduction", "age_at_gluten_introduction",
      "breastmilk_and_formula", "formula", "c_section", "hla_low_risk",
      "hla_high_risk", "female", "celiac:time_to_onset", "celiac",
      "time_to_onset"
    )

    # Filter for metric and exclude intercept
    df_results_plt <- df_results %>%
      filter(metric == metric_name, predictor != "(Intercept)") %>%
      # Convert predictor to factor with manual ordering
      mutate(predictor = factor(predictor, levels = predictor_order)) %>%
      mutate(predictor_label = PREDICTOR_LABELS[as.character(predictor)])

    # Create plot
    plt_single <- ggplot(
      df_results_plt,
      aes(x = estimate, y = predictor, color = cohort)
    ) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_errorbarh(
        aes(
          xmin = estimate - 1.96 * std_error,
          xmax = estimate + 1.96 * std_error
        ),
        width = 0.5, size = 0.8, position = position_dodge(0.6)
      ) +
      geom_point(size = 2, position = position_dodge(0.6)) +
      scale_color_manual(values = COLORS_COHORT) +
      scale_y_discrete(labels = df_results_plt$predictor_label) +
      labs(
        x = "Effect Size (95% CI)",
        y = "Predictor",
        title = metric_titles[[metric_name]]
      ) +
      theme_bw(base_size = BASE_SIZE) +
      theme(
        axis.text = element_text(size = BASE_SIZE),
        plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5)
      )

    # Add significance asterisks using geom_text for better positioning
    if (any(df_results_plt$sig_adj != "")) {
      sig_data <- df_results_plt[df_results_plt$sig_adj != "", ]

      # Calculate exact y positions for dodged text
      sig_data$y_pos <- as.numeric(sig_data$predictor)
      # Adjust y position based on cohort (USA = -0.2, Italy = +0.2)
      sig_data$y_pos <- sig_data$y_pos + ifelse(sig_data$cohort == "USA", 0.3, -0.3)

      plt_single <- plt_single + geom_text(
        data = sig_data,
        aes(
          x = estimate + 1.96 * std_error + 0.05,
          y = y_pos,
          label = sig_adj,
          color = cohort
        ),
        size = 8,
        inherit.aes = FALSE,
        show.legend = FALSE,
        hjust = 0, # Left-align the text
        vjust = 0.5 # Center vertically
      )
    }
    return(plt_single)
  }
  # Create individual forest plots for each metric
  l_plt <- list()
  for (metric in v_metrics) {
    l_plt[[metric]] <- create_one_forest_plot(df_results, metric)
  }

  # Extract legend from one of the plots
  legend_plot <- l_plt[["shannon"]] + theme(legend.position = "top")
  shared_legend <- cowplot::get_legend(legend_plot)
  
  # Create combined plot using patchwork but remove individual legends
  # Only left column gets y-axis title; only bottom row gets x-axis title
  p_shannon  <- l_plt[["shannon"]] + theme(axis.title.x = element_blank(), legend.position = "none")
  p_simpson  <- l_plt[["simpson"]] + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
  p_observed <- l_plt[["observed"]] + theme(legend.position = "none")
  p_chao1    <- l_plt[["chao1"]] + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

  # Create the 2x2 grid without legend
  panels_grid <- (p_shannon | p_simpson) / (p_observed | p_chao1) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "D")),
      theme = theme(
        plot.tag.position = c(0.05, 0.95),
        plot.tag = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    )
  
  # Create main title
  main_title <- cowplot::ggdraw() + 
    cowplot::draw_label("Alpha Diversity - Forest Plots", fontface = "bold", size = 16)
  
  # Combine: title at top, then legend, then panels
  plt_combined <- cowplot::plot_grid(
    main_title,
    shared_legend, 
    panels_grid,
    ncol = 1,
    rel_heights = c(0.05, 0.05, 0.87)
  )
  
  plt_combined
}

# Function to create boxplot plots for all metrics
create_boxplots <- function(df_results, v_metrics, df_alpha_usa, df_alpha_italy) {
  # Helper function to create single boxplot for a metric and cohort
  create_single_boxplot <- function(df_alpha, metric_name, cohort_name, df_results) {
    # Get adjusted p-value from model results for subtitle
    p_adj <- df_results$p_value_adj[
      df_results$metric == metric_name &
        df_results$cohort == cohort_name &
        df_results$predictor == "celiac"
    ]

    # Reorder disease_status so CELIAC appears first
    df_alpha$disease_status_ordered <- factor(df_alpha$disease_status, levels = c("CELIAC", "CONTROL"))
    
    # Create boxplot
    p <- ggplot(df_alpha, aes_string(x = "disease_status_ordered", y = metric_name, fill = "disease_status_ordered")) +
      geom_boxplot(alpha = BOXPLOT_WIDTH, outlier.alpha = BOXPLOT_OUTLIER_ALPHA) +
      scale_fill_manual(values = COLORS_DISEASE) +
        labs(
          subtitle = paste("Adj p-value :", round(p_adj, 2), " (", cohort_name, ")"),
          y = metric_name,
          x = NULL
        ) +
      theme_bw(base_size = BASE_SIZE) +
      theme(
        plot.title = element_text(size = TITLE_SIZE, hjust = 0.5),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5)
      )

    return(p)
  }

  # Create USA plots
  cat("Creating USA boxplots...\n")
  usa_plots <- list()
  for (metric in v_metrics) {
    usa_plots[[metric]] <- create_single_boxplot(df_alpha_usa, metric, "USA", df_results)
  }

  # Create Italy plots
  cat("Creating Italy boxplots...\n")
  italy_plots <- list()
  for (metric in v_metrics) {
    italy_plots[[metric]] <- create_single_boxplot(df_alpha_italy, metric, "Italy", df_results)
  }

  # Apply axis modifications for the 2x4 grid layout
  # Top row (Shannon): USA has y-axis, Italy has no y-axis, both have no x-ticks
  usa_shannon <- usa_plots[["shannon"]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  italy_shannon <- italy_plots[["shannon"]] + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # Second row (Simpson): USA has y-axis, Italy has no y-axis, both have no x-ticks
  usa_simpson <- usa_plots[["simpson"]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  italy_simpson <- italy_plots[["simpson"]] + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # Third row (Observed): USA has y-axis, Italy has no y-axis, both have no x-ticks
  usa_observed <- usa_plots[["observed"]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  italy_observed <- italy_plots[["observed"]] + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # Bottom row (Chao1): USA has both axes, Italy has no y-axis, both keep x-ticks
  usa_chao1 <- usa_plots[["chao1"]]
  italy_chao1 <- italy_plots[["chao1"]] + theme(axis.title.y = element_blank())

  # Create combined plot with USA and Italy side by side for each metric
  combined_plot <- plot_grid(
    # Shannon row: USA | Italy
    usa_shannon, italy_shannon,
    # Simpson row: USA | Italy  
    usa_simpson, italy_simpson,
    # Observed row: USA | Italy
    usa_observed, italy_observed,
    # Chao1 row: USA | Italy
    usa_chao1, italy_chao1,
    ncol = 2,
    nrow = 4,
    labels = c("A", "B", 
               "C", "D",
               "E", "F", 
               "G", "H"),
    label_size = 12
  )

  # Add main title
  main_title <- ggdraw() +
    draw_label("Alpha Diversity - Boxplots for disease status",
      fontface = "bold", size = 16, x = 0.5, hjust = 0.5
    )

  final_plot <- plot_grid(main_title, combined_plot,
    ncol = 1, rel_heights = c(0.1, 0.9)
  )

  # Print the combined plot
  print(final_plot)
}

create_trajectory_plots <- function(df_results, v_metrics, df_alpha_usa, df_alpha_italy) {
  # Get unique time_to_onset values for each cohort separately
  usa_time_values <- sort(unique(as.numeric(as.character(df_alpha_usa$time_to_onset))))
  usa_time_values <- usa_time_values[!is.na(usa_time_values)]
  
  italy_time_values <- sort(unique(as.numeric(as.character(df_alpha_italy$time_to_onset))))
  italy_time_values <- italy_time_values[!is.na(italy_time_values)]
  
  # Helper function to create single trajectory plot for a metric and cohort
  create_single_trajectory_plot <- function(df_alpha, metric_name, cohort_name, df_results, time_breaks) {
    # Prepare data and ensure time_to_onset is numeric
    df_alpha$time_to_onset <- as.numeric(as.character(df_alpha$time_to_onset))

    # Create trajectory plot
    p <- ggplot(df_alpha, aes_string(x = "time_to_onset", y = metric_name, color = "disease_status")) +
      geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
      geom_smooth(method = TRAJECTORY_SMOOTH_METHOD, se = TRAJECTORY_CONFIDENCE_INTERVAL, size = TRAJECTORY_LINE_SIZE) +
      scale_color_manual(values = COLORS_DISEASE) +
      scale_x_continuous(breaks = time_breaks, limits = c(min(time_breaks), max(time_breaks))) +
      labs(
        x = "Time to onset",
        y = metric_name,
        color = "Disease Status"
      ) +
      theme_bw(base_size = BASE_SIZE) +
      theme(
        plot.title = element_text(size = TITLE_SIZE, hjust = 0.5),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5),
        legend.position = LEGEND_POSITION
      )

    return(p)
  }

  # Create USA plots
  cat("Creating USA trajectory plots...\n")
  usa_plots <- list()
  for (metric in v_metrics) {
    usa_plots[[metric]] <- create_single_trajectory_plot(df_alpha_usa, metric, "USA", df_results, usa_time_values)
  }

  # Create Italy plots
  cat("Creating Italy trajectory plots...\n")
  italy_plots <- list()
  for (metric in v_metrics) {
    italy_plots[[metric]] <- create_single_trajectory_plot(df_alpha_italy, metric, "Italy", df_results, italy_time_values)
  }

  # Apply axis/legend modifications for a 2x4 grid (USA | Italy per row)
  # Top row (Shannon): hide x-axis title on both; hide y-axis title on Italy; no legends
  usa_shannon <- usa_plots[["shannon"]] + theme(legend.position = "none",
                                               axis.title.x = element_blank())
  italy_shannon <- italy_plots[["shannon"]] + theme(legend.position = "none",
                                                   axis.title.y = element_blank(),
                                                   axis.title.x = element_blank())

  # Second row (Simpson)
  usa_simpson <- usa_plots[["simpson"]] + theme(legend.position = "none",
                                               axis.title.x = element_blank())
  italy_simpson <- italy_plots[["simpson"]] + theme(legend.position = "none",
                                                   axis.title.y = element_blank(),
                                                   axis.title.x = element_blank())

  # Third row (Observed)
  usa_observed <- usa_plots[["observed"]] + theme(legend.position = "none",
                                                 axis.title.x = element_blank())
  italy_observed <- italy_plots[["observed"]] + theme(legend.position = "none",
                                                     axis.title.y = element_blank(),
                                                     axis.title.x = element_blank())

  # Bottom row (Chao1): keep x-axis title; hide y on Italy; no legends
  usa_chao1 <- usa_plots[["chao1"]] + theme(legend.position = "none")
  italy_chao1 <- italy_plots[["chao1"]] + theme(legend.position = "none",
                                               axis.title.y = element_blank())

  # Combine into a single 2x4 grid with labels A-H
  combined_plot <- plot_grid(
    usa_shannon, italy_shannon,
    usa_simpson, italy_simpson,
    usa_observed, italy_observed,
    usa_chao1, italy_chao1,
    ncol = 2,
    nrow = 4,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    label_size = 12
  )

  # Add a main title and print
  main_title <- ggdraw() +
    draw_label("Alpha Diversity - Trajector plots (celiac:time_to_onset)",
      fontface = "bold", size = 16, x = 0.5, hjust = 0.5
    )

  final_plot <- plot_grid(main_title, combined_plot,
    ncol = 1, rel_heights = c(0.08, 0.92)
  )

  print(final_plot)
}

setup <- function() {
  # Load required libraries
  library(phyloseq)
  library(vegan)
  library(lmerTest)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(openxlsx)
  
  # Set random seed for reproducibility
  set.seed(42)
  
  # Set up paths
  p_dir_wd <<- Sys.getenv(
    "PROJ_CELIAC_PATH",
    "/Users/dhohaabid/Documents/proj_celiac"
  )
  p_dir_out <<- file.path(p_dir_wd, "output")
  p_dir_fig <<- file.path(p_dir_wd, "figures")
  p_dir_res <<- file.path(p_dir_wd, "results")
  
  # Define input files
  p_f_ps_usa <<- file.path(p_dir_out, "ps_gg2_cleaned_usa.RDS")
  p_f_ps_italy <<- file.path(p_dir_out, "ps_gg2_cleaned_italy.RDS")
  
  # Define output files
  p_f_figs <<- file.path(p_dir_fig, "02_alpha_diversity.pdf")
  p_f_excel <<- file.path(p_dir_res, "02_alpha_diversity.xlsx")
  
  # Source utility functions
  source(file.path(p_dir_wd, "code", "utils", "config.R"))
  source(file.path(p_dir_wd, "code", "utils", "utils.R"))
}


# Section: Execute... ----
cat("
╔══════════════════════════════════════════════════════════════╗
║                    ALPHA DIVERSITY ANALYSIS                  ║
║                        STARTING NOW                          ║
╚══════════════════════════════════════════════════════════════╝
")
setup()
main()