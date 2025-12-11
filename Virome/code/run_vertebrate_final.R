## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)





## -----------------------------------------------------------------------------

filter_and_save_wide_data <- function(genesCoverm, 
                                      covered_fraction_threshold = 0.75, 
                                      output_path = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs") {
  
  # Set column names to match the original structure
  colnames(genesCoverm) <- c("sampleID", "ORFID", "RPKM", "ReadCount", "Variance", "Mean", "covered_fraction", "covered_bases")
  
  # Apply coverage filter
  genesCoverm2 <- genesCoverm %>%
    mutate(
      ReadCount_modified = ifelse(covered_fraction >= covered_fraction_threshold, ReadCount, 0),
      Changed = ifelse(ReadCount != ReadCount_modified, 1, 0)
    ) %>%
    mutate(ReadCount = ReadCount_modified) %>%
    select(-ReadCount_modified)
  
  # Count the number of changes
  num_changes <- sum(genesCoverm2$Changed)
  print(paste("Number of changes:", num_changes))
  
  # Convert to wide format
  wide_rpkm_genes <- reshape2::dcast(genesCoverm2, ORFID ~ sampleID, value.var = "ReadCount")
  colnames(wide_rpkm_genes) <- gsub("_stats", "", colnames(wide_rpkm_genes))
  
  # Filter out rows with all zeroes
  filtered_wide_rpkm_genes <- wide_rpkm_genes %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  
  # Save the filtered wide format data
  rds_filename_wide <- paste0("wide_rpkm_genes_", covered_fraction_threshold * 100, "Cov.rds")
  
  saveRDS(filtered_wide_rpkm_genes, file = file.path(output_path, rds_filename_wide))
  
  return(filtered_wide_rpkm_genes)
}





# ==========================
# Auto Heatmaps (Complete, tailored to your labels)
# ==========================
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tibble)
library(data.table)
library(edgeR)
library(glmmTMB)
library(emmeans)
library(purrr)
library(broom.mixed)


# ---- Catalog of timepoint-based families matching YOUR labels ----
.family_catalog <- tibble::tribble(
  ~family,                    ~match_regex,                                            ~tp_regex,                                                                 ~note,
  "TP (CONTROL)",             "^onset_timeline_combined",                              "onset_timeline_combined(t0(?:-(?:\\d+|over\\d+))?)$",                     "main TP effects in CONTROL",
  "TP (CELIAC)",  "^TP \\(CELIAC\\):",  "TP \\(CELIAC\\):\\s*(t0(?:-(?:\\d+|over\\d+))?)\\s+vs\\s+t0",  "CELIAC vs t0",
  "Between groups",           "^CONTROL vs CELIAC @\\s*",                              "^CONTROL vs CELIAC @\\s*(t0(?:-(?:\\d+|over\\d+))?)$",                     "between-group at each timepoint",
  "Dx × TP interactions",     "^Dx\\.Status\\w+:onset_timeline_combined",              "onset_timeline_combined(t0(?:-(?:\\d+|over\\d+))?)$",                     "interaction rows carrying a TP token"
  # If you later add adjacent labels, you can re-enable/add patterns here.
)

# ---- Detect which families exist in df$term ----
detect_timepoint_families <- function(df) {
  stopifnot("term" %in% names(df))
  .family_catalog %>%
    rowwise() %>%
    mutate(has = any(grepl(match_regex, df$term))) %>%
    ungroup() %>%
    filter(has) %>%
    pull(family)
}



# ---- Build a single family's heatmap (significant cells only) ----
make_family_heatmap <- function(
    df, family, tp_levels,
    fdr = 0.05,
    scale = c("log", "natural"),
    min_contigs = 1,
    analysis_type = c("PA", "abundance")
  ) {
    scale <- match.arg(scale)
    analysis_type <- match.arg(analysis_type)
    stopifnot(all(c("vertebrate_contig_id","term","estimate","p.adj") %in% names(df)))

    spec <- .family_catalog %>% filter(family == !!family)
    if (nrow(spec) == 0) stop("Unknown family: ", family)

    # 1) Filter rows for this family
    sel <- df %>% filter(str_detect(term, spec$match_regex[1]))
    if (!nrow(sel)) {
      return(
        ggplot() + theme_void() +
          annotate("text", x = .5, y = .5,
                   label = paste0("term family '", family, "' doesn't have any rows"),
                   size = 5) + xlim(0,1) + ylim(0,1)
      )
    }

    # 2) Extract timepoint token using family-specific regex
    tp_mat <- stringr::str_match(sel$term, spec$tp_regex[1])
    # For families with one captured TP group, it's in column 2
    if (ncol(tp_mat) >= 2) {
      sel$TP <- tp_mat[,2]
    } else {
      # Fallback for CONTROL: sometimes simpler to strip prefix if regex failed
      if (family == "TP (CONTROL)") {
        sel$TP <- sub("^onset_timeline_combined", "", sel$term)
      } else {
        sel$TP <- NA_character_
      }
    }

    # 3) Keep significant cells only - FIXED LINE BELOW
    sel <- sel %>%
      mutate(
        TP        = factor(TP, levels = tp_levels),
        sig       = !is.na(p.adj) & p.adj < fdr,
        value_sig = if (scale == "natural") exp(estimate) else estimate  # FIXED: regular if/else
      ) %>%
      filter(sig, !is.na(TP))

    if (!nrow(sel)) {
      return(
        ggplot() + theme_void() +
          annotate("text", x = .5, y = .5,
                   label = paste0("term family '", family, "' doesn't have any significant contigs"),
                   size = 5) + xlim(0,1) + ylim(0,1)
      )
    }

    # 4) Keep contigs with ≥ min_contigs significant cells in this family
    keep_ids <- sel %>% count(vertebrate_contig_id, name = "hits") %>% filter(hits >= min_contigs) %>% pull(vertebrate_contig_id)
    sel <- sel %>% filter(vertebrate_contig_id %in% keep_ids)
    if (!nrow(sel)) {
      return(
        ggplot() + theme_void() +
          annotate("text", x = .5, y = .5,
                   label = paste0("term family '", family, "': no contigs pass min_contigs = ", min_contigs),
                   size = 5) + xlim(0,1) + ylim(0,1)
      )
    }

    # 5) Order contigs by max |effect|
    ord <- sel %>%
      group_by(vertebrate_contig_id) %>%
      summarise(max_abs = max(abs(value_sig), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_abs)) %>%
      pull(vertebrate_contig_id)

    sel$contig_id <- factor(sel$vertebrate_contig_id, levels = ord)
    sel$TP        <- droplevels(sel$TP)

    # 6) Plot - ANALYSIS-SPECIFIC LABELS
    if (analysis_type == "PA") {
      # For presence/absence (logistic regression)
      fill_scale <- if (scale == "log") {
        scale_fill_gradient2(
          name = "log-odds",
          low = "#2b8cbe",
          mid = "white",
          high = "#d7301f",
          midpoint = 0,
          limits = c(-700, 500), # comment this out when making heatmaps for multiple significant cells
          oob = scales::squish # comment this out when making heatmaps for multiple significant cells
        )
      } else {
        scale_fill_gradient2(
          name = "OR",
          low = "#2b8cbe",
          mid = "white",
          high = "#d7301f",
          midpoint = 1,
          trans = "log10",
          limits = c(exp(-700), exp(500)), # comment this out when making heatmaps for multiple significant cells
          oob = scales::squish # comment this out when making heatmaps for multiple significant cells
        )
      }
    } else {
      # For abundance (negative binomial regression)
      fill_scale <- if (scale == "log") {
        scale_fill_gradient2(
          name = "log fold-change",
          low = "#2b8cbe",
          mid = "white",
          high = "#d7301f",
          midpoint = 0,
          limits = c(-700, 500), # comment this out when making heatmaps for multiple significant cells
          oob = scales::squish # comment this out when making heatmaps for multiple significant cells
        )
      } else {
        scale_fill_gradient2(
          name = "Fold-change",
          low = "#2b8cbe",
          mid = "white",
          high = "#d7301f",
          midpoint = 1,
          trans = "log10",
          limits = c(exp(-700), exp(500)), # comment this out when making heatmaps for multiple significant cells
          oob = scales::squish # comment this out when making heatmaps for multiple significant cells
        )
      }
    }

    # Determine subtitle text based on analysis type
    scale_text <- if (analysis_type == "PA") {
      if (scale == "log") "log-odds" else "OR"
    } else {
      if (scale == "log") "log fold-change" else "fold-change"
    }

    ggplot(sel, aes(x = TP, y = vertebrate_contig_id, fill = value_sig)) +
      geom_tile(
        width = 0.95,
        height = 0.95,
        colour = "black",
        linewidth = 0.3
      ) +
      fill_scale +
      scale_x_discrete(drop = TRUE) +
      scale_y_discrete(drop = TRUE) +
      labs(
        title    = paste0("Heatmap — ", family),
        subtitle = paste0("Only significant cells (BH FDR < ", fdr, ")"),
        x = "Timepoint",
        y = "contig"
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.position = "right"
      )
  }


auto_heatmaps_detected <- function(
    df, tp_levels,
    fdr = 0.05,
    scale = c("log", "natural"),
    min_contigs = 1,
    analysis_type = c("PA", "abundance")
  ) {
    scale <- match.arg(scale)
    analysis_type <- match.arg(analysis_type)
    fams  <- detect_timepoint_families(df)
    if (length(fams) == 0) {
      message("No timepoint-based families detected in `term`.")
      return(list())
    }
    out <- setNames(vector("list", length(fams)), fams)
    for (fam in fams) {
      out[[fam]] <- make_family_heatmap(
        df = df,
        family = fam,
        tp_levels = tp_levels,
        fdr = fdr,
        scale = scale,
        min_contigs = min_contigs,
        analysis_type = analysis_type
      )
    }
    out
  }




# Extract timepoint-specific interaction rows and add per-gene FDR
extract_timepoint_results <- function(model_results,
                                      model_name,
                                      gene_col = "vertebrate_contig_id",
                                      adjust_method = "BH") {
  if (is.null(model_results) || nrow(model_results) == 0) {
    return(data.frame())
  }

  timepoint_results <- model_results %>%
    # keep interaction terms of interest
    filter(grepl("Dx.StatusCELIAC:onset_timeline_combined", term)) %>%
    mutate(
      timepoint = gsub("Dx.StatusCELIAC:onset_timeline_combined", "", term),
      model     = model_name,
      # safer than get(): pull from the data frame by name
      gene      = .data[[gene_col]]
    ) %>%
    select(gene, timepoint, estimate, std.error, p.value, conf.low, conf.high, model) %>%
    # per gene x model, adjust p across that gene's timepoints
    group_by(gene, model) %>%
    mutate(
      p.adj = {
        pv <- p.value
        ok <- !is.na(pv)
        out <- rep(NA_real_, length(pv))
        if (any(ok)) out[ok] <- p.adjust(pv[ok], method = adjust_method)
        out
      },
      n_tests_gene = n()  # handy to audit how many timepoints were tested for this gene
    ) %>%
    ungroup()

  timepoint_results
}




# Apply to all result objects you keep in the workspace for a cohort
extract_all_timepoints <- function(cohort_name, adjust_method = "BH") {
  results_list <- list()

  # PA models
  try({
    obj <- paste0(cohort_name, "_glmmTMB_logit_results")
    if (exists(obj)) {
      results_list[["PA_model1"]] <- extract_timepoint_results(
        get(obj), "PA_glmmTMB_logit", gene_col = "vertebrate_contig_id", adjust_method = adjust_method
      )
    }
  }, silent = TRUE)

  # Abundance models
  try({
    obj <- paste0(cohort_name, "_nb_results")
    if (exists(obj)) {
      results_list[["abundance_model1"]] <- extract_timepoint_results(
        get(obj), "abundance_nbinom", gene_col = "vertebrate_contig_id", adjust_method = adjust_method
      )
    }
  }, silent = TRUE)

  if (length(results_list) > 0) {
    do.call(rbind, results_list)
  } else {
    data.frame()
  }
}




## -----------------------------------------------------------------------------

metadata <- read.csv("~/Handley Lab Dropbox/16S/Celiac/metadata/Updated_Metadata_with_Onset_Timeline.csv") %>%
            column_to_rownames("X")




## -----------------------------------------------------------------------------


# # load the raw count table
# contig_count_table <- fread("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/")
# 
# # load contig annotatin table
# contig_annotation_table <- fread("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/contigAnnotations.tsv")
# 
# # identify phage families
# phage_families <- c(
#   "Siphoviridae",
#   "Myoviridae",
#   "Podoviridae",
#   "Herelleviridae",
#   "Drexlerviridae",
#   "Ackermannviridae",
#   "Autographiviridae",
#   "unclassified Caudovirales family",
#   "Microviridae",
#   "Inoviridae"
# )
# 
# # extract vertebrate contig IDs
# vertebrate_contig_IDs <- contig_annotation_table %>%
#                           filter(!family %in% phage_families) %>%
#                           pull(contigID)
# 
# 
# # filter contig count table to only keep vertebrate contigs
# vertebrate_contig_count_table <- contig_count_table %>%
#                                  filter(V2 %in% vertebrate_contig_IDs)
# 
# 
# #write.table(vertebrate_contig_count_table,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/Coverm_concatenated_vertebrate_contigs.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 
# # filter it by 75% coverage rate
# filter_and_save_wide_data(vertebrate_contig_count_table, 
#                           covered_fraction_threshold = 0.75, 
#                           output_path = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total") 



## -----------------------------------------------------------------------------

# load the 75% mapping rate filtered contig table
vertebrate_contig_processed <- readRDS("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/wide_rpkm_genes_75Cov.rds") %>%
                 column_to_rownames("ORFID")


# then let's check colSums and ensure there is no sample have all 0s across orfs
sample_all_0 <- names(colSums(vertebrate_contig_processed)[colSums(vertebrate_contig_processed) == 0])
# "NovaSeq_N966_I12886_39602_Celiac_Leonard_Stool_01_GEMM_119_30M"       "NovaSeq_N966_I12920_39636_Celiac_Leonard_Stool_01_GEMM_154_30M"      
# "NovaSeq_N978_I13149_39572_Celiac_Leonard_Stool_01_GEMM_093_24M"       "NovaSeq_N978_I13150_39573_Celiac_Leonard_Stool_01_GEMM_093_30M"      
# "NovaSeq_N978_I13151_39574_Celiac_Leonard_Stool_01_GEMM_093_36M"       "NovaSeq_N978_I13152_39575_Celiac_Leonard_Stool_01_GEMM_093_48M"      
# "NovaSeq_N978_I13166_39589_Celiac_Leonard_Stool_01_GEMM_106_24M"       "NovaSeq_N980_I13250_39529_Leonard_Human_Stool_Celiac_01_GEMM_057_60M"
# "NovaSeq_N981_I13300_39814_Celiac_Leonard_Stool_02_GEMM_021_36M"       "NovaSeq_N981_I13317_39831_Celiac_Leonard_Stool_02_GEMM_057_24M"      
# "NovaSeq_N982_I13323_39837_Celiac_Leonard_Stool_02_GEMM_077_12M"       "NovaSeq_N983_I13413_39587_Celiac_Leonard_Stool_01_GEMM_106_12M"      
# "NovaSeq_N983_I13417_39799_Leonard_Human_stool_Celiac_02_GEMM_006_60M"


# let's remove the 6 0 counts samples
vertebrate_contig_processed <- vertebrate_contig_processed %>%
                               select(-all_of(sample_all_0))
dim(vertebrate_contig_processed)
# 31264   327


# assign vertebrate_contig_processed to vertebrate_contig_count, so vertebrate_contig_processed won't be changed in downstream analysis
vertebrate_contig_count <- vertebrate_contig_processed

colnames(vertebrate_contig_count) <- gsub("NovaSeq_N966_|NovaSeq_N978_|NovaSeq_N979_|NovaSeq_N980_|NovaSeq_N981_|NovaSeq_N982_|NovaSeq_N983_|_stats","",colnames(vertebrate_contig_count))

colnames(vertebrate_contig_count) <- str_replace_all(colnames(vertebrate_contig_count),"Leonard_Human_stool_Celiac","Celiac_Leonard_Stool")
colnames(vertebrate_contig_count) <- sapply(strsplit(colnames(vertebrate_contig_count),"_Stool_|_Stool_Celiac_"),"[",2)
colnames(vertebrate_contig_count) <- ifelse(colnames(vertebrate_contig_count) %like% "_6Y",gsub("_6Y","_72M",colnames(vertebrate_contig_count)),
                                           ifelse(colnames(vertebrate_contig_count) %like% "_7Y",gsub("_7Y","_84M",colnames(vertebrate_contig_count)),colnames(vertebrate_contig_count)))

vertebrate_contig_count.table_0.75 <- t(vertebrate_contig_count) 
dim(vertebrate_contig_count.table_0.75)
# 327 31264


vertebrate_contig_count.table_0.75 <- merge(metadata,vertebrate_contig_count.table_0.75,by = 0) %>% # 315 samples left
                                      filter(feeding_first_year != "Unknown") %>%  # 299 samples left
                                      filter(Age.at.Gluten.Introduction..months. != "Unknown") # 295 sample left
dim(vertebrate_contig_count.table_0.75)
# 295 samples x  31291 metadata+contigs
# 31237 contigs


all_0_contigs <- names(colSums(vertebrate_contig_count.table_0.75[28:31291])[colSums(vertebrate_contig_count.table_0.75[28:31291]) == 0])
# 975 contigs become 0 through all samples after the 32 samples were removed from the above step
# so we need to remove those contigs

vertebrate_contig_count.table_0.75 <- vertebrate_contig_count.table_0.75 %>%
                            select(-all_of(all_0_contigs))
dim(vertebrate_contig_count.table_0.75)
# 295 samples x 30316 metadata+orfs
# 30289 contigs


vertebrate_contig_abundance.table <- vertebrate_contig_count.table_0.75 %>%
                       select(c("Row.names",28:30316)) %>%
                       column_to_rownames("Row.names") %>% 
                       t(.) %>%
                       data.frame(.)
colnames(vertebrate_contig_abundance.table) <- gsub("X","",colnames(vertebrate_contig_abundance.table))


# check prevelence filtering
prevalence <- rowSums(vertebrate_contig_abundance.table > 0) / ncol(vertebrate_contig_abundance.table)
table(prevalence > 0.03)  # 3% threshold
# FALSE  TRUE 
# 26778  3511 

# Filter Contigs with prevalence > 3%
vertebrate_contig_abundance.clean <- vertebrate_contig_abundance.table[prevalence > 0.03, ]
      
dim(vertebrate_contig_abundance.clean)
# 3511  295 (3% prev)

# check if there are any samples become 0 counts
names(colSums(vertebrate_contig_abundance.clean)[colSums(vertebrate_contig_abundance.clean) == 0])
# "01_GEMM_006_60M" "01_GEMM_060_18M" "01_GEMM_102_12M" "01_GEMM_109_48M" "02_GEMM_006_30M" "02_GEMM_034_12M"

# remove the two 0 counts samples
# we will have to remove sample 02_GEMM_021_36M because it becomes 0 counts at this step
vertebrate_contig_abundance.clean <- vertebrate_contig_abundance.clean %>%
                       select(-c("01_GEMM_006_60M","01_GEMM_060_18M","01_GEMM_102_12M","01_GEMM_109_48M","02_GEMM_006_30M","02_GEMM_034_12M"))
dim(vertebrate_contig_abundance.clean)
# 3511  289

total_vertebrate_contig_abundance_0.75_0.03 <- vertebrate_contig_abundance.clean

#write.csv(vertebrate_contig_abundance.clean,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_abundance_0.75_0.03.csv",quote = FALSE,row.names = TRUE)


metadata.clean <- vertebrate_contig_count.table_0.75[1:27] %>%
                  filter(!Row.names %in% c("01_GEMM_006_60M","01_GEMM_060_18M","01_GEMM_102_12M","01_GEMM_109_48M","02_GEMM_006_30M","02_GEMM_034_12M")) %>%
                  column_to_rownames("Row.names")
dim(metadata.clean)
# 289  26

total_vertebrate_metadata_0.75_0.03 <- metadata.clean

#write.csv(metadata.clean,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_metadata_0.75_0.03.csv",quote = FALSE,row.names = TRUE)


# set up categorical variable levels
total_vertebrate_metadata_0.75_0.03$feeding_first_year <- factor(total_vertebrate_metadata_0.75_0.03$feeding_first_year,levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
total_vertebrate_metadata_0.75_0.03$HLA.Category <- factor(total_vertebrate_metadata_0.75_0.03$HLA.Category,levels = c("Standard Risk","High Risk","Low/No Risk"))
total_vertebrate_metadata_0.75_0.03$Sex <- factor(total_vertebrate_metadata_0.75_0.03$Sex,levels = c("Female","Male"))
total_vertebrate_metadata_0.75_0.03$Delivery.Mode <- factor(total_vertebrate_metadata_0.75_0.03$Delivery.Mode,levels = c("Vaginal","C-Section"))
total_vertebrate_metadata_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(total_vertebrate_metadata_0.75_0.03$Age.at.Gluten.Introduction..months.)
total_vertebrate_metadata_0.75_0.03$Dx.Status <- factor(total_vertebrate_metadata_0.75_0.03$Dx.Status,levels = c("CONTROL","CELIAC"))
total_vertebrate_metadata_0.75_0.03$onset_timeline_combined <- factor(total_vertebrate_metadata_0.75_0.03$onset_timeline_combined,levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))




# check if the column names of abundance table match the order of the row names of metadata table
all(rownames(metadata.clean) == colnames(vertebrate_contig_abundance.clean))
# TRUE

count_matrix <- as.matrix(vertebrate_contig_abundance.clean[, sapply(vertebrate_contig_abundance.clean, is.numeric)])
sparsity <- sum(count_matrix == 0) / length(count_matrix)
# 0.8886535





## -----------------------------------------------------------------------------
library(data.table)

# # filter contig count table to only keep vertebrate contigs
# US_vertebrate_contig_count_table <- contig_count_table %>%
#                                  filter(V2 %in% vertebrate_contig_IDs) %>%
#                                  filter(grepl("_01_GEMM_", V1))
# 
# 
# write.table(US_vertebrate_contig_count_table,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_Coverm_concatenated_vertebrate_contigs.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 
# 
# # filter it by 75% coverage rate
# filter_and_save_wide_data(US_vertebrate_contig_count_table, 
#                           covered_fraction_threshold = 0.75, 
#                           output_path = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US") 
# 


# load the 75% mapping rate filtered contig table
US_vertebrate_contig_processed <- readRDS("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/wide_rpkm_genes_75Cov.rds") %>%
                 column_to_rownames("ORFID")


# then let's check colSums and ensure there is no sample have all 0s across orfs
# sample_all_0 <- names(colSums(US_vertebrate_contig_processed)[colSums(US_vertebrate_contig_processed) == 0])
# [1] "NovaSeq_N966_I12886_39602_Celiac_Leonard_Stool_01_GEMM_119_30M"       "NovaSeq_N966_I12920_39636_Celiac_Leonard_Stool_01_GEMM_154_30M"      
# [3] "NovaSeq_N978_I13149_39572_Celiac_Leonard_Stool_01_GEMM_093_24M"       "NovaSeq_N978_I13150_39573_Celiac_Leonard_Stool_01_GEMM_093_30M"      
# [5] "NovaSeq_N978_I13151_39574_Celiac_Leonard_Stool_01_GEMM_093_36M"       "NovaSeq_N978_I13152_39575_Celiac_Leonard_Stool_01_GEMM_093_48M"      
# [7] "NovaSeq_N978_I13166_39589_Celiac_Leonard_Stool_01_GEMM_106_24M"       "NovaSeq_N980_I13250_39529_Leonard_Human_Stool_Celiac_01_GEMM_057_60M"
# [9] "NovaSeq_N983_I13413_39587_Celiac_Leonard_Stool_01_GEMM_106_12M" 


# let's remove the 6 0 counts samples
US_vertebrate_contig_processed <- US_vertebrate_contig_processed %>%
                               select(-all_of(sample_all_0))
dim(US_vertebrate_contig_processed)
# 23366   199


# assign vertebrate_contig_processed to vertebrate_contig_count, so vertebrate_contig_processed won't be changed in downstream analysis
US_vertebrate_contig_count <- US_vertebrate_contig_processed

colnames(US_vertebrate_contig_count) <- gsub("NovaSeq_N966_|NovaSeq_N978_|NovaSeq_N979_|NovaSeq_N980_|NovaSeq_N981_|NovaSeq_N982_|NovaSeq_N983_|_stats","",colnames(US_vertebrate_contig_count))

colnames(US_vertebrate_contig_count) <- str_replace_all(colnames(US_vertebrate_contig_count),"Leonard_Human_stool_Celiac","Celiac_Leonard_Stool")
colnames(US_vertebrate_contig_count) <- sapply(strsplit(colnames(US_vertebrate_contig_count),"_Stool_|_Stool_Celiac_"),"[",2)
colnames(US_vertebrate_contig_count) <- ifelse(colnames(US_vertebrate_contig_count) %like% "_6Y",gsub("_6Y","_72M",colnames(US_vertebrate_contig_count)),
                                           ifelse(colnames(US_vertebrate_contig_count) %like% "_7Y",gsub("_7Y","_84M",colnames(US_vertebrate_contig_count)),colnames(US_vertebrate_contig_count)))

US_vertebrate_contig_count.table_0.75 <- t(US_vertebrate_contig_count) 
dim(US_vertebrate_contig_count.table_0.75)
# 199 23366



US_vertebrate_contig_count.table_0.75 <- merge(metadata,US_vertebrate_contig_count.table_0.75,by = 0) %>% # 189 samples left
                                      filter(feeding_first_year != "Unknown") %>%  # 189 samples left
                                      filter(Age.at.Gluten.Introduction..months. != "Unknown") # 185 sample left
dim(US_vertebrate_contig_count.table_0.75)
# 185 samples x  23393 metadata+contigs
# 31237 contigs


all_0_contigs <- names(colSums(US_vertebrate_contig_count.table_0.75[28:23393])[colSums(US_vertebrate_contig_count.table_0.75[28:23393]) == 0])
# 654 contigs become 0 through all samples after the 14 samples were removed from the above step
# so we need to remove those contigs

US_vertebrate_contig_count.table_0.75 <- US_vertebrate_contig_count.table_0.75 %>%
                            select(-all_of(all_0_contigs))
dim(US_vertebrate_contig_count.table_0.75)
# 185 samples x 22739 metadata+contigs
# 22712 contigs



US_vertebrate_contig_abundance.table <- US_vertebrate_contig_count.table_0.75 %>%
                       select(c("Row.names",28:22739)) %>%
                       column_to_rownames("Row.names") %>% 
                       t(.) %>%
                       data.frame(.)
colnames(US_vertebrate_contig_abundance.table) <- gsub("X","",colnames(US_vertebrate_contig_abundance.table))


# check prevelence filtering
prevalence <- rowSums(US_vertebrate_contig_abundance.table > 0) / ncol(US_vertebrate_contig_abundance.table)
table(prevalence > 0.03)  # 3% threshold
# FALSE  TRUE 
# 19334  3378

# Filter Contigs with prevalence > 3%
US_vertebrate_contig_abundance.clean <- US_vertebrate_contig_abundance.table[prevalence > 0.03, ]
      
dim(US_vertebrate_contig_abundance.clean)
# 3378  185 (3% prev)

# check if there are any samples become 0 counts
names(colSums(US_vertebrate_contig_abundance.clean)[colSums(US_vertebrate_contig_abundance.clean) == 0])
# "01_GEMM_006_60M" "01_GEMM_016_12M" "01_GEMM_060_18M" "01_GEMM_102_12M" "01_GEMM_109_48M"

# remove the two 0 counts samples
# we will have to remove sample 02_GEMM_021_36M because it becomes 0 counts at this step
US_vertebrate_contig_abundance.clean <- US_vertebrate_contig_abundance.clean %>%
                       select(-c("01_GEMM_006_60M","01_GEMM_016_12M","01_GEMM_060_18M","01_GEMM_102_12M","01_GEMM_109_48M"))
dim(US_vertebrate_contig_abundance.clean)
# 3378  180


US_vertebrate_contig_abundance_0.75_0.03 <- US_vertebrate_contig_abundance.clean

# write.csv(US_vertebrate_contig_abundance.clean,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_abundance_0.75_0.03.csv",quote = FALSE,row.names = TRUE)
colnames(US_vertebrate_contig_abundance.clean[1:27])

US_metadata.clean <- US_vertebrate_contig_count.table_0.75[1:27] %>%
                  filter(!Row.names %in% c("01_GEMM_006_60M","01_GEMM_016_12M","01_GEMM_060_18M","01_GEMM_102_12M","01_GEMM_109_48M")) %>%
                  column_to_rownames("Row.names")
dim(US_metadata.clean)
# 180  26

US_vertebrate_metadata_0.75_0.03 <- US_metadata.clean


#write.csv(US_metadata.clean,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_metadata_0.75_0.03.txt",quote = FALSE,row.names = TRUE)



# set up categorical variable levels
US_vertebrate_metadata_0.75_0.03$feeding_first_year <- factor(US_vertebrate_metadata_0.75_0.03$feeding_first_year,levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
US_vertebrate_metadata_0.75_0.03$HLA.Category <- factor(US_vertebrate_metadata_0.75_0.03$HLA.Category,levels = c("Standard Risk","High Risk","Low/No Risk"))
US_vertebrate_metadata_0.75_0.03$Sex <- factor(US_vertebrate_metadata_0.75_0.03$Sex,levels = c("Female","Male"))
US_vertebrate_metadata_0.75_0.03$Delivery.Mode <- factor(US_vertebrate_metadata_0.75_0.03$Delivery.Mode,levels = c("Vaginal","C-Section"))
US_vertebrate_metadata_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(US_vertebrate_metadata_0.75_0.03$Age.at.Gluten.Introduction..months.)
US_vertebrate_metadata_0.75_0.03$Dx.Status <- factor(US_vertebrate_metadata_0.75_0.03$Dx.Status,levels = c("CONTROL","CELIAC"))
US_vertebrate_metadata_0.75_0.03$onset_timeline_combined <- factor(US_vertebrate_metadata_0.75_0.03$onset_timeline_combined,levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))


# check if the column names of abundance table match the order of the row names of metadata table
all(rownames(US_metadata.clean) == colnames(US_vertebrate_contig_abundance.clean))
# TRUE

count_matrix <- as.matrix(US_vertebrate_contig_abundance.clean[, sapply(US_vertebrate_contig_abundance.clean, is.numeric)])
sparsity <- sum(count_matrix == 0) / length(count_matrix)
# 0.9003503




## -----------------------------------------------------------------------------
# library(data.table)
# 
# # filter contig count table to only keep vertebrate contigs
# Italy_vertebrate_contig_count_table <- contig_count_table %>%
#                                  filter(V2 %in% vertebrate_contig_IDs) %>%
#                                  filter(grepl("_02_GEMM_", V1))
# 
# 
# write.table(Italy_vertebrate_contig_count_table,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_Coverm_concatenated_vertebrate_contigs.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 
# 
# # filter it by 75% coverage rate
# filter_and_save_wide_data(Italy_vertebrate_contig_count_table, 
#                           covered_fraction_threshold = 0.75, 
#                           output_path = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy") 
# 
# 

# load the 75% mapping rate filtered contig table
Italy_vertebrate_contig_processed <- readRDS("~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/wide_rpkm_genes_75Cov.rds") %>%
                 column_to_rownames("ORFID")


# then let's check colSums and ensure there is no sample have all 0s across orfs
sample_all_0 <- names(colSums(Italy_vertebrate_contig_processed)[colSums(Italy_vertebrate_contig_processed) == 0])
# [1] "NovaSeq_N981_I13300_39814_Celiac_Leonard_Stool_02_GEMM_021_36M"       "NovaSeq_N981_I13317_39831_Celiac_Leonard_Stool_02_GEMM_057_24M"      
# [3] "NovaSeq_N982_I13323_39837_Celiac_Leonard_Stool_02_GEMM_077_12M"       "NovaSeq_N983_I13417_39799_Leonard_Human_stool_Celiac_02_GEMM_006_60M"


# let's remove the 6 0 counts samples
Italy_vertebrate_contig_processed <- Italy_vertebrate_contig_processed %>%
                               select(-all_of(sample_all_0))
dim(Italy_vertebrate_contig_processed)
# 17174   128


# assign vertebrate_contig_processed to vertebrate_contig_count, so vertebrate_contig_processed won't be changed in downstream analysis
Italy_vertebrate_contig_count <- Italy_vertebrate_contig_processed

colnames(Italy_vertebrate_contig_count) <- gsub("NovaSeq_N966_|NovaSeq_N978_|NovaSeq_N979_|NovaSeq_N980_|NovaSeq_N981_|NovaSeq_N982_|NovaSeq_N983_|_stats","",colnames(Italy_vertebrate_contig_count))

colnames(Italy_vertebrate_contig_count) <- str_replace_all(colnames(Italy_vertebrate_contig_count),"Leonard_Human_stool_Celiac","Celiac_Leonard_Stool")
colnames(Italy_vertebrate_contig_count) <- sapply(strsplit(colnames(Italy_vertebrate_contig_count),"_Stool_|_Stool_Celiac_"),"[",2)
colnames(Italy_vertebrate_contig_count) <- ifelse(colnames(Italy_vertebrate_contig_count) %like% "_6Y",gsub("_6Y","_72M",colnames(Italy_vertebrate_contig_count)),
                                           ifelse(colnames(Italy_vertebrate_contig_count) %like% "_7Y",gsub("_7Y","_84M",colnames(Italy_vertebrate_contig_count)),colnames(Italy_vertebrate_contig_count)))

Italy_vertebrate_contig_count.table_0.75 <- t(Italy_vertebrate_contig_count) 
dim(Italy_vertebrate_contig_count.table_0.75)
# 128 17174


Italy_vertebrate_contig_count.table_0.75 <- merge(metadata,Italy_vertebrate_contig_count.table_0.75,by = 0) %>% # 126 samples left
                                      filter(feeding_first_year != "Unknown") %>%  # 110 samples left
                                      filter(Age.at.Gluten.Introduction..months. != "Unknown") # 110 sample left
dim(Italy_vertebrate_contig_count.table_0.75)
# 110 samples x  17201 metadata+contigs
# 17174 contigs


all_0_contigs <- names(colSums(Italy_vertebrate_contig_count.table_0.75[28:17201])[colSums(Italy_vertebrate_contig_count.table_0.75[28:17201]) == 0])
# 779 contigs become 0 through all samples after the 18 samples were removed from the above step
# so we need to remove those contigs

Italy_vertebrate_contig_count.table_0.75 <- Italy_vertebrate_contig_count.table_0.75 %>%
                            select(-all_of(all_0_contigs))
dim(Italy_vertebrate_contig_count.table_0.75)
# 110 samples x 16422 metadata+contigs
# 16395 contigs


Italy_vertebrate_contig_abundance.table <- Italy_vertebrate_contig_count.table_0.75 %>%
                       select(c("Row.names",28:16422)) %>%
                       column_to_rownames("Row.names") %>% 
                       t(.) %>%
                       data.frame(.)
colnames(Italy_vertebrate_contig_abundance.table) <- gsub("X","",colnames(Italy_vertebrate_contig_abundance.table))


# check prevelence filtering
prevalence <- rowSums(Italy_vertebrate_contig_abundance.table > 0) / ncol(Italy_vertebrate_contig_abundance.table)
table(prevalence > 0.03)  # 3% threshold
# FALSE  TRUE 
# 12973  3422

# Filter Contigs with prevalence > 3%
Italy_vertebrate_contig_abundance.clean <- Italy_vertebrate_contig_abundance.table[prevalence > 0.03, ]
      
dim(Italy_vertebrate_contig_abundance.clean)
# 3422  110 (3% prev)

# check if there are any samples become 0 counts
names(colSums(Italy_vertebrate_contig_abundance.clean)[colSums(Italy_vertebrate_contig_abundance.clean) == 0])
# "02_GEMM_006_30M"

# remove the two 0 counts samples
# we will have to remove sample 02_GEMM_021_36M because it becomes 0 counts at this step
Italy_vertebrate_contig_abundance.clean <- Italy_vertebrate_contig_abundance.clean %>%
                       select(-c("02_GEMM_006_30M"))
dim(Italy_vertebrate_contig_abundance.clean)
# 3422  109

Italy_vertebrate_contig_abundance_0.75_0.03 <- Italy_vertebrate_contig_abundance.clean


#write.csv(Italy_vertebrate_contig_abundance.clean,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_abundance_0.75_0.03.csv",quote = FALSE,row.names = TRUE)


Italy_metadata.clean <- Italy_vertebrate_contig_count.table_0.75[1:27] %>%
                  filter(!Row.names %in% c("02_GEMM_006_30M")) %>%
                  column_to_rownames("Row.names")
dim(Italy_metadata.clean)
# 109  26

Italy_vertebrate_metadata_0.75_0.03 <- Italy_metadata.clean

#write.csv(Italy_metadata.clean,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_metadata_0.75_0.03.csv",quote = FALSE,row.names = TRUE)


# set up categorical variable levels
Italy_vertebrate_metadata_0.75_0.03$feeding_first_year <- factor(Italy_vertebrate_metadata_0.75_0.03$feeding_first_year,levels = c("Breast_fed","Formula","Breastmilk_and_formula"))
Italy_vertebrate_metadata_0.75_0.03$HLA.Category <- factor(Italy_vertebrate_metadata_0.75_0.03$HLA.Category,levels = c("Standard Risk","High Risk","Low/No Risk"))
Italy_vertebrate_metadata_0.75_0.03$Sex <- factor(Italy_vertebrate_metadata_0.75_0.03$Sex,levels = c("Female","Male"))
Italy_vertebrate_metadata_0.75_0.03$Delivery.Mode <- factor(Italy_vertebrate_metadata_0.75_0.03$Delivery.Mode,levels = c("Vaginal","C-Section"))
Italy_vertebrate_metadata_0.75_0.03$Age.at.Gluten.Introduction..months. <- as.numeric(Italy_vertebrate_metadata_0.75_0.03$Age.at.Gluten.Introduction..months.)
Italy_vertebrate_metadata_0.75_0.03$Dx.Status <- factor(Italy_vertebrate_metadata_0.75_0.03$Dx.Status,levels = c("CONTROL","CELIAC"))
Italy_vertebrate_metadata_0.75_0.03$onset_timeline_combined <- factor(Italy_vertebrate_metadata_0.75_0.03$onset_timeline_combined,levels = c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42"))


# check if the column names of abundance table match the order of the row names of metadata table
all(rownames(Italy_metadata.clean) == colnames(Italy_vertebrate_contig_abundance.clean))
# TRUE

count_matrix <- as.matrix(Italy_vertebrate_contig_abundance.clean[, sapply(Italy_vertebrate_contig_abundance.clean, is.numeric)])
sparsity <- sum(count_matrix == 0) / length(count_matrix)
# 0.8580126




## ----generate total PA table--------------------------------------------------

total.vertebrate_contig.target_cpm      <- 0.5   # e.g., 0.25–1.0 for very sparse data
total.vertebrate_contig.min_reads_floor <- 3L    # guard against singletons

## ---- 1) Align samples ------------------------------------------------------
stopifnot(all(colnames(total_vertebrate_contig_abundance_0.75_0.03) %in% rownames(total_vertebrate_metadata_0.75_0.03)))
# Reorder metadata to match counts columns
total_vertebrate_metadata_0.75_0.03 <- total_vertebrate_metadata_0.75_0.03[colnames(total_vertebrate_contig_abundance_0.75_0.03), , drop = FALSE]


## ---- 2) Build DGEList & normalize (TMMwsp) --------------------------------

total.vertebrate_contig.y <- DGEList(counts = total_vertebrate_contig_abundance_0.75_0.03, samples = total_vertebrate_metadata_0.75_0.03)
# Optional: drop genes with all zeros (purely empty rows)
total.vertebrate_contig.y <- total.vertebrate_contig.y[rowSums(total.vertebrate_contig.y$counts) > 0, , keep.lib.sizes = FALSE]

# Library-size normalization (robust for sparse data)
total.vertebrate_contig.y <- calcNormFactors(total.vertebrate_contig.y, method = "TMMwsp") 


## ---- 3) CPM-anchored presence/absence -------------------------------------
# Normalized CPMs (uses effective library sizes = lib.size * norm.factors)
total.vertebrate_contig.cpm_mat <- edgeR::cpm(total.vertebrate_contig.y, normalized.lib.sizes = TRUE)


# Presence rule: CPM >= target_cpm AND raw count >= min_reads_floor
total.vertebrate_contig.present_logical <- (total.vertebrate_contig.cpm_mat >= total.vertebrate_contig.target_cpm) & (total.vertebrate_contig.y$counts >= total.vertebrate_contig.min_reads_floor)


# Convert to 0/1 matrix with dimnames preserved
total.vertebrate_contig.PA <- matrix(
  as.integer(total.vertebrate_contig.present_logical),
  nrow = nrow(total.vertebrate_contig.present_logical),
  ncol = ncol(total.vertebrate_contig.present_logical),
  dimnames = dimnames(total.vertebrate_contig.present_logical)
)

dim(total.vertebrate_contig.PA)
# 3511  289





## ----generate US PA table-----------------------------------------------------

US.vertebrate_contig.target_cpm  <- 0.5   # e.g., 0.25–1.0 for very sparse data
US.vertebrate_contig.min_reads_floor <- 3L    # guard against singletons

## ---- 1) Align samples ------------------------------------------------------
stopifnot(all(colnames(US_vertebrate_contig_abundance_0.75_0.03) %in% rownames(US_vertebrate_metadata_0.75_0.03)))
# Reorder metadata to match counts columns
US_vertebrate_metadata_0.75_0.03 <- US_vertebrate_metadata_0.75_0.03[colnames(US_vertebrate_contig_abundance_0.75_0.03), , drop = FALSE]
sum(US_vertebrate_contig_abundance_0.75_0.03 > 0)


## ---- 2) Build DGEList & normalize (TMMwsp) --------------------------------

US.vertebrate_contig.y <- DGEList(counts = US_vertebrate_contig_abundance_0.75_0.03, samples = US_vertebrate_metadata_0.75_0.03)
# Optional: drop genes with all zeros (purely empty rows)
US.vertebrate_contig.y <- US.vertebrate_contig.y[rowSums(US.vertebrate_contig.y$counts) > 0, , keep.lib.sizes = FALSE]


# Library-size normalization (robust for sparse data)
US.vertebrate_contig.y <- calcNormFactors(US.vertebrate_contig.y, method = "TMMwsp")


## ---- 3) CPM-anchored presence/absence -------------------------------------
# Normalized CPMs (uses effective library sizes = lib.size * norm.factors)
US.vertebrate_contig.cpm_mat <- edgeR::cpm(US.vertebrate_contig.y, normalized.lib.sizes = TRUE)


# Presence rule: CPM >= target_cpm AND raw count >= min_reads_floor
US.vertebrate_contig.present_logical <- (US.vertebrate_contig.cpm_mat >= US.vertebrate_contig.target_cpm) & (US.vertebrate_contig.y$counts >= US.vertebrate_contig.min_reads_floor)
table(US.vertebrate_contig.present_logical)
#  FALSE   TRUE 
# 547452  60588 

# every non-zero count passed

# Convert to 0/1 matrix with dimnames preserved
US.vertebrate_contig.PA <- matrix(
  as.integer(US.vertebrate_contig.present_logical),
  nrow = nrow(US.vertebrate_contig.present_logical),
  ncol = ncol(US.vertebrate_contig.present_logical),
  dimnames = dimnames(US.vertebrate_contig.present_logical)
)

dim(US.vertebrate_contig.PA)
# 3378  180




## ----generate Italy PA table--------------------------------------------------


Italy.vertebrate_contig.target_cpm      <- 0.5   # e.g., 0.25–1.0 for very sparse data
Italy.vertebrate_contig.min_reads_floor <- 3L    # guard against singletons

## ---- 1) Align samples ------------------------------------------------------
stopifnot(all(colnames(Italy_vertebrate_contig_abundance_0.75_0.03) %in% rownames(Italy_vertebrate_metadata_0.75_0.03)))
# Reorder metadata to match counts columns
Italy_vertebrate_metadata_0.75_0.03 <- Italy_vertebrate_metadata_0.75_0.03[colnames(Italy_vertebrate_contig_abundance_0.75_0.03), , drop = FALSE]

## ---- 2) Build DGEList & normalize (TMMwsp) --------------------------------

Italy.vertebrate_contig.y <- DGEList(counts = Italy_vertebrate_contig_abundance_0.75_0.03, samples = Italy_vertebrate_metadata_0.75_0.03)
# Optional: drop genes with all zeros (purely empty rows)
Italy.vertebrate_contig.y <- Italy.vertebrate_contig.y[rowSums(Italy.vertebrate_contig.y$counts) > 0, , keep.lib.sizes = FALSE]

# Library-size normalization (robItalyt for sparse data)
Italy.vertebrate_contig.y <- calcNormFactors(Italy.vertebrate_contig.y, method = "TMMwsp")

## ---- 3) CPM-anchored presence/absence -------------------------------------
# Normalized CPMs (Italyes effective library sizes = lib.size * norm.factors)
Italy.vertebrate_contig.cpm_mat <- edgeR::cpm(Italy.vertebrate_contig.y, normalized.lib.sizes = TRUE)

# Presence rule: CPM >= target_cpm AND raw count >= min_reads_floor
Italy.vertebrate_contig.present_logical <- (Italy.vertebrate_contig.cpm_mat >= Italy.vertebrate_contig.target_cpm) & (Italy.vertebrate_contig.y$counts >= Italy.vertebrate_contig.min_reads_floor)

# Convert to 0/1 matrix with dimnames preserved
Italy.vertebrate_contig.PA <- matrix(
  as.integer(Italy.vertebrate_contig.present_logical),
  nrow = nrow(Italy.vertebrate_contig.present_logical),
  ncol = ncol(Italy.vertebrate_contig.present_logical),
  dimnames = dimnames(Italy.vertebrate_contig.present_logical)
)

dim(Italy.vertebrate_contig.PA)
# 3422  109




## ----PA_glmmTMB_logit_total---------------------------------------------------


# Create long-format data for total cohort
vertebrate_contig_total_pa_long <- total.vertebrate_contig.PA %>%
  as.data.frame() %>%
  rownames_to_column("vertebrate_contig_id") %>%
  pivot_longer(cols = -vertebrate_contig_id, names_to = "sample_id", values_to = "PA") %>%
  left_join(total_vertebrate_metadata_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id")

table(vertebrate_contig_total_pa_long$PA)
#      0      1 
#  901699 112980


# Fit model for each ORF
vertebrate_contig_total_glmmTMB_logit_results <- vertebrate_contig_total_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + Country + HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                family = binomial(link = "logit"), data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(vertebrate_contig_id, tidy_results) %>%
  unnest(tidy_results)



vertebrate_contig_total_PA.adjusted.results <- vertebrate_contig_total_glmmTMB_logit_results %>%
  filter(!is.na(p.value)) %>%
  group_by(term) %>%   # often you adjust per term of interest
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

dim(vertebrate_contig_total_PA.adjusted.results %>% filter(p.adj < 0.05))
# 228




## -----------------------------------------------------------------------------

# ---- add right after your vertebrate_contig_total_PA.adjusted.results code ----
library(emmeans)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(glmmTMB)

# 0) Define timepoint order (must include "t0" first)
tp_levels <- c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")
stopifnot("t0" %in% tp_levels)


# 1) Refit per-orf models (just for emmeans; your original table is untouched)
vertebrate_contig_total_PA_models_for_emm <- vertebrate_contig_total_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(
        PA ~ Dx.Status * onset_timeline_combined + Country + Sex +
             Age.at.Gluten.Introduction..months. +
             HLA.Category + feeding_first_year + Delivery.Mode +
             (1 | patientID),
        family = binomial(link = "logit"),
        data   = .x
      ),
      error = function(e) NULL
    ))
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  select(vertebrate_contig_id, model)



# 2) Helper: get CELIAC TP vs t0 contrasts from emmeans (with SE/CI/p)
# --- Safe helper: CELIAC TP vs t0 with emmeans (handles z.ratio/t.ratio/none) ---
# SAFE helper: CELIAC TP vs t0 with emmeans (uses base::summary, robust CIs/statistic)
extract_celiac_tp_vs_t0 <- function(mod, tp_levels, dx_level = "CELIAC") {
  
  print("test0")
  if (is.null(mod) || !("t0" %in% tp_levels)) return(tibble())

  # emmeans across timepoints within Dx
  emm <- tryCatch(

    emmeans(
      mod, ~ onset_timeline_combined | Dx.Status,
      at = list(onset_timeline_combined = tp_levels),
      type = "link"  # log-odds scale
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())


  # contrasts: each TP vs t0 (within each Dx)
  ref_idx <- which(tp_levels == "t0")[1]
  con <- tryCatch(
    contrast(emm, method = "trt.vs.ctrl", ref = ref_idx),
    error = function(e) NULL
  )

  if (is.null(con)) return(tibble())


  # IMPORTANT: use base::summary(), not emmeans::summary()
  s <- tibble::as_tibble(summary(con, infer = c(TRUE, TRUE)))


  # unify 'statistic' column (z.ratio OR t.ratio OR NA)
  if (!("statistic" %in% names(s))) {

    zr <- if ("z.ratio" %in% names(s)) s$z.ratio else NULL
    tr <- if ("t.ratio" %in% names(s)) s$t.ratio else NULL
    s$statistic <- dplyr::coalesce(zr, tr, rep(NA_real_, nrow(s)))
  }
  

  # ensure CI columns
  if (!all(c("lower.CL","upper.CL") %in% names(s))) {

    ci <- tryCatch(tibble::as_tibble(confint(con)), error = function(e) NULL)
    if (!is.null(ci)) {

      byk <- intersect(names(ci), names(s))
      s <- dplyr::left_join(
        s,
        dplyr::select(ci, dplyr::any_of(c(byk, "lower.CL","upper.CL","asymp.LCL","asymp.UCL"))),
        by = byk
      )
      if (!"lower.CL" %in% names(s) && "asymp.LCL" %in% names(s)) s$lower.CL <- s$asymp.LCL
      if (!"upper.CL" %in% names(s) && "asymp.UCL" %in% names(s)) s$upper.CL <- s$asymp.UCL
    }
  }
  if (!"lower.CL" %in% names(s)) s$lower.CL <- NA_real_
  if (!"upper.CL" %in% names(s)) s$upper.CL <- NA_real_


  s %>%
    dplyr::filter(Dx.Status == dx_level) %>%
    dplyr::mutate(
      TP        = stringr::str_remove(contrast, " - t0$"),
      term      = paste0("TP (", dx_level, "): ", TP, " vs t0"),
      OR        = exp(estimate)
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error = SE,
      statistic,
      p.value,
      conf.low  = lower.CL,
      conf.high = upper.CL,
      OR
    ) %>%
    dplyr::filter(!grepl("^TP \\(.+\\):\\s*t0\\s+vs\\s+t0$", term))
}


# --- Build the CELIAC rows and append to your original table -------------------
vertebrate_contig_total_PA_celiac_tp_rows <- vertebrate_contig_total_PA_models_for_emm %>%
  dplyr::mutate(out = purrr::map(model, ~ extract_celiac_tp_vs_t0(.x, tp_levels))) %>%
  dplyr::select(vertebrate_contig_id, out) %>%
  tidyr::unnest(out) %>%
  # same per-term FDR rule you used for the original table
  dplyr::group_by(term) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

vertebrate_contig_total_PA.with_celiac <- dplyr::bind_rows(
  vertebrate_contig_total_PA.adjusted.results,  # your original model-term stats
  vertebrate_contig_total_PA_celiac_tp_rows                     # added: CELIAC TP vs t0 stats from emmeans
)


dim(vertebrate_contig_total_PA.with_celiac %>% filter(p.adj < 0.05))
# 228  13





## -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom.mixed)
library(emmeans)
library(tibble)

# ---------------------------------------------------------------
# 1) Fit models (you already have this part)
# ---------------------------------------------------------------
total_vertebrate_contig_PA_models <- vertebrate_contig_total_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(
        PA ~ Dx.Status * onset_timeline_combined + Country + Sex + Age.at.Gluten.Introduction..months. +
          HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
        family = binomial(link = "logit"), data = .x
      ),
      error = function(e) NULL
    ))
  )

# ---------------------------------------------------------------
# 2) Extract between-group contrasts (CELIAC - CONTROL → flip)
# ---------------------------------------------------------------
extract_contrasts <- function(fit) {
  if (is.null(fit)) return(tibble())

  emm <- tryCatch(
    emmeans(fit, ~ Dx.Status | onset_timeline_combined),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  contr <- tryCatch(
    contrast(emm, method = "revpairwise", by = "onset_timeline_combined", adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(contr)) return(tibble())

  s  <- as_tibble(summary(contr))
  ci <- as_tibble(confint(contr)) %>%
        select(any_of(c("contrast","onset_timeline_combined","lower.CL","upper.CL")))

  keys <- intersect(names(s), names(ci))
  if (length(keys) > 0) {
    s <- left_join(s, ci, by = keys)
  }
  s
}

total_vertebrate_contig_PA_contrast_results <- total_vertebrate_contig_PA_models %>%
  mutate(contrasts = map(model, extract_contrasts)) %>%
  select(vertebrate_contig_id, contrasts) %>%
  tidyr::unnest(contrasts)

# adjust within each timepoint
total_vertebrate_contig_PA_contrast_results <- total_vertebrate_contig_PA_contrast_results %>%
  filter(!is.na(p.value)) %>%
  group_by(onset_timeline_combined) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ---------------------------------------------------------------
# 3) Flip CELIAC - CONTROL to CONTROL vs CELIAC
# ---------------------------------------------------------------
between_df <- total_vertebrate_contig_PA_contrast_results %>%
  filter(str_detect(contrast, "^CELIAC\\s*-\\s*CONTROL$"))

has_z  <- "z.ratio"  %in% names(between_df)
has_t  <- "t.ratio"  %in% names(between_df)
has_SE <- "SE"       %in% names(between_df)
has_l  <- "lower.CL" %in% names(between_df)
has_u  <- "upper.CL" %in% names(between_df)

between_tp_std <- between_df %>%
  mutate(
    est_flip   = -estimate,
    se_use     = if (has_SE) SE else NA_real_,
    stat_use   = if (has_z) z.ratio else if (has_t) t.ratio else NA_real_,
    low_flip   = if (has_l) -upper.CL else NA_real_,
    high_flip  = if (has_u) -lower.CL else NA_real_,
    OR         = exp(est_flip),
    term       = paste0("CONTROL vs CELIAC @ ", onset_timeline_combined),
    source     = "between_groups"
  ) %>%
  transmute(
    vertebrate_contig_id,
    term,
    estimate  = est_flip,
    std.error = se_use,
    statistic = stat_use,
    p.value,
    p.adj,
    conf.low  = low_flip,
    conf.high = high_flip,
    OR,
    source
  )

# ---------------------------------------------------------------
# 4) Harmonize with your main table (already created earlier)
#    Replace with the object that contains tidy terms + CELIAC contrasts
# ---------------------------------------------------------------
common_spec <- list(
  vertebrate_contig_id = NA_character_,
  term      = NA_character_,
  estimate  = NA_real_,
  std.error = NA_real_,
  statistic = NA_real_,
  p.value   = NA_real_,
  p.adj     = NA_real_,
  conf.low  = NA_real_,
  conf.high = NA_real_,
  OR        = NA_real_,
  source    = NA_character_
)

add_missing_cols_typed <- function(df, spec) {
  for (nm in names(spec)) if (!nm %in% names(df)) df[[nm]] <- spec[[nm]]
  df
}


# Use your actual main results object here
main_tbl <- vertebrate_contig_total_PA.with_celiac   # <-- adjust if named differently

main_std <- main_tbl %>%
  add_missing_cols_typed(common_spec) %>%
  mutate(
    source = as.character(source),
    source = ifelse(is.na(source) | source == "", "model_or_emmeans", source)
  ) %>%
  select(names(common_spec))

if (nrow(between_tp_std) == 0) {
  between_std <- tibble(
    vertebrate_contig_id = character(), term = character(),
    estimate = numeric(), std.error = numeric(), statistic = numeric(),
    p.value = numeric(), p.adj = numeric(),
    conf.low = numeric(), conf.high = numeric(),
    OR = numeric(), source = character()
  )
} else {
  between_std <- between_tp_std %>%
    add_missing_cols_typed(common_spec) %>%
    mutate(
      source = as.character(source),
      source = ifelse(is.na(source) | source == "", "between_groups", source)
    ) %>%
    select(names(common_spec))
}

# ---------------------------------------------------------------
# 5) Final combined table
# ---------------------------------------------------------------
vertebrate_contig_total_PA.everything <- bind_rows(main_std, between_std) %>%
  arrange(vertebrate_contig_id, term)

dim(vertebrate_contig_total_PA.everything %>% filter(p.adj < 0.05))
# 228  11


write.csv(vertebrate_contig_total_PA.everything,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_PA_results_everything.csv",row.names = FALSE)



## -----------------------------------------------------------------------------

library(dplyr)
library(stringr)
library(ggplot2)



tp_levels <- c("t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")  # add "t0-over42" if present


vertebrate_contig_total_PA_heatmaps <- auto_heatmaps_detected(
    df            = vertebrate_contig_total_PA.everything,
    tp_levels     = tp_levels,
    fdr           = 0.05,
    scale         = "log",        # Use log-odds scale
    min_contigs   = 1,
    analysis_type = "PA"          # Logistic regression labels
  )


ggsave(vertebrate_contig_total_PA_heatmaps$`TP (CONTROL)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_PA_results_everything_TP_CONTROL.pdf",dpi = 600, width = 6, height = 6)


ggsave(vertebrate_contig_total_PA_heatmaps$`TP (CELIAC)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_PA_results_everything_TP_CELIAC.pdf",dpi = 600, width = 6, height = 2.5)

ggsave(vertebrate_contig_total_PA_heatmaps$`Between groups`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_PA_results_everything_Between_groups.pdf",dpi = 600, width = 6, height = 6)

ggsave(vertebrate_contig_total_PA_heatmaps$`Dx × TP interactions`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_PA_results_everything_interactions.pdf",dpi = 600, width = 6, height = 6)



## ----PA_glmmTMB_logit_US------------------------------------------------------
# Create long-format data for US cohort
vertebrate_contig_US_pa_long <- US.vertebrate_contig.PA %>%
  as.data.frame() %>%
  rownames_to_column("vertebrate_contig_id") %>%
  pivot_longer(cols = -vertebrate_contig_id, names_to = "sample_id", values_to = "PA") %>%
  left_join(US_vertebrate_metadata_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id")

# Fit model for each contig
vertebrate_contig_US_PA_glmmTMB_logit_results <- vertebrate_contig_US_pa_long %>%
    group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                family = binomial(link = "logit"), data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(vertebrate_contig_id, tidy_results) %>%
  unnest(tidy_results)


vertebrate_contig_US_PA.adjusted.results <- vertebrate_contig_US_PA_glmmTMB_logit_results %>%
  filter(!is.na(p.value)) %>%
  group_by(term) %>%   # often you adjust per term of interest
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

dim(vertebrate_contig_US_PA.adjusted.results %>% filter(p.adj < 0.05))
# 3




## -----------------------------------------------------------------------------

# ---- add right after your orf_total_PA.adjusted.results code ----


# 0) Define timepoint order (must include "t0" first)
tp_levels <- c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")
stopifnot("t0" %in% tp_levels)

# 1) Refit per-orf models (just for emmeans; your original table is untouched)
vertebrate_contig_US_PA_models_for_emm <- vertebrate_contig_US_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(
        PA ~ Dx.Status * onset_timeline_combined + Sex +
             Age.at.Gluten.Introduction..months. +
             HLA.Category + feeding_first_year + Delivery.Mode +
             (1 | patientID),
        family = binomial(link = "logit"),
        data   = .x
      ),
      error = function(e) NULL
    ))
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  select(vertebrate_contig_id, model)

# 2) Helper: get CELIAC TP vs t0 contrasts from emmeans (with SE/CI/p)
# --- Safe helper: CELIAC TP vs t0 with emmeans (handles z.ratio/t.ratio/none) ---
# SAFE helper: CELIAC TP vs t0 with emmeans (uses base::summary, robust CIs/statistic)
extract_celiac_tp_vs_t0 <- function(mod, tp_levels, dx_level = "CELIAC") {
  if (is.null(mod) || !("t0" %in% tp_levels)) return(tibble())

  # emmeans across timepoints within Dx
  emm <- tryCatch(
    emmeans(
      mod, ~ onset_timeline_combined | Dx.Status,
      at = list(onset_timeline_combined = tp_levels),
      type = "link"  # log-odds scale
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  # contrasts: each TP vs t0 (within each Dx)
  ref_idx <- which(tp_levels == "t0")[1]
  con <- tryCatch(
    contrast(emm, method = "trt.vs.ctrl", ref = ref_idx),
    error = function(e) NULL
  )
  if (is.null(con)) return(tibble())

  # IMPORTANT: use base::summary(), not emmeans::summary()
  s <- tibble::as_tibble(summary(con, infer = c(TRUE, TRUE)))

  # unify 'statistic' column (z.ratio OR t.ratio OR NA)
  if (!("statistic" %in% names(s))) {
    zr <- if ("z.ratio" %in% names(s)) s$z.ratio else NULL
    tr <- if ("t.ratio" %in% names(s)) s$t.ratio else NULL
    s$statistic <- dplyr::coalesce(zr, tr, rep(NA_real_, nrow(s)))
  }

  # ensure CI columns
  if (!all(c("lower.CL","upper.CL") %in% names(s))) {
    ci <- tryCatch(tibble::as_tibble(confint(con)), error = function(e) NULL)
    if (!is.null(ci)) {
      byk <- intersect(names(ci), names(s))
      s <- dplyr::left_join(
        s,
        dplyr::select(ci, dplyr::any_of(c(byk, "lower.CL","upper.CL","asymp.LCL","asymp.UCL"))),
        by = byk
      )
      if (!"lower.CL" %in% names(s) && "asymp.LCL" %in% names(s)) s$lower.CL <- s$asymp.LCL
      if (!"upper.CL" %in% names(s) && "asymp.UCL" %in% names(s)) s$upper.CL <- s$asymp.UCL
    }
  }
  if (!"lower.CL" %in% names(s)) s$lower.CL <- NA_real_
  if (!"upper.CL" %in% names(s)) s$upper.CL <- NA_real_

  s %>%
    dplyr::filter(Dx.Status == dx_level) %>%
    dplyr::mutate(
      TP        = stringr::str_remove(contrast, " - t0$"),
      term      = paste0("TP (", dx_level, "): ", TP, " vs t0"),
      OR        = exp(estimate)
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error = SE,
      statistic,
      p.value,
      conf.low  = lower.CL,
      conf.high = upper.CL,
      OR
    ) %>%
    dplyr::filter(!grepl("^TP \\(.+\\):\\s*t0\\s+vs\\s+t0$", term))
}


# --- Build the CELIAC rows and append to your original table -------------------
vertebrate_contig_US_PA_celiac_tp_rows <- vertebrate_contig_US_PA_models_for_emm %>%
  dplyr::mutate(out = purrr::map(model, ~ extract_celiac_tp_vs_t0(.x, tp_levels))) %>%
  dplyr::select(vertebrate_contig_id, out) %>%
  tidyr::unnest(out) %>%
  # same per-term FDR rule you used for the original table
  dplyr::group_by(term) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

vertebrate_contig_US_PA.with_celiac <- dplyr::bind_rows(
  vertebrate_contig_US_PA.adjusted.results,  # your original model-term stats
  vertebrate_contig_US_PA_celiac_tp_rows                     # added: CELIAC TP vs t0 stats from emmeans
)


dim(vertebrate_contig_US_PA.with_celiac %>% filter(p.adj < 0.05))
# 3


# write.csv(orf_US_PA.with_celiac, "~/Handley Lab Dropbox/16S/Celiac/Phage/phage_detection_pipeline_new_assembly/Virome/Orf_orf_Phrog_compositional/orf_orf_results_figures/orf/US_orf_PA_model1_glmmTMB_logit_adjusted.csv", row.names = FALSE)



## -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom.mixed)
library(emmeans)
library(tibble)

# ---------------------------------------------------------------
# 1) Fit models (you already have this part)
# ---------------------------------------------------------------
US_vertebrate_contig_PA_models <- vertebrate_contig_US_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(
        PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. +
          HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
        family = binomial(link = "logit"), data = .x
      ),
      error = function(e) NULL
    ))
  )


# ---------------------------------------------------------------
# 2) Extract between-group contrasts (CELIAC - CONTROL → flip)
# ---------------------------------------------------------------
extract_contrasts <- function(fit) {
  if (is.null(fit)) return(tibble())

  emm <- tryCatch(
    emmeans(fit, ~ Dx.Status | onset_timeline_combined),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  contr <- tryCatch(
    contrast(emm, method = "revpairwise", by = "onset_timeline_combined", adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(contr)) return(tibble())

  s  <- as_tibble(summary(contr))
  ci <- as_tibble(confint(contr)) %>%
        select(any_of(c("contrast","onset_timeline_combined","lower.CL","upper.CL")))

  keys <- intersect(names(s), names(ci))
  if (length(keys) > 0) {
    s <- left_join(s, ci, by = keys)
  }
  s
}

US_vertebrate_contig_PA_contrast_results <- US_vertebrate_contig_PA_models %>%
  mutate(contrasts = map(model, extract_contrasts)) %>%
  select(vertebrate_contig_id, contrasts) %>%
  tidyr::unnest(contrasts)

# adjust within each timepoint
US_vertebrate_contig_PA_contrast_results <- US_vertebrate_contig_PA_contrast_results %>%
  filter(!is.na(p.value)) %>%
  group_by(onset_timeline_combined) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()


# ---------------------------------------------------------------
# 3) Flip CELIAC - CONTROL to CONTROL vs CELIAC
# ---------------------------------------------------------------
between_df <- US_vertebrate_contig_PA_contrast_results %>%
  filter(str_detect(contrast, "^CELIAC\\s*-\\s*CONTROL$"))

has_z  <- "z.ratio"  %in% names(between_df)
has_t  <- "t.ratio"  %in% names(between_df)
has_SE <- "SE"       %in% names(between_df)
has_l  <- "lower.CL" %in% names(between_df)
has_u  <- "upper.CL" %in% names(between_df)

between_tp_std <- between_df %>%
  mutate(
    est_flip   = -estimate,
    se_use     = if (has_SE) SE else NA_real_,
    stat_use   = if (has_z) z.ratio else if (has_t) t.ratio else NA_real_,
    low_flip   = if (has_l) -upper.CL else NA_real_,
    high_flip  = if (has_u) -lower.CL else NA_real_,
    OR         = exp(est_flip),
    term       = paste0("CONTROL vs CELIAC @ ", onset_timeline_combined),
    source     = "between_groups"
  ) %>%
  transmute(
    vertebrate_contig_id,
    term,
    estimate  = est_flip,
    std.error = se_use,
    statistic = stat_use,
    p.value,
    p.adj,
    conf.low  = low_flip,
    conf.high = high_flip,
    OR,
    source
  )

# ---------------------------------------------------------------
# 4) Harmonize with your main table (already created earlier)
#    Replace with the object that contains tidy terms + CELIAC contrasts
# ---------------------------------------------------------------
common_spec <- list(
  vertebrate_contig_id = NA_character_,
  term      = NA_character_,
  estimate  = NA_real_,
  std.error = NA_real_,
  statistic = NA_real_,
  p.value   = NA_real_,
  p.adj     = NA_real_,
  conf.low  = NA_real_,
  conf.high = NA_real_,
  OR        = NA_real_,
  source    = NA_character_
)

add_missing_cols_typed <- function(df, spec) {
  for (nm in names(spec)) if (!nm %in% names(df)) df[[nm]] <- spec[[nm]]
  df
}


# Use your actual main results object here
main_tbl <- vertebrate_contig_US_PA.with_celiac   # <-- adjust if named differently

main_std <- main_tbl %>%
  add_missing_cols_typed(common_spec) %>%
  mutate(
    source = as.character(source),
    source = ifelse(is.na(source) | source == "", "model_or_emmeans", source)
  ) %>%
  select(names(common_spec))

if (nrow(between_tp_std) == 0) {
  between_std <- tibble(
    vertebrate_contig_id = character(), term = character(),
    estimate = numeric(), std.error = numeric(), statistic = numeric(),
    p.value = numeric(), p.adj = numeric(),
    conf.low = numeric(), conf.high = numeric(),
    OR = numeric(), source = character()
  )
} else {
  between_std <- between_tp_std %>%
    add_missing_cols_typed(common_spec) %>%
    mutate(
      source = as.character(source),
      source = ifelse(is.na(source) | source == "", "between_groups", source)
    ) %>%
    select(names(common_spec))
}

# ---------------------------------------------------------------
# 5) Final combined table
# ---------------------------------------------------------------
vertebrate_contig_US_PA.everything <- bind_rows(main_std, between_std) %>%
  arrange(vertebrate_contig_id, term)


dim(vertebrate_contig_US_PA.everything %>% filter(p.adj < 0.05))
# 3


write.csv(vertebrate_contig_US_PA.everything,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_PA_results_everything.csv",row.names = FALSE)



## -----------------------------------------------------------------------------

tp_levels <- c("t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")  # add "t0-over42" if present


vertebrate_contig_US_PA_heatmaps <- auto_heatmaps_detected(
    df            = vertebrate_contig_US_PA.everything,
    tp_levels     = tp_levels,
    fdr           = 0.05,
    scale         = "log",        # Use log-odds scale
    min_contigs      = 1,
    analysis_type = "PA"          # Logistic regression labels
  )


ggsave(vertebrate_contig_US_PA_heatmaps$`TP (CONTROL)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_PA_results_everything_TP_CONTROL.pdf",dpi = 600, width = 6, height = 6)


ggsave(vertebrate_contig_US_PA_heatmaps$`TP (CELIAC)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_PA_results_everything_TP_CELIAC.pdf",dpi = 600, width = 6, height = 2.5)

ggsave(vertebrate_contig_US_PA_heatmaps$`Between groups`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_PA_results_everything_Between_groups.pdf",dpi = 600, width = 6, height = 6)

ggsave(vertebrate_contig_US_PA_heatmaps$`Dx × TP interactions`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_PA_results_everything_interactions.pdf",dpi = 600, width = 6, height = 6)



## ----PA_glmmTMB_logit_Italy---------------------------------------------------
# Create long-format data for Italy cohort
vertebrate_contig_Italy_pa_long <- Italy.vertebrate_contig.PA %>%
  as.data.frame() %>%
  rownames_to_column("vertebrate_contig_id") %>%
  pivot_longer(cols = -vertebrate_contig_id, names_to = "sample_id", values_to = "PA") %>%
  left_join(Italy_vertebrate_metadata_0.75_0.03 %>% rownames_to_column("sample_id"), by = "sample_id")

table(vertebrate_contig_Italy_pa_long$PA)

# Fit model for each contig
vertebrate_contig_Italy_glmmTMB_logit_results <- vertebrate_contig_Italy_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
                family = binomial(link = "logit"), data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(vertebrate_contig_id, tidy_results) %>%
  unnest(tidy_results)


vertebrate_contig_Italy_PA.adjusted.results <- vertebrate_contig_Italy_glmmTMB_logit_results %>%
  filter(!is.na(p.value)) %>%
  group_by(term) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

dim(vertebrate_contig_Italy_PA.adjusted.results %>% filter(p.adj < 0.05))




## -----------------------------------------------------------------------------

# 0) Define timepoint order (must include "t0" first)
tp_levels <- c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")
stopifnot("t0" %in% tp_levels)

# 1) Refit per-contig models (just for emmeans; your original table is untouched)
vertebrate_contig_Italy_PA_models_for_emm <- vertebrate_contig_Italy_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(
        PA ~ Dx.Status * onset_timeline_combined + Sex +
             Age.at.Gluten.Introduction..months. +
             HLA.Category + feeding_first_year + Delivery.Mode +
             (1 | patientID),
        family = binomial(link = "logit"),
        data   = .x
      ),
      error = function(e) NULL
    ))
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  select(vertebrate_contig_id, model)

# 2) Helper: get CELIAC TP vs t0 contrasts from emmeans (with SE/CI/p)
extract_celiac_tp_vs_t0 <- function(mod, tp_levels, dx_level = "CELIAC") {
  if (is.null(mod) || !("t0" %in% tp_levels)) return(tibble())

  # emmeans across timepoints within Dx
  emm <- tryCatch(
    emmeans(
      mod, ~ onset_timeline_combined | Dx.Status,
      at = list(onset_timeline_combined = tp_levels),
      type = "link"  # log-odds scale
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  # contrasts: each TP vs t0 (within each Dx)
  ref_idx <- which(tp_levels == "t0")[1]
  con <- tryCatch(
    contrast(emm, method = "trt.vs.ctrl", ref = ref_idx),
    error = function(e) NULL
  )
  if (is.null(con)) return(tibble())

  # IMPORTANT: use base::summary(), not emmeans::summary()
  s <- tibble::as_tibble(summary(con, infer = c(TRUE, TRUE)))

  # unify 'statistic' column (z.ratio OR t.ratio OR NA)
  if (!("statistic" %in% names(s))) {
    zr <- if ("z.ratio" %in% names(s)) s$z.ratio else NULL
    tr <- if ("t.ratio" %in% names(s)) s$t.ratio else NULL
    s$statistic <- dplyr::coalesce(zr, tr, rep(NA_real_, nrow(s)))
  }

  # ensure CI columns
  if (!all(c("lower.CL","upper.CL") %in% names(s))) {
    ci <- tryCatch(tibble::as_tibble(confint(con)), error = function(e) NULL)
    if (!is.null(ci)) {
      byk <- intersect(names(ci), names(s))
      s <- dplyr::left_join(
        s,
        dplyr::select(ci, dplyr::any_of(c(byk, "lower.CL","upper.CL","asymp.LCL","asymp.UCL"))),
        by = byk
      )
      if (!"lower.CL" %in% names(s) && "asymp.LCL" %in% names(s)) s$lower.CL <- s$asymp.LCL
      if (!"upper.CL" %in% names(s) && "asymp.UCL" %in% names(s)) s$upper.CL <- s$asymp.UCL
    }
  }
  if (!"lower.CL" %in% names(s)) s$lower.CL <- NA_real_
  if (!"upper.CL" %in% names(s)) s$upper.CL <- NA_real_

  s %>%
    dplyr::filter(Dx.Status == dx_level) %>%
    dplyr::mutate(
      TP        = stringr::str_remove(contrast, " - t0$"),
      term      = paste0("TP (", dx_level, "): ", TP, " vs t0"),
      OR        = exp(estimate)
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error = SE,
      statistic,
      p.value,
      conf.low  = lower.CL,
      conf.high = upper.CL,
      OR
    ) %>%
    dplyr::filter(!grepl("^TP \\(.+\\):\\s*t0\\s+vs\\s+t0$", term))
}


# --- Build the CELIAC rows and append to your original table -------------------
vertebrate_contig_Italy_PA_celiac_tp_rows <- vertebrate_contig_Italy_PA_models_for_emm %>%
  dplyr::mutate(out = purrr::map(model, ~ extract_celiac_tp_vs_t0(.x, tp_levels))) %>%
  dplyr::select(vertebrate_contig_id, out) %>%
  tidyr::unnest(out) %>%
  dplyr::group_by(term) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

vertebrate_contig_Italy_PA.with_celiac <- dplyr::bind_rows(
  vertebrate_contig_Italy_PA.adjusted.results,
  vertebrate_contig_Italy_PA_celiac_tp_rows
)

dim(vertebrate_contig_Italy_PA.with_celiac %>% filter(p.adj < 0.05))




## -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom.mixed)
library(emmeans)
library(tibble)

# ---------------------------------------------------------------
# 1) Fit models
# ---------------------------------------------------------------
Italy_vertebrate_contig_PA_models <- vertebrate_contig_Italy_pa_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(
        PA ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. +
          HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID),
        family = binomial(link = "logit"), data = .x
      ),
      error = function(e) NULL
    ))
  )

# ---------------------------------------------------------------
# 2) Extract between-group contrasts (CELIAC - CONTROL → flip)
# ---------------------------------------------------------------
extract_contrasts <- function(fit) {
  if (is.null(fit)) return(tibble())

  emm <- tryCatch(
    emmeans(fit, ~ Dx.Status | onset_timeline_combined),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  contr <- tryCatch(
    contrast(emm, method = "revpairwise", by = "onset_timeline_combined", adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(contr)) return(tibble())

  s  <- as_tibble(summary(contr))
  ci <- as_tibble(confint(contr)) %>%
        select(any_of(c("contrast","onset_timeline_combined","lower.CL","upper.CL")))

  keys <- intersect(names(s), names(ci))
  if (length(keys) > 0) {
    s <- left_join(s, ci, by = keys)
  }
  s
}

Italy_vertebrate_contig_PA_contrast_results <- Italy_vertebrate_contig_PA_models %>%
  mutate(contrasts = map(model, extract_contrasts)) %>%
  select(vertebrate_contig_id, contrasts) %>%
  tidyr::unnest(contrasts)

# adjust within each timepoint
Italy_vertebrate_contig_PA_contrast_results <- Italy_vertebrate_contig_PA_contrast_results %>%
  filter(!is.na(p.value)) %>%
  group_by(onset_timeline_combined) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ---------------------------------------------------------------
# 3) Flip CELIAC - CONTROL to CONTROL vs CELIAC
# ---------------------------------------------------------------
between_df <- Italy_vertebrate_contig_PA_contrast_results %>%
  filter(str_detect(contrast, "^CELIAC\\s*-\\s*CONTROL$"))

has_z  <- "z.ratio"  %in% names(between_df)
has_t  <- "t.ratio"  %in% names(between_df)
has_SE <- "SE"       %in% names(between_df)
has_l  <- "lower.CL" %in% names(between_df)
has_u  <- "upper.CL" %in% names(between_df)

between_tp_std <- between_df %>%
  mutate(
    est_flip   = -estimate,
    se_use     = if (has_SE) SE else NA_real_,
    stat_use   = if (has_z) z.ratio else if (has_t) t.ratio else NA_real_,
    low_flip   = if (has_l) -upper.CL else NA_real_,
    high_flip  = if (has_u) -lower.CL else NA_real_,
    OR         = exp(est_flip),
    term       = paste0("CONTROL vs CELIAC @ ", onset_timeline_combined),
    source     = "between_groups"
  ) %>%
  transmute(
    vertebrate_contig_id,
    term,
    estimate  = est_flip,
    std.error = se_use,
    statistic = stat_use,
    p.value,
    p.adj,
    conf.low  = low_flip,
    conf.high = high_flip,
    OR,
    source
  )

# ---------------------------------------------------------------
# 4) Harmonize with your main table
# ---------------------------------------------------------------
common_spec <- list(
  vertebrate_contig_id = NA_character_,
  term      = NA_character_,
  estimate  = NA_real_,
  std.error = NA_real_,
  statistic = NA_real_,
  p.value   = NA_real_,
  p.adj     = NA_real_,
  conf.low  = NA_real_,
  conf.high = NA_real_,
  OR        = NA_real_,
  source    = NA_character_
)

add_missing_cols_typed <- function(df, spec) {
  for (nm in names(spec)) if (!nm %in% names(df)) df[[nm]] <- spec[[nm]]
  df
}

main_tbl <- vertebrate_contig_Italy_PA.with_celiac

main_std <- main_tbl %>%
  add_missing_cols_typed(common_spec) %>%
  mutate(
    source = as.character(source),
    source = ifelse(is.na(source) | source == "", "model_or_emmeans", source)
  ) %>%
  select(names(common_spec))

if (nrow(between_tp_std) == 0) {
  between_std <- tibble(
    vertebrate_contig_id = character(), term = character(),
    estimate = numeric(), std.error = numeric(), statistic = numeric(),
    p.value = numeric(), p.adj = numeric(),
    conf.low = numeric(), conf.high = numeric(),
    OR = numeric(), source = character()
  )
} else {
  between_std <- between_tp_std %>%
    add_missing_cols_typed(common_spec) %>%
    mutate(
      source = as.character(source),
      source = ifelse(is.na(source) | source == "", "between_groups", source)
    ) %>%
    select(names(common_spec))
}

# ---------------------------------------------------------------
# 5) Final combined table
# ---------------------------------------------------------------
vertebrate_contig_Italy_PA.everything <- bind_rows(main_std, between_std) %>%
  arrange(vertebrate_contig_id, term)

dim(vertebrate_contig_Italy_PA.everything %>% filter(p.adj < 0.05))


write.csv(vertebrate_contig_Italy_PA.everything,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_PA_results_everything.csv",row.names = FALSE)



## -----------------------------------------------------------------------------

tp_levels <- c("t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")


vertebrate_contig_Italy_PA_heatmaps <- auto_heatmaps_detected(
    df            = vertebrate_contig_Italy_PA.everything,
    tp_levels     = tp_levels,
    fdr           = 0.05,
    scale         = "log",
    min_contigs   = 1,
    analysis_type = "PA"
  )


ggsave(vertebrate_contig_Italy_PA_heatmaps$`TP (CONTROL)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_PA_results_everything_TP_CONTROL.pdf",dpi = 600, width = 6, height = 6)


ggsave(vertebrate_contig_Italy_PA_heatmaps$`TP (CELIAC)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_PA_results_everything_TP_CELIAC.pdf",dpi = 600, width = 6, height = 2.5)

ggsave(vertebrate_contig_Italy_PA_heatmaps$`Between groups`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_PA_results_everything_Between_groups.pdf",dpi = 600, width = 6, height = 6)

ggsave(vertebrate_contig_Italy_PA_heatmaps$`Dx × TP interactions`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_PA_results_everything_interactions.pdf",dpi = 600, width = 6, height = 6)



## ----abundance_NB_total-------------------------------------------------------
library(edgeR)
library(glmmTMB)

# Prepare DGEList and normalization for total cohort
total.vertebrate_y <- DGEList(counts = total_vertebrate_contig_abundance_0.75_0.03, samples = total_vertebrate_metadata_0.75_0.03)
total.vertebrate_y <- calcNormFactors(total.vertebrate_y, method = "TMMwsp")
total.vertebrate_lib_eff <- with(total.vertebrate_y$samples, lib.size * norm.factors)

# Create long-format data for abundance modeling
vertebrate_contig_total_abundance_long <- total_vertebrate_contig_abundance_0.75_0.03 %>%
  as.data.frame() %>%
  rownames_to_column("vertebrate_contig_id") %>%
  pivot_longer(cols = -vertebrate_contig_id, names_to = "sample_id", values_to = "count") %>%
  left_join(total.vertebrate_y$samples %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  dplyr::mutate(offset = log(lib.size * norm.factors))

# Fit NB model for each contig
vertebrate_contig_total_abundance_nb_results <- vertebrate_contig_total_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(count ~ Dx.Status * onset_timeline_combined + Country + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(vertebrate_contig_id, tidy_results) %>%
  unnest(tidy_results)


vertebrate_contig_total_abundance.adjusted.results <- vertebrate_contig_total_abundance_nb_results %>%
  filter(!is.na(p.value)) %>%
  group_by(term) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

dim(vertebrate_contig_total_abundance.adjusted.results %>% filter(p.adj < 0.05))




## -----------------------------------------------------------------------------

library(emmeans)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(glmmTMB)

# 0) Define timepoint order (must include "t0" first)
tp_levels <- c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")
stopifnot("t0" %in% tp_levels)

# 1) Refit per-contig models (just for emmeans; your original table is untouched)
vertebrate_contig_total_abundance_models_for_emm <- vertebrate_contig_total_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(count ~ Dx.Status * onset_timeline_combined + Country + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x),
      error = function(e) NULL
    ))
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  select(vertebrate_contig_id, model)

# 2) Helper: get CELIAC TP vs t0 contrasts from emmeans (with SE/CI/p)
extract_celiac_tp_vs_t0 <- function(mod, tp_levels, dx_level = "CELIAC") {
  if (is.null(mod) || !("t0" %in% tp_levels)) return(tibble())

  # emmeans across timepoints within Dx
  emm <- tryCatch(
    emmeans(
      mod, ~ onset_timeline_combined | Dx.Status,
      at = list(onset_timeline_combined = tp_levels),
      type = "link"
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  # contrasts: each TP vs t0 (within each Dx)
  ref_idx <- which(tp_levels == "t0")[1]
  con <- tryCatch(
    contrast(emm, method = "trt.vs.ctrl", ref = ref_idx),
    error = function(e) NULL
  )
  if (is.null(con)) return(tibble())

  # IMPORTANT: use base::summary(), not emmeans::summary()
  s <- tibble::as_tibble(summary(con, infer = c(TRUE, TRUE)))

  # unify 'statistic' column (z.ratio OR t.ratio OR NA)
  if (!("statistic" %in% names(s))) {
    zr <- if ("z.ratio" %in% names(s)) s$z.ratio else NULL
    tr <- if ("t.ratio" %in% names(s)) s$t.ratio else NULL
    s$statistic <- dplyr::coalesce(zr, tr, rep(NA_real_, nrow(s)))
  }

  # ensure CI columns
  if (!all(c("lower.CL","upper.CL") %in% names(s))) {
    ci <- tryCatch(tibble::as_tibble(confint(con)), error = function(e) NULL)
    if (!is.null(ci)) {
      byk <- intersect(names(ci), names(s))
      s <- dplyr::left_join(
        s,
        dplyr::select(ci, dplyr::any_of(c(byk, "lower.CL","upper.CL","asymp.LCL","asymp.UCL"))),
        by = byk
      )
      if (!"lower.CL" %in% names(s) && "asymp.LCL" %in% names(s)) s$lower.CL <- s$asymp.LCL
      if (!"upper.CL" %in% names(s) && "asymp.UCL" %in% names(s)) s$upper.CL <- s$asymp.UCL
    }
  }
  if (!"lower.CL" %in% names(s)) s$lower.CL <- NA_real_
  if (!"upper.CL" %in% names(s)) s$upper.CL <- NA_real_

  s %>%
    dplyr::filter(Dx.Status == dx_level) %>%
    dplyr::mutate(
      TP        = stringr::str_remove(contrast, " - t0$"),
      term      = paste0("TP (", dx_level, "): ", TP, " vs t0"),
      OR        = exp(estimate)
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error = SE,
      statistic,
      p.value,
      conf.low  = lower.CL,
      conf.high = upper.CL,
      OR
    ) %>%
    dplyr::filter(!grepl("^TP \\(.+\\):\\s*t0\\s+vs\\s+t0$", term))
}


# --- Build the CELIAC rows and append to your original table -------------------
vertebrate_contig_total_abundance_celiac_tp_rows <- vertebrate_contig_total_abundance_models_for_emm %>%
  dplyr::mutate(out = purrr::map(model, ~ extract_celiac_tp_vs_t0(.x, tp_levels))) %>%
  dplyr::select(vertebrate_contig_id, out) %>%
  tidyr::unnest(out) %>%
  dplyr::group_by(term) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

vertebrate_contig_total_abundance.with_celiac <- dplyr::bind_rows(
  vertebrate_contig_total_abundance.adjusted.results,
  vertebrate_contig_total_abundance_celiac_tp_rows
)

dim(vertebrate_contig_total_abundance.with_celiac %>% filter(p.adj < 0.05))



## -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom.mixed)
library(emmeans)
library(tibble)

# ---------------------------------------------------------------
# 1) Fit models
# ---------------------------------------------------------------
total_vertebrate_contig_abundance_models <- vertebrate_contig_total_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(count ~ Dx.Status * onset_timeline_combined + Country + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x),
      error = function(e) NULL
    ))
  )

# ---------------------------------------------------------------
# 2) Extract between-group contrasts (CELIAC - CONTROL → flip)
# ---------------------------------------------------------------
extract_contrasts <- function(fit) {
  if (is.null(fit)) return(tibble())

  emm <- tryCatch(
    emmeans(fit, ~ Dx.Status | onset_timeline_combined),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  contr <- tryCatch(
    contrast(emm, method = "revpairwise", by = "onset_timeline_combined", adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(contr)) return(tibble())

  s  <- as_tibble(summary(contr))
  ci <- as_tibble(confint(contr)) %>%
        select(any_of(c("contrast","onset_timeline_combined","lower.CL","upper.CL")))

  keys <- intersect(names(s), names(ci))
  if (length(keys) > 0) {
    s <- left_join(s, ci, by = keys)
  }
  s
}

total_vertebrate_contig_abundance_contrast_results <- total_vertebrate_contig_abundance_models %>%
  mutate(contrasts = map(model, extract_contrasts)) %>%
  select(vertebrate_contig_id, contrasts) %>%
  tidyr::unnest(contrasts)

# adjust within each timepoint
total_vertebrate_contig_abundance_contrast_results <- total_vertebrate_contig_abundance_contrast_results %>%
  filter(!is.na(p.value)) %>%
  group_by(onset_timeline_combined) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ---------------------------------------------------------------
# 3) Flip CELIAC - CONTROL to CONTROL vs CELIAC
# ---------------------------------------------------------------
between_df <- total_vertebrate_contig_abundance_contrast_results %>%
  filter(str_detect(contrast, "^CELIAC\\s*-\\s*CONTROL$"))

has_z  <- "z.ratio"  %in% names(between_df)
has_t  <- "t.ratio"  %in% names(between_df)
has_SE <- "SE"       %in% names(between_df)
has_l  <- "lower.CL" %in% names(between_df)
has_u  <- "upper.CL" %in% names(between_df)

between_tp_std <- between_df %>%
  mutate(
    est_flip   = -estimate,
    se_use     = if (has_SE) SE else NA_real_,
    stat_use   = if (has_z) z.ratio else if (has_t) t.ratio else NA_real_,
    low_flip   = if (has_l) -upper.CL else NA_real_,
    high_flip  = if (has_u) -lower.CL else NA_real_,
    OR         = exp(est_flip),
    term       = paste0("CONTROL vs CELIAC @ ", onset_timeline_combined),
    source     = "between_groups"
  ) %>%
  transmute(
    vertebrate_contig_id,
    term,
    estimate  = est_flip,
    std.error = se_use,
    statistic = stat_use,
    p.value,
    p.adj,
    conf.low  = low_flip,
    conf.high = high_flip,
    OR,
    source
  )

# ---------------------------------------------------------------
# 4) Harmonize with your main table
# ---------------------------------------------------------------
common_spec <- list(
  vertebrate_contig_id = NA_character_,
  term      = NA_character_,
  estimate  = NA_real_,
  std.error = NA_real_,
  statistic = NA_real_,
  p.value   = NA_real_,
  p.adj     = NA_real_,
  conf.low  = NA_real_,
  conf.high = NA_real_,
  OR        = NA_real_,
  source    = NA_character_
)

add_missing_cols_typed <- function(df, spec) {
  for (nm in names(spec)) if (!nm %in% names(df)) df[[nm]] <- spec[[nm]]
  df
}

main_tbl <- vertebrate_contig_total_abundance.with_celiac

main_std <- main_tbl %>%
  add_missing_cols_typed(common_spec) %>%
  mutate(
    source = as.character(source),
    source = ifelse(is.na(source) | source == "", "model_or_emmeans", source)
  ) %>%
  select(names(common_spec))

if (nrow(between_tp_std) == 0) {
  between_std <- tibble(
    vertebrate_contig_id = character(), term = character(),
    estimate = numeric(), std.error = numeric(), statistic = numeric(),
    p.value = numeric(), p.adj = numeric(),
    conf.low = numeric(), conf.high = numeric(),
    OR = numeric(), source = character()
  )
} else {
  between_std <- between_tp_std %>%
    add_missing_cols_typed(common_spec) %>%
    mutate(
      source = as.character(source),
      source = ifelse(is.na(source) | source == "", "between_groups", source)
    ) %>%
    select(names(common_spec))
}

# ---------------------------------------------------------------
# 5) Final combined table
# ---------------------------------------------------------------
vertebrate_contig_total_abundance.everything <- bind_rows(main_std, between_std) %>%
  arrange(vertebrate_contig_id, term)

dim(vertebrate_contig_total_abundance.everything %>% filter(p.adj < 0.05))


write.csv(vertebrate_contig_total_abundance.everything,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_abundance_results_everything.csv",row.names = FALSE)



## -----------------------------------------------------------------------------

tp_levels <- c("t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")


vertebrate_contig_total_abundance_heatmaps <- auto_heatmaps_detected(
    df            = vertebrate_contig_total_abundance.everything,
    tp_levels     = tp_levels,
    fdr           = 0.05,
    scale         = "log",
    min_contigs   = 1,
    analysis_type = "abundance"
  )


ggsave(vertebrate_contig_total_abundance_heatmaps$`TP (CONTROL)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_abundance_results_everything_TP_CONTROL.pdf",dpi = 600, width = 6, height = 6)


ggsave(vertebrate_contig_total_abundance_heatmaps$`TP (CELIAC)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_abundance_results_everything_TP_CELIAC.pdf",dpi = 600, width = 6, height = 2.5)

ggsave(vertebrate_contig_total_abundance_heatmaps$`Between groups`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_abundance_results_everything_Between_groups.pdf",dpi = 600, width = 6, height = 6)

ggsave(vertebrate_contig_total_abundance_heatmaps$`Dx × TP interactions`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/total/total_vertebrate_contig_abundance_results_everything_interactions.pdf",dpi = 600, width = 6, height = 6)



## ----abundance_NB_US----------------------------------------------------------
library(edgeR)
library(glmmTMB)

# Prepare DGEList and normalization for US cohort
US.vertebrate_y <- DGEList(counts = US_vertebrate_contig_abundance_0.75_0.03, samples = US_vertebrate_metadata_0.75_0.03)
US.vertebrate_y <- calcNormFactors(US.vertebrate_y, method = "TMMwsp")
US.vertebrate_lib_eff <- with(US.vertebrate_y$samples, lib.size * norm.factors)

# Create long-format data for abundance modeling
vertebrate_contig_US_abundance_long <- US_vertebrate_contig_abundance_0.75_0.03 %>%
  as.data.frame() %>%
  rownames_to_column("vertebrate_contig_id") %>%
  pivot_longer(cols = -vertebrate_contig_id, names_to = "sample_id", values_to = "count") %>%
  left_join(US.vertebrate_y$samples %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  dplyr::mutate(offset = log(lib.size * norm.factors))

# Fit NB model for each contig
vertebrate_contig_US_abundance_nb_results <- vertebrate_contig_US_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(vertebrate_contig_id, tidy_results) %>%
  unnest(tidy_results)


vertebrate_contig_US_abundance.adjusted.results <- vertebrate_contig_US_abundance_nb_results %>%
  filter(!is.na(p.value)) %>%
  group_by(term) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

dim(vertebrate_contig_US_abundance.adjusted.results %>% filter(p.adj < 0.05))




## -----------------------------------------------------------------------------

library(emmeans)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(glmmTMB)

# 0) Define timepoint order (must include "t0" first)
tp_levels <- c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")
stopifnot("t0" %in% tp_levels)

# 1) Refit per-contig models (just for emmeans; your original table is untouched)
vertebrate_contig_US_abundance_models_for_emm <- vertebrate_contig_US_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x),
      error = function(e) NULL
    ))
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  select(vertebrate_contig_id, model)

# 2) Helper: get CELIAC TP vs t0 contrasts from emmeans (with SE/CI/p)
extract_celiac_tp_vs_t0 <- function(mod, tp_levels, dx_level = "CELIAC") {
  if (is.null(mod) || !("t0" %in% tp_levels)) return(tibble())

  # emmeans across timepoints within Dx
  emm <- tryCatch(
    emmeans(
      mod, ~ onset_timeline_combined | Dx.Status,
      at = list(onset_timeline_combined = tp_levels),
      type = "link"
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  # contrasts: each TP vs t0 (within each Dx)
  ref_idx <- which(tp_levels == "t0")[1]
  con <- tryCatch(
    contrast(emm, method = "trt.vs.ctrl", ref = ref_idx),
    error = function(e) NULL
  )
  if (is.null(con)) return(tibble())

  # IMPORTANT: use base::summary(), not emmeans::summary()
  s <- tibble::as_tibble(summary(con, infer = c(TRUE, TRUE)))

  # unify 'statistic' column (z.ratio OR t.ratio OR NA)
  if (!("statistic" %in% names(s))) {
    zr <- if ("z.ratio" %in% names(s)) s$z.ratio else NULL
    tr <- if ("t.ratio" %in% names(s)) s$t.ratio else NULL
    s$statistic <- dplyr::coalesce(zr, tr, rep(NA_real_, nrow(s)))
  }

  # ensure CI columns
  if (!all(c("lower.CL","upper.CL") %in% names(s))) {
    ci <- tryCatch(tibble::as_tibble(confint(con)), error = function(e) NULL)
    if (!is.null(ci)) {
      byk <- intersect(names(ci), names(s))
      s <- dplyr::left_join(
        s,
        dplyr::select(ci, dplyr::any_of(c(byk, "lower.CL","upper.CL","asymp.LCL","asymp.UCL"))),
        by = byk
      )
      if (!"lower.CL" %in% names(s) && "asymp.LCL" %in% names(s)) s$lower.CL <- s$asymp.LCL
      if (!"upper.CL" %in% names(s) && "asymp.UCL" %in% names(s)) s$upper.CL <- s$asymp.UCL
    }
  }
  if (!"lower.CL" %in% names(s)) s$lower.CL <- NA_real_
  if (!"upper.CL" %in% names(s)) s$upper.CL <- NA_real_

  s %>%
    dplyr::filter(Dx.Status == dx_level) %>%
    dplyr::mutate(
      TP        = stringr::str_remove(contrast, " - t0$"),
      term      = paste0("TP (", dx_level, "): ", TP, " vs t0"),
      OR        = exp(estimate)
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error = SE,
      statistic,
      p.value,
      conf.low  = lower.CL,
      conf.high = upper.CL,
      OR
    ) %>%
    dplyr::filter(!grepl("^TP \\(.+\\):\\s*t0\\s+vs\\s+t0$", term))
}


# --- Build the CELIAC rows and append to your original table -------------------
vertebrate_contig_US_abundance_celiac_tp_rows <- vertebrate_contig_US_abundance_models_for_emm %>%
  dplyr::mutate(out = purrr::map(model, ~ extract_celiac_tp_vs_t0(.x, tp_levels))) %>%
  dplyr::select(vertebrate_contig_id, out) %>%
  tidyr::unnest(out) %>%
  dplyr::group_by(term) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

vertebrate_contig_US_abundance.with_celiac <- dplyr::bind_rows(
  vertebrate_contig_US_abundance.adjusted.results,
  vertebrate_contig_US_abundance_celiac_tp_rows
)

dim(vertebrate_contig_US_abundance.with_celiac %>% filter(p.adj < 0.05))



## -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom.mixed)
library(emmeans)
library(tibble)

# ---------------------------------------------------------------
# 1) Fit models
# ---------------------------------------------------------------
US_vertebrate_contig_abundance_models <- vertebrate_contig_US_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x),
      error = function(e) NULL
    ))
  )

# ---------------------------------------------------------------
# 2) Extract between-group contrasts (CELIAC - CONTROL → flip)
# ---------------------------------------------------------------
extract_contrasts <- function(fit) {
  if (is.null(fit)) return(tibble())

  emm <- tryCatch(
    emmeans(fit, ~ Dx.Status | onset_timeline_combined),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  contr <- tryCatch(
    contrast(emm, method = "revpairwise", by = "onset_timeline_combined", adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(contr)) return(tibble())

  s  <- as_tibble(summary(contr))
  ci <- as_tibble(confint(contr)) %>%
        select(any_of(c("contrast","onset_timeline_combined","lower.CL","upper.CL")))

  keys <- intersect(names(s), names(ci))
  if (length(keys) > 0) {
    s <- left_join(s, ci, by = keys)
  }
  s
}

US_vertebrate_contig_abundance_contrast_results <- US_vertebrate_contig_abundance_models %>%
  mutate(contrasts = map(model, extract_contrasts)) %>%
  select(vertebrate_contig_id, contrasts) %>%
  tidyr::unnest(contrasts)

# adjust within each timepoint
US_vertebrate_contig_abundance_contrast_results <- US_vertebrate_contig_abundance_contrast_results %>%
  filter(!is.na(p.value)) %>%
  group_by(onset_timeline_combined) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ---------------------------------------------------------------
# 3) Flip CELIAC - CONTROL to CONTROL vs CELIAC
# ---------------------------------------------------------------
between_df <- US_vertebrate_contig_abundance_contrast_results %>%
  filter(str_detect(contrast, "^CELIAC\\s*-\\s*CONTROL$"))

has_z  <- "z.ratio"  %in% names(between_df)
has_t  <- "t.ratio"  %in% names(between_df)
has_SE <- "SE"       %in% names(between_df)
has_l  <- "lower.CL" %in% names(between_df)
has_u  <- "upper.CL" %in% names(between_df)

between_tp_std <- between_df %>%
  mutate(
    est_flip   = -estimate,
    se_use     = if (has_SE) SE else NA_real_,
    stat_use   = if (has_z) z.ratio else if (has_t) t.ratio else NA_real_,
    low_flip   = if (has_l) -upper.CL else NA_real_,
    high_flip  = if (has_u) -lower.CL else NA_real_,
    OR         = exp(est_flip),
    term       = paste0("CONTROL vs CELIAC @ ", onset_timeline_combined),
    source     = "between_groups"
  ) %>%
  transmute(
    vertebrate_contig_id,
    term,
    estimate  = est_flip,
    std.error = se_use,
    statistic = stat_use,
    p.value,
    p.adj,
    conf.low  = low_flip,
    conf.high = high_flip,
    OR,
    source
  )

# ---------------------------------------------------------------
# 4) Harmonize with your main table
# ---------------------------------------------------------------
common_spec <- list(
  vertebrate_contig_id = NA_character_,
  term      = NA_character_,
  estimate  = NA_real_,
  std.error = NA_real_,
  statistic = NA_real_,
  p.value   = NA_real_,
  p.adj     = NA_real_,
  conf.low  = NA_real_,
  conf.high = NA_real_,
  OR        = NA_real_,
  source    = NA_character_
)

add_missing_cols_typed <- function(df, spec) {
  for (nm in names(spec)) if (!nm %in% names(df)) df[[nm]] <- spec[[nm]]
  df
}

main_tbl <- vertebrate_contig_US_abundance.with_celiac

main_std <- main_tbl %>%
  add_missing_cols_typed(common_spec) %>%
  mutate(
    source = as.character(source),
    source = ifelse(is.na(source) | source == "", "model_or_emmeans", source)
  ) %>%
  select(names(common_spec))

if (nrow(between_tp_std) == 0) {
  between_std <- tibble(
    vertebrate_contig_id = character(), term = character(),
    estimate = numeric(), std.error = numeric(), statistic = numeric(),
    p.value = numeric(), p.adj = numeric(),
    conf.low = numeric(), conf.high = numeric(),
    OR = numeric(), source = character()
  )
} else {
  between_std <- between_tp_std %>%
    add_missing_cols_typed(common_spec) %>%
    mutate(
      source = as.character(source),
      source = ifelse(is.na(source) | source == "", "between_groups", source)
    ) %>%
    select(names(common_spec))
}

# ---------------------------------------------------------------
# 5) Final combined table
# ---------------------------------------------------------------
vertebrate_contig_US_abundance.everything <- bind_rows(main_std, between_std) %>%
  arrange(vertebrate_contig_id, term)

dim(vertebrate_contig_US_abundance.everything %>% filter(p.adj < 0.05))


write.csv(vertebrate_contig_US_abundance.everything,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_abundance_results_everything.csv",row.names = FALSE)



## -----------------------------------------------------------------------------

tp_levels <- c("t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")


vertebrate_contig_US_abundance_heatmaps <- auto_heatmaps_detected(
    df            = vertebrate_contig_US_abundance.everything,
    tp_levels     = tp_levels,
    fdr           = 0.05,
    scale         = "log",
    min_contigs   = 1,
    analysis_type = "abundance"
  )


ggsave(vertebrate_contig_US_abundance_heatmaps$`TP (CONTROL)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_abundance_results_everything_TP_CONTROL.pdf",dpi = 600, width = 6, height = 6)


ggsave(vertebrate_contig_US_abundance_heatmaps$`TP (CELIAC)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_abundance_results_everything_TP_CELIAC.pdf",dpi = 600, width = 6, height = 2.5)

ggsave(vertebrate_contig_US_abundance_heatmaps$`Between groups`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_abundance_results_everything_Between_groups.pdf",dpi = 600, width = 6, height = 6)

ggsave(vertebrate_contig_US_abundance_heatmaps$`Dx × TP interactions`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/US/US_vertebrate_contig_abundance_results_everything_interactions.pdf",dpi = 600, width = 6, height = 6)



## ----abundance_NB_Italy-------------------------------------------------------
library(edgeR)
library(glmmTMB)

# Prepare DGEList and normalization for Italy cohort
Italy.vertebrate_y <- DGEList(counts = Italy_vertebrate_contig_abundance_0.75_0.03, samples = Italy_vertebrate_metadata_0.75_0.03)
Italy.vertebrate_y <- calcNormFactors(Italy.vertebrate_y, method = "TMMwsp")
Italy.vertebrate_lib_eff <- with(Italy.vertebrate_y$samples, lib.size * norm.factors)

# Create long-format data for abundance modeling
vertebrate_contig_Italy_abundance_long <- Italy_vertebrate_contig_abundance_0.75_0.03 %>%
  as.data.frame() %>%
  rownames_to_column("vertebrate_contig_id") %>%
  pivot_longer(cols = -vertebrate_contig_id, names_to = "sample_id", values_to = "count") %>%
  left_join(Italy.vertebrate_y$samples %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  dplyr::mutate(offset = log(lib.size * norm.factors))

# Fit NB model for each contig
vertebrate_contig_Italy_abundance_nb_results <- vertebrate_contig_Italy_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ {
      tryCatch({
        glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x)
      }, error = function(e) NULL)
    }),
    tidy_results = map(model, ~ {
      if (!is.null(.x)) {
        tryCatch(tidy(.x, conf.int = TRUE), error = function(e) NULL)
      } else NULL
    })
  ) %>%
  select(vertebrate_contig_id, tidy_results) %>%
  unnest(tidy_results)


vertebrate_contig_Italy_abundance.adjusted.results <- vertebrate_contig_Italy_abundance_nb_results %>%
  filter(!is.na(p.value)) %>%
  group_by(term) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

dim(vertebrate_contig_Italy_abundance.adjusted.results %>% filter(p.adj < 0.05))




## -----------------------------------------------------------------------------

library(emmeans)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(glmmTMB)

# 0) Define timepoint order (must include "t0" first)
tp_levels <- c("t0","t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")
stopifnot("t0" %in% tp_levels)

# 1) Refit per-contig models (just for emmeans; your original table is untouched)
vertebrate_contig_Italy_abundance_models_for_emm <- vertebrate_contig_Italy_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x),
      error = function(e) NULL
    ))
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  select(vertebrate_contig_id, model)

# 2) Helper: get CELIAC TP vs t0 contrasts from emmeans (with SE/CI/p)
extract_celiac_tp_vs_t0 <- function(mod, tp_levels, dx_level = "CELIAC") {
  if (is.null(mod) || !("t0" %in% tp_levels)) return(tibble())

  # emmeans across timepoints within Dx
  emm <- tryCatch(
    emmeans(
      mod, ~ onset_timeline_combined | Dx.Status,
      at = list(onset_timeline_combined = tp_levels),
      type = "link"
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  # contrasts: each TP vs t0 (within each Dx)
  ref_idx <- which(tp_levels == "t0")[1]
  con <- tryCatch(
    contrast(emm, method = "trt.vs.ctrl", ref = ref_idx),
    error = function(e) NULL
  )
  if (is.null(con)) return(tibble())

  # IMPORTANT: use base::summary(), not emmeans::summary()
  s <- tibble::as_tibble(summary(con, infer = c(TRUE, TRUE)))

  # unify 'statistic' column (z.ratio OR t.ratio OR NA)
  if (!("statistic" %in% names(s))) {
    zr <- if ("z.ratio" %in% names(s)) s$z.ratio else NULL
    tr <- if ("t.ratio" %in% names(s)) s$t.ratio else NULL
    s$statistic <- dplyr::coalesce(zr, tr, rep(NA_real_, nrow(s)))
  }

  # ensure CI columns
  if (!all(c("lower.CL","upper.CL") %in% names(s))) {
    ci <- tryCatch(tibble::as_tibble(confint(con)), error = function(e) NULL)
    if (!is.null(ci)) {
      byk <- intersect(names(ci), names(s))
      s <- dplyr::left_join(
        s,
        dplyr::select(ci, dplyr::any_of(c(byk, "lower.CL","upper.CL","asymp.LCL","asymp.UCL"))),
        by = byk
      )
      if (!"lower.CL" %in% names(s) && "asymp.LCL" %in% names(s)) s$lower.CL <- s$asymp.LCL
      if (!"upper.CL" %in% names(s) && "asymp.UCL" %in% names(s)) s$upper.CL <- s$asymp.UCL
    }
  }
  if (!"lower.CL" %in% names(s)) s$lower.CL <- NA_real_
  if (!"upper.CL" %in% names(s)) s$upper.CL <- NA_real_

  s %>%
    dplyr::filter(Dx.Status == dx_level) %>%
    dplyr::mutate(
      TP        = stringr::str_remove(contrast, " - t0$"),
      term      = paste0("TP (", dx_level, "): ", TP, " vs t0"),
      OR        = exp(estimate)
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error = SE,
      statistic,
      p.value,
      conf.low  = lower.CL,
      conf.high = upper.CL,
      OR
    ) %>%
    dplyr::filter(!grepl("^TP \\(.+\\):\\s*t0\\s+vs\\s+t0$", term))
}


# --- Build the CELIAC rows and append to your original table -------------------
vertebrate_contig_Italy_abundance_celiac_tp_rows <- vertebrate_contig_Italy_abundance_models_for_emm %>%
  dplyr::mutate(out = purrr::map(model, ~ extract_celiac_tp_vs_t0(.x, tp_levels))) %>%
  dplyr::select(vertebrate_contig_id, out) %>%
  tidyr::unnest(out) %>%
  dplyr::group_by(term) %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

vertebrate_contig_Italy_abundance.with_celiac <- dplyr::bind_rows(
  vertebrate_contig_Italy_abundance.adjusted.results,
  vertebrate_contig_Italy_abundance_celiac_tp_rows
)

dim(vertebrate_contig_Italy_abundance.with_celiac %>% filter(p.adj < 0.05))



## -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom.mixed)
library(emmeans)
library(tibble)

# ---------------------------------------------------------------
# 1) Fit models
# ---------------------------------------------------------------
Italy_vertebrate_contig_abundance_models <- vertebrate_contig_Italy_abundance_long %>%
  group_by(vertebrate_contig_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      glmmTMB(count ~ Dx.Status * onset_timeline_combined + Sex + Age.at.Gluten.Introduction..months. + 
                HLA.Category + feeding_first_year + Delivery.Mode + (1 | patientID) + offset(offset),
                family = nbinom2, data = .x),
      error = function(e) NULL
    ))
  )

# ---------------------------------------------------------------
# 2) Extract between-group contrasts (CELIAC - CONTROL → flip)
# ---------------------------------------------------------------
extract_contrasts <- function(fit) {
  if (is.null(fit)) return(tibble())

  emm <- tryCatch(
    emmeans(fit, ~ Dx.Status | onset_timeline_combined),
    error = function(e) NULL
  )
  if (is.null(emm)) return(tibble())

  contr <- tryCatch(
    contrast(emm, method = "revpairwise", by = "onset_timeline_combined", adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(contr)) return(tibble())

  s  <- as_tibble(summary(contr))
  ci <- as_tibble(confint(contr)) %>%
        select(any_of(c("contrast","onset_timeline_combined","lower.CL","upper.CL")))

  keys <- intersect(names(s), names(ci))
  if (length(keys) > 0) {
    s <- left_join(s, ci, by = keys)
  }
  s
}

Italy_vertebrate_contig_abundance_contrast_results <- Italy_vertebrate_contig_abundance_models %>%
  mutate(contrasts = map(model, extract_contrasts)) %>%
  select(vertebrate_contig_id, contrasts) %>%
  tidyr::unnest(contrasts)

# adjust within each timepoint
Italy_vertebrate_contig_abundance_contrast_results <- Italy_vertebrate_contig_abundance_contrast_results %>%
  filter(!is.na(p.value)) %>%
  group_by(onset_timeline_combined) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ---------------------------------------------------------------
# 3) Flip CELIAC - CONTROL to CONTROL vs CELIAC
# ---------------------------------------------------------------
between_df <- Italy_vertebrate_contig_abundance_contrast_results %>%
  filter(str_detect(contrast, "^CELIAC\\s*-\\s*CONTROL$"))

has_z  <- "z.ratio"  %in% names(between_df)
has_t  <- "t.ratio"  %in% names(between_df)
has_SE <- "SE"       %in% names(between_df)
has_l  <- "lower.CL" %in% names(between_df)
has_u  <- "upper.CL" %in% names(between_df)

between_tp_std <- between_df %>%
  mutate(
    est_flip   = -estimate,
    se_use     = if (has_SE) SE else NA_real_,
    stat_use   = if (has_z) z.ratio else if (has_t) t.ratio else NA_real_,
    low_flip   = if (has_l) -upper.CL else NA_real_,
    high_flip  = if (has_u) -lower.CL else NA_real_,
    OR         = exp(est_flip),
    term       = paste0("CONTROL vs CELIAC @ ", onset_timeline_combined),
    source     = "between_groups"
  ) %>%
  transmute(
    vertebrate_contig_id,
    term,
    estimate  = est_flip,
    std.error = se_use,
    statistic = stat_use,
    p.value,
    p.adj,
    conf.low  = low_flip,
    conf.high = high_flip,
    OR,
    source
  )

# ---------------------------------------------------------------
# 4) Harmonize with your main table
# ---------------------------------------------------------------
common_spec <- list(
  vertebrate_contig_id = NA_character_,
  term      = NA_character_,
  estimate  = NA_real_,
  std.error = NA_real_,
  statistic = NA_real_,
  p.value   = NA_real_,
  p.adj     = NA_real_,
  conf.low  = NA_real_,
  conf.high = NA_real_,
  OR        = NA_real_,
  source    = NA_character_
)

add_missing_cols_typed <- function(df, spec) {
  for (nm in names(spec)) if (!nm %in% names(df)) df[[nm]] <- spec[[nm]]
  df
}

main_tbl <- vertebrate_contig_Italy_abundance.with_celiac

main_std <- main_tbl %>%
  add_missing_cols_typed(common_spec) %>%
  mutate(
    source = as.character(source),
    source = ifelse(is.na(source) | source == "", "model_or_emmeans", source)
  ) %>%
  select(names(common_spec))

if (nrow(between_tp_std) == 0) {
  between_std <- tibble(
    vertebrate_contig_id = character(), term = character(),
    estimate = numeric(), std.error = numeric(), statistic = numeric(),
    p.value = numeric(), p.adj = numeric(),
    conf.low = numeric(), conf.high = numeric(),
    OR = numeric(), source = character()
  )
} else {
  between_std <- between_tp_std %>%
    add_missing_cols_typed(common_spec) %>%
    mutate(
      source = as.character(source),
      source = ifelse(is.na(source) | source == "", "between_groups", source)
    ) %>%
    select(names(common_spec))
}

# ---------------------------------------------------------------
# 5) Final combined table
# ---------------------------------------------------------------
vertebrate_contig_Italy_abundance.everything <- bind_rows(main_std, between_std) %>%
  arrange(vertebrate_contig_id, term)

dim(vertebrate_contig_Italy_abundance.everything %>% filter(p.adj < 0.05))


write.csv(vertebrate_contig_Italy_abundance.everything,"~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_abundance_results_everything.csv",row.names = FALSE)



## -----------------------------------------------------------------------------

tp_levels <- c("t0-6","t0-12","t0-18","t0-24","t0-30","t0-36","t0-over42")


vertebrate_contig_Italy_abundance_heatmaps <- auto_heatmaps_detected(
    df            = vertebrate_contig_Italy_abundance.everything,
    tp_levels     = tp_levels,
    fdr           = 0.05,
    scale         = "log",
    min_contigs   = 1,
    analysis_type = "abundance"
  )


ggsave(vertebrate_contig_Italy_abundance_heatmaps$`TP (CONTROL)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_abundance_results_everything_TP_CONTROL.pdf",dpi = 600, width = 6, height = 6)


ggsave(vertebrate_contig_Italy_abundance_heatmaps$`TP (CELIAC)`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_abundance_results_everything_TP_CELIAC.pdf",dpi = 600, width = 6, height = 2.5)

ggsave(vertebrate_contig_Italy_abundance_heatmaps$`Between groups`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_abundance_results_everything_Between_groups.pdf",dpi = 600, width = 6, height = 6)

ggsave(vertebrate_contig_Italy_abundance_heatmaps$`Dx × TP interactions`,file = "~/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs/Italy/Italy_vertebrate_contig_abundance_results_everything_interactions.pdf",dpi = 600, width = 6, height = 6)


