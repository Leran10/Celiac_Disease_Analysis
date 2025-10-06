# Section: load libraries... ----
library(ggplot2)
library(dada2)
packageVersion("dada2")
library(tidyverse)
packageVersion("tidyverse")
library(ggpubr)
packageVersion("ggpubr")
library(phyloseq)
packageVersion("phyloseq")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(DECIPHER)
packageVersion("DECIPHER")
library(phangorn)
packageVersion("phangorn")
library(openxlsx)
packageVersion("openxlsx")

# Set global ggplot settings
theme_set(theme_bw(base_size = 10))


# Section: Utility functions... ----
# These functions are used throught out the code in this file.
get_samples_ids <- function (f_name) {
  # Get the sample id from the file name; these ids are
  # unique and they have to match the metadata file

  # Remove the extension from the file name
  f_base_name <- tools::file_path_sans_ext(f_name)
  if (startsWith(f_base_name, "M951_Leonard_CELIACS")){
    # Extract parts of the base name
    parts <- strsplit(f_base_name, "_")[[1]]
    # Extract parts of the sample id
    country <- parts[4]
    gmm <- parts[5]
    patient_id <- parts[6]
    age <- parts[7]
    if (endsWith(age, "Y")){
      age <- paste0(as.numeric(substr(age, 1, nchar(age) - 1)) * 12, "M")
    }
    # assemble sample id
    sample_id <- paste0(country, "_", gmm, "_", patient_id, "_", age)
  } else if (startsWith(f_base_name, "M951_Water_Water")){
    sample_id <- strsplit(f_base_name, "_")[[1]][5]
  } else {
    stop("Unrecognized filename pattern: ", f_base_name)
  }
  sample_id
}


# Section: set paths... ----
cat("=== Section: set paths and variables .. ===\n")
cat("setting paths..\n")
p_wd <- Sys.getenv("PROJ_CELIAC_PATH", "/scratch/sahlab/dabid/proj_celiac")
p_wd_data <- file.path(p_wd, "data")
p_wd_data_miseq <- file.path(p_wd_data, "MiSeq_M951_16S")
p_annot_gg2 <- file.path(p_wd_data, "annotations_16S", "2022.10.backbone.full-length.dada2.train.fna")  # recommended by Scott
p_wd_output <- file.path(p_wd, "output")
p_ps_gg2 = file.path(p_wd_output, "ps_gg2.RDS")

cat("set variables..\n")
cpus_per_task <- Sys.getenv("SLURM_CPUS_PER_TASK", "5")
cpus_per_task <- as.numeric(cpus_per_task)
cat("# of cpus per task is: ", cpus_per_task, "\n")

# Section: read and filter forward and reverse reads ----
# Forward and reverse fastq f_names have format: 
# SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz
cat("=== Section: read & filter F and R reads .. ===\n")
# get paths of forward and reverse reads
l_path_for <- sort(list.files(p_wd_data_miseq, pattern = "R1", full.names = TRUE))
l_path_rev <- sort(list.files(p_wd_data_miseq, pattern = "R2", full.names = TRUE))
# get names of forwards and reverse reads
l_name_for <- sort(list.files(p_wd_data_miseq, pattern = "R1", full.names = FALSE))
l_name_rev <- sort(list.files(p_wd_data_miseq, pattern = "R2", full.names = FALSE))

# get unique samples ids from file names
l_sample_ids_for <- mapply(get_samples_ids, l_name_for)
l_sample_ids_rev <- mapply(get_samples_ids, l_name_rev)

cat("# of samples in forward reads: ", length(l_sample_ids_for), "\n")
cat("# of samples in reverse reads: ", length(l_sample_ids_rev), "\n")

# create filtered directory if it doesn't exist
p_wd_data_miseq_filtered <- file.path(p_wd_data_miseq, "filtered")
if (!dir.exists(p_wd_data_miseq_filtered)) {
  dir.create(p_wd_data_miseq_filtered, recursive = TRUE)
}
# set paths for filtered reads
l_path_for_fltr <- file.path(p_wd_data_miseq, "filtered", paste0(l_sample_ids_for, "_F_fltr.fastq"))
l_path_rev_fltr <- file.path(p_wd_data_miseq, "filtered", paste0(l_sample_ids_rev, "_R_fltr.fastq"))

# filter and trim reads
qc <- filterAndTrim(fwd = l_path_for,
                    filt = l_path_for_fltr,
                    rev = l_path_rev,
                    filt.rev = l_path_rev_fltr,
                    trimLeft = c(0, 0),
                    truncLen = c(250, 230),
                    maxEE = 2,
                    maxN = 0,
                    rm.phix = TRUE,
                    compress = TRUE,
                    verbose = TRUE,
                    multithread = cpus_per_task)
df_qc = as.data.frame(qc)
cat("total reads:", nrow(df_qc), "\n")
cat("# of samples with 100% retention:", sum(df_qc$reads.out == df_qc$reads.in), "\n")

# Identify samples with no reads after filtering
samples_with_no_reads <- rownames(df_qc)[df_qc$reads.out == 0]
cat("Samples with no reads after filtering:", length(samples_with_no_reads), "\n")
if(length(samples_with_no_reads) > 0) {
  cat("Sample names with no reads:\n")
  print(samples_with_no_reads)
}

# Create filtered lists excluding samples with no reads
l_kept_samples <- df_qc$reads.out > 0
l_path_for_fltr_kept <- l_path_for_fltr[l_kept_samples]
l_path_rev_fltr_kept <- l_path_rev_fltr[l_kept_samples]
l_sample_ids_for_kept <- l_sample_ids_for[l_kept_samples]
l_sample_ids_rev_kept <- l_sample_ids_rev[l_kept_samples]

cat("Samples kept forward reads:", length(l_sample_ids_for_kept), "\n")
cat("Samples kept reverse reads:", length(l_sample_ids_rev_kept), "\n")


# Section: dereplicate, denoise, merge reads & create ASV table... ----
cat("=== Section: dereplicate, denoise, merge reads & create ASV table .. ===\n")
cat("dereplicate..\n")
# dereplicate: reduce the number of reads to the unique sequences
l_path_for_fltr_derep <- derepFastq(l_path_for_fltr_kept)
l_path_rev_fltr_derep <- derepFastq(l_path_rev_fltr_kept)
# name them by the sample names
names(l_path_for_fltr_derep) <- basename(l_path_for_fltr_kept)
names(l_path_rev_fltr_derep) <- basename(l_path_rev_fltr_kept)
cat("denoise..\n")
# learn the error rates for the forward and reverse reads
error_rate_for <- learnErrors(l_path_for_fltr_kept, multithread = cpus_per_task) 
error_rate_rev <- learnErrors(l_path_rev_fltr_kept, multithread = cpus_per_task)
cat("run dada2..\n")
# denoise, identify ASVs, and estimate abundance
l_path_for_fltr_derep_denoise <- dada(l_path_for_fltr_derep,
                                      err = error_rate_for, multithread = cpus_per_task)
l_path_rev_fltr_derep_denoise <- dada(l_path_rev_fltr_derep,
                                       err = error_rate_rev, multithread = cpus_per_task)
cat("merge reads..\n")
merged_reads <- mergePairs(l_path_for_fltr_derep_denoise,
                           l_path_for_fltr_derep,
                           l_path_rev_fltr_derep_denoise,
                           l_path_rev_fltr_derep,
                           verbose=TRUE)
cat("create ASV table..\n")
tab_asv <- makeSequenceTable(merged_reads)
cat("ASV table dimensions: ", dim(tab_asv), "\n")


# Section: remove chimeras... ----
cat("=== Section: remove chimeras .. ===\n")
tab_asv_no_chimeras <- removeBimeraDenovo(tab_asv, method = "consensus", verbose = TRUE)
cat("total reads: ", ncol(tab_asv), "\n")
cat("total reads after removing chimeras: ", ncol(tab_asv_no_chimeras), "\n")


# Section: Summary of reads across steps: filtering, de-noising,and chimera removal... ----
cat("=== Section: summary of reads across steps.. ===\n")
# summary 1: all samples (filtering stage only using filterAndTrim command)
df_summary_fltr <- data.frame(
  sample_id = l_sample_ids_for,
  total_reads = df_qc$reads.in,
  filtered_reads = df_qc$reads.out,
  filtering_retention = df_qc$reads.out / df_qc$reads.in,
  status = ifelse(df_qc$reads.out > 0, "kept", "filtered_out"))

# summary 2: kept samples only (all stages)
df_summary_kept <- data.frame(sample_id = l_sample_ids_for_kept,
                              total_reads = df_qc$reads.in[l_kept_samples],
                              filtered_reads = df_qc$reads.out[l_kept_samples],
                              denoised_reads_f = sapply(l_path_for_fltr_derep_denoise,
                                                        function(x) sum(getUniques(x))),
                              denoised_reads_r = sapply(l_path_rev_fltr_derep_denoise,
                                                        function(x) sum(getUniques(x))),
                              merged = sapply(merged_reads, function(x) sum(getUniques(x))),
                              no_chimeras = rowSums(tab_asv_no_chimeras))

# Calculate retention rates for kept samples
df_summary_kept <- df_summary_kept %>%
  mutate(filtering_retention = filtered_reads / total_reads,
         denoising_retention_f = denoised_reads_f / filtered_reads,
         denoising_retention_r = denoised_reads_r / filtered_reads,
         merging_retention = merged / filtered_reads,
         chimera_retention = ifelse(merged > 0, no_chimeras / merged, NA_real_),
         overall_retention = no_chimeras / total_reads)

# Display summary statistics
cat("Total samples processed:", nrow(df_summary_fltr), "\n")
cat("Samples kept after filtering:", sum(df_summary_fltr$status == "kept"), "\n")
cat("Samples filtered out:", sum(df_summary_fltr$status == "filtered_out"), "\n")
cat("Mean filtering retention (all samples):",
  round(mean(df_summary_fltr$filtering_retention, na.rm = TRUE), 3), "\n")
cat("Mean filtering retention (kept samples):",
  round(mean(df_summary_kept$filtering_retention), 3), "\n")
cat("Mean denoising retention forward reads:",
  round(mean(df_summary_kept$denoising_retention_f), 3), "\n")
cat("Mean denoising retention reverse reads:",
  round(mean(df_summary_kept$denoising_retention_r), 3), "\n")
cat("Mean merging retention:", round(mean(df_summary_kept$merging_retention), 3), "\n")
cat("Mean chimera retention:", round(mean(df_summary_kept$chimera_retention), 3), "\n")
cat("Mean overall retention:", round(mean(df_summary_kept$overall_retention), 3), "\n")

# Write summaries to Excel file
excel_file <- file.path(p_wd_output, "reads_summary.xlsx")
cat("Saving read summaries to:", excel_file, "\n")

# Create workbook and add sheets
wb <- createWorkbook()
addWorksheet(wb, "all_samples")
addWorksheet(wb, "kept_samples")

# Write data to sheets
writeData(wb, "all_samples", df_summary_fltr)
writeData(wb, "kept_samples", df_summary_kept)

# Save workbook
saveWorkbook(wb, excel_file, overwrite = TRUE)


# Section: assign taxonomy... ----
cat("=== Section: assign taxonomy and create a phyloseq object .. ===\n")
taxa_gg2 <- assignTaxonomy(tab_asv_no_chimeras, p_annot_gg2, multithread = TRUE)
taxa_gg2 <- apply(taxa_gg2, 2, function(x) gsub("[dpcofgs]__", "", x))
cat("create a phyloseq object..\n")
ps_gg2 <- phyloseq(otu_table(tab_asv_no_chimeras, taxa_are_rows = FALSE), tax_table(taxa_gg2))

# Clean sample names by removing _F_fltr.fastq suffix
cat("Cleaning sample names...\n")
sample_names(ps_gg2) <- gsub("_F_fltr\\.fastq$", "", sample_names(ps_gg2))
cat("Sample names after cleaning:\n")
print(head(sample_names(ps_gg2)))


# Section: Phylogenetic tree... ----
cat("=== Section: create a phylogenetic tree .. ===\n")
# extract asv sequences adn set names to sequences (matching dada2 columns)
l_seq <- getSequences(tab_asv_no_chimeras)
names(l_seq) <- l_seq
# align using DECIPHER; fast and robust
align <- AlignSeqs(DNAStringSet(l_seq), anchor = NA)
# build ML tree under GTR model
align_phydat <- phyDat(as.matrix(align), type = "DNA")
mat_dist <- dist.ml(align_phydat)
tree_nj <- NJ(mat_dist)
fit <- pml(tree_nj, data = align_phydat)
fit_gtr <- optim.pml(update(fit, k = 4, inv = 0.2),
                     model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))
tree <- phangorn::midpoint(fit_gtr$tree)  # optional but helpful
# attach the tree to the phyloseq object
stopifnot(identical(taxa_names(ps_gg2), names(l_seq)))
phy_tree(ps_gg2) <- tree
ps_gg2
# save the phyloseq object
saveRDS(ps_gg2, file = p_ps_gg2)