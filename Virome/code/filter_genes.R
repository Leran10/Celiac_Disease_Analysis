#Part 1: Load, filter, and save the wide format table
filter_and_save_wide_data <- function(genesCoverm, 
                                      covered_fraction_threshold = 0.75, 
                                      output_path = "~/Handley Lab Dropbox/Manuscripts/RC2/phage/02_phageMicrobialAnalysis/geneBasedAnalysis/results") {
  
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
