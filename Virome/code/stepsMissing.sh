
#Step 7: Run checkV on the passing_ViralContigs file
# mkdir -p "${QUALPRED_WD}"

# conda activate /ref/sahlab/software/anaconda3/envs/checkV 
# checkv end_to_end "${FILTERING_WD}/passing_Viralcontigs.fasta" \
#        "${QUALPRED_WD}/passing_Viralcontigs_checkV" \
#        -d /ref/sahlab/data/viral_analysis_DBs/checkV_DB/checkv-db-v1.5 \
#        -t ${SLURM_CPUS_PER_TASK}

# echo "checkV passing_ViralContigs done" 

#step 8: do the parsing using all the prediction informations and functional annotation and create the final phage set
mkdir -p "${PHAGE_INFO_WD}"
conda activate R

Rscript "${SCRIPT_DIR}/02_phagePrediction.R" \
"${PHAGEPRED_WD}/phold/pholdCompare/phold_per_cds_predictions.tsv" \
"${PHAGEPRED_WD}/jaeger/passing_Viralcontigs_default_jaeger.tsv" \
"${PHAGEPRED_WD}/geNomad/passing_Viralcontigs_summary/passing_Viralcontigs_virus_summary.tsv" \
"${QUALPRED_WD}/passing_Viralcontigs_checkV/quality_summary.tsv" \
"${PHAGE_INFO_WD}"

seqkit grep -f "${PHAGE_INFO_WD}/contig_ids.txt" "${FILTERING_WD}/passing_Viralcontigs.fasta" > "${PHAGE_INFO_WD}/phageContigs.fasta"

echo "phageContigs identification done" 

#final checkV checking 
conda activate /ref/sahlab/software/anaconda3/envs/checkV 
checkv end_to_end "${PHAGE_INFO_WD}/phageContigs.fasta" \
       "${QUALPRED_WD}/phageContigs_checkV" \
       -d /ref/sahlab/data/viral_analysis_DBs/checkV_DB/checkv-db-v1.5 \
       -t ${SLURM_CPUS_PER_TASK}

echo "checkV phageContigs done" 