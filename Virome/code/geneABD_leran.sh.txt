#!/bin/bash

#Remove redundancy of the genes

# --mem=300GB
# --cpus-per-task=24
# --time=300:00:00
# --output=geneABD_workflow.o%j

source /ref/sahlab/software/miniforge3/bin/activate
conda activate pholdENV

# phold proteins-predict -i ${WD}/crassvirales_linclust_095_rep_seq.fasta \
# -o ${WD}/crassvirales_pholdPredict --cpu -t 24 \
# -d /ref/sahlab/data/protein_Structure_DBs/phold

# conda activate mmseqs

# GENESEQS=/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/post_assembly_real/14_pharokka_phages_picobirna/prodigal-gv.ffn
# WD=/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/post_assembly_real/18_genAbdEstimation

# mmseqs easy-linclust ${GENESEQS} \
# ${WD}/linclust_099 \
# /scratch/sahlab/Luis/temporal/ \
# --min-seq-id 0.99 \
# --cov-mode 0 \
# -c 0.9 \
# --threads 24
# #--cov-mode 0: coverage of query and target


#Mapping non-redundant gene file 

#SBATCH --mem=50GB
#SBATCH --cpus-per-task=24
#SBATCH --time=300:00:00
#SBATCH --array=1-1150%20
#SBATCH--output=geneABD_workflow_%A_%a.out

ID=$( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/post_assembly_real/lookup_mappingForCoverage.txt )

source /ref/sahlab/software/anaconda3/bin/activate
conda activate coverm
WD=/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/post_assembly_real/18_genAbdEstimation
OUT=${WD}/mappingStats

coverm contig -1 /lts/sahlab/data4/luis/RC2_analysisProgress/01_cleanAllFastqForMappingCov/${ID}_R1.fastq \
-2 /lts/sahlab/data4/luis/RC2_analysisProgress/01_cleanAllFastqForMappingCov/${ID}_R2.fastq \
-r ${WD}/linclust_099_rep_seq.fasta \
--mapper minimap2-sr \
--threads ${SLURM_CPUS_PER_TASK} \
--methods rpkm count variance mean covered_fraction covered_bases > ${OUT}/${ID}_stats.txt

#merge all stats 
#for file in *.txt; do NAME=$(basename ${file} .txt);  tail -n +2 "$file" | awk -v fname="$NAME" '{print fname "\t" $0}' >> ../Coverm_concatenated.txt; done
