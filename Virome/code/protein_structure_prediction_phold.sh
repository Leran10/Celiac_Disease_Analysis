#!/bin/bash

#SBATCH --mem=100GB
#SBATCH --cpus-per-task=24
#SBATCH --time=300:00:00
#SBATCH --output=geneABD_workflow.o%j

source /ref/sahlab/software/miniforge3/bin/activate

# #Get gene sequence file from the vOTU_repSeqs.fasta file 
# conda activate pharokka_env

VOTU=/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac/06_phages_cluster/vOTU_repSeqs.fasta
GENBASEDDIR=/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac/09_genAbdEstimation
PROTSEQS=${GENBASEDDIR}/vOTU_repSeqs_proteins.faa
GENESEQS=${GENBASEDDIR}/vOTU_repSeqs_genes.ffn

# mkdir -p "${GENBASEDDIR}"
# mkdir -p "${GENBASEDDIR}/mappingStats"

# pyrodigal -a ${PROTSEQS} -d ${GENESEQS} -i ${VOTU} -j ${SLURM_CPUS_PER_TASK} -p meta

#Phold protein annotation

conda activate pholdENV

phold proteins-predict -i ${PROTSEQS} \
-o ${GENBASEDDIR}/vOTU_repSeqs_proteins_pholdPredict \
--cpu \
-t ${SLURM_CPUS_PER_TASK} \
-d /ref/sahlab/data/protein_Structure_DBs/phold

phold proteins-compare -i ${PROTSEQS} \
--predictions_dir ${GENBASEDDIR}/vOTU_repSeqs_proteins_pholdPredict \
-o ${GENBASEDDIR}/vOTU_repSeqs_proteins_pholdCompare \
-t ${SLURM_CPUS_PER_TASK} \
-d /ref/sahlab/data/protein_Structure_DBs/phold

# #Clustering genes by identity to remove redundacy 
# conda activate mmseqs2_v15-6f452

# mmseqs easy-linclust ${GENESEQS} \
# ${GENBASEDDIR}/linclust_099 \
# /scratch/sahlab/Luis/temporal/ \
# --min-seq-id 0.99 \
# --cov-mode 0 \
# -c 0.9 \
# --threads 24
#--cov-mode 0: coverage of query and target


#Mapping non-redundant gene file 

# --mem=20GB
# --cpus-per-task=6
# --time=300:00:00
# --array=1-1150%40
# --output=geneABD_workflow_%A_%a.out

# ID=$( sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/post_assembly_real/lookup_mappingForCoverage.txt )

# source /ref/sahlab/software/anaconda3/bin/activate
# conda activate coverm

# GENBASEDDIR=/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac/09_genAbdEstimation
# READS=/lts/sahlab/data4/luis/RC2_analysisProgress/01_cleanAllFastqForMappingCov
# OUT=${GENBASEDDIR}/mappingStats

# coverm contig -1 ${READS}/${ID}_R1.fastq \
# -2 ${READS}/${ID}_R2.fastq \
# -r ${GENBASEDDIR}/linclust_099_rep_seq.fasta \
# --mapper minimap2-sr \
# --threads ${SLURM_CPUS_PER_TASK} \
# --methods rpkm count variance mean covered_fraction covered_bases > ${OUT}/${ID}_stats.txt

#for file in *.txt; do NAME=$(basename ${file} .txt);  tail -n +2 "$file" | awk -v fname="$NAME" '{print fname "\t" $0}' >> ../Coverm_concatenated.txt; done