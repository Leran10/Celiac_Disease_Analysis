#!/bin/bash

#SBATCH --mem=100GB
#SBATCH --cpus-per-task=24
#SBATCH --time=500:00:00
#SBATCH --output=phage_New_WF.o%j 

#Exit immediately if a command exits with a non-zero status
set -e

#Activate environment
source /ref/sahlab/software/miniforge3/bin/activate

#Define main paths for easier maintenance
BASE_DIR="/scratch/sahlab/Luis/hecOUT_Celiac_total_data"
SCRIPT_DIR="/home/luisalberto/Scripts"
REF_DIR="/ref/sahlab/data"

#Input files and directories
HECATOMBASSEMBLY="${BASE_DIR}/processing/assembly/CONTIG_DICTIONARY/FLYE/assembly_graph.gfa"
READS="${BASE_DIR}/processing/Reads"
DB="${REF_DIR}/nr/nr_mmseqs15_DB"
GENOMAD_DB="${REF_DIR}/genomad_db"
PHOLD_DB="${REF_DIR}/protein_Structure_DBs/phold"

#Output directories
RENEO_OUT="${BASE_DIR}/postHec/01_reneo"
MMSEQS_OUT="${BASE_DIR}/postHec/02_mmseqs"
TMP_DIR="${BASE_DIR}/temp"
FILTERING_WD="${BASE_DIR}/postHec/03_selectingViruses"
PHAGEPRED_WD="${BASE_DIR}/postHec/04_phagePred"

Step 1: Binning using reneo. 
IMPORTANT: For large datasets, make sure to have the reneo and Koverage config files in the directory with proper resources. 
conda activate reneo
reneo run --input "${HECATOMBASSEMBLY}" \
          --reads "${READS}" \
          --minlength 1000 \
          --output "${RENEO_OUT}" \
          --threads ${SLURM_CPUS_PER_TASK}

Step 2: mmseqs LCA for selecting root, NAs, and viral contigs
mkdir -p "${MMSEQS_OUT}"
conda activate mmseqs2_v15-6f452
mmseqs easy-taxonomy "${RENEO_OUT}/genomes_and_unresolved_edges.fasta" \
       "${DB}" \
       "${MMSEQS_OUT}/genomes_and_unresolved_edges_mmseqs" \
       "${TMP_DIR}" \
       --min-length 30 \
       -e 1e-15 \
       --search-type 2 \
       -s 4.0 \
       --shuffle 0 \
       --lca-mode 2  \
       -a \
       --tax-lineage 2 \
       --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
       --threads ${SLURM_CPUS_PER_TASK} \
       --split-mode 0 \
       --orf-filter 1

Step 3: Selecting contigs and generating viral fasta file
mkdir -p "${FILTERING_WD}"

python "${SCRIPT_DIR}/filter_mmseqsLCA.py" --mmseqs_LCA_table "${MMSEQS_OUT}/genomes_and_unresolved_edges_mmseqs_lca.tsv" \
       --o_filtered_LCA_table "${FILTERING_WD}/filtered_output.txt" \
       --o_passing_contig_ids "${FILTERING_WD}/passing_contig_ids.txt" \
       --contigs "${RENEO_OUT}/genomes_and_unresolved_edges.fasta"

#Step 4: Get fasta for the contigs
seqkit grep -f "${FILTERING_WD}/passing_contig_ids.txt" "${RENEO_OUT}/genomes_and_unresolved_edges.fasta" > "${FILTERING_WD}/passing_Viralcontigs.fasta"

#Step 5: Phage prediction
mkdir -p "${PHAGEPRED_WD}"

conda activate jaeger
Jaeger -i "${FILTERING_WD}/passing_Viralcontigs.fasta"  \
       -o "${PHAGEPRED_WD}/jaeger" \
       -s 2.5 \
       --fsize 1000 \
       --stride 1000

conda activate /ref/sahlab/software/anaconda3/envs/genomad
genomad end-to-end --min-score 0.6 \
       --cleanup \
       --threads ${SLURM_CPUS_PER_TASK} \
       "${FILTERING_WD}/passing_Viralcontigs.fasta" \
       "${PHAGEPRED_WD}/geNomad" \
       "${GENOMAD_DB}"

#Step 6: Phold from phage prediction
conda activate pholdENV
phold run -i "${FILTERING_WD}/passing_Viralcontigs.fasta" \
       -o "${PHAGEPRED_WD}/phold" \
       -d "${PHOLD_DB}" \
       -t ${SLURM_CPUS_PER_TASK} --cpu --force