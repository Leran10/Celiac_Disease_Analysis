#!/bin/bash
#SBATCH --mem=40G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH -o /scratch/sahlab/dabid/slurm_logs/proj_celiac/00_create_phyloseq_obj_%J.out
#SBATCH -e /scratch/sahlab/dabid/slurm_logs/proj_celiac/00_create_phyloseq_obj_%J.err
#SBATCH -J 00_create_phyloseq_obj

# load modules using spack
eval $(/ref/sahlab/software/spack-0.18/bin/spack load --sh r@4.1.3%gcc@8.5.0 r-jpeg@0.1-10 \
openmpi@4.1.3 r-rmpi@0.7-2.1 r-rhtslib r-genomicalignments r-magrittr@2.0.3 r-ggpubr \
r-systemfonts r-openxlsx r-vegan r-lme4 r-gridextra r-cowplot)

# set up environment variables
export PROJ_CELIAC_PATH="/scratch/sahlab/dabid/proj_celiac"

# run the script
Rscript --no-save --vanilla ${PROJ_CELIAC_PATH}/code/00_create_phyloseq_obj.R 