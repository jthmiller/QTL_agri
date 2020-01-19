#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a_perms.out

echo "${SLURM_NPROCS}"

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/04_final_norm_imp_perms.R  --vanilla "${1}" "${2}"
