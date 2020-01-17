#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p med
#SBATCH --array=1-24
#SBATCH --mem=20G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a_qtl.out

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/02_ELR_tspmap_with_missing.R --vanilla "$SLURM_ARRAY_TASK_ID" "${1}"
