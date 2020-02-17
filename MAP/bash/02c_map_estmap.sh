#!/bin/bash -l
#SBATCH -t 13:00:00
#SBATCH -p low
#SBATCH --array=1-24
#SBATCH --mem=20G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a_est_map.out

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/mapping/03_estmap.R --vanilla "$SLURM_ARRAY_TASK_ID" "${1}"
