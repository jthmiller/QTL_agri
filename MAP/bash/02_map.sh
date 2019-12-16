#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p med
#SBATCH --array=1-24
#SBATCH --mem=30G

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/02_${1}_tspmap.R --vanilla "$SLURM_ARRAY_TASK_ID"
