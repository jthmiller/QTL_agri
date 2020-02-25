#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a_perms.out

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/07_NBH.R  --vanilla "${1}" "${2}"
