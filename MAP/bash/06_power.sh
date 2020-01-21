#!/bin/bash -l
#SBATCH -t 12:00:00
#SBATCH -p low
#SBATCH --mem=10G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a_power.out

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/models/power.R  --vanilla "${1}"
