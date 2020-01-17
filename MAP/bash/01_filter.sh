#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a_filter.out

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/01_${1}_filter.R
