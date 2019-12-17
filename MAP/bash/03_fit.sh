#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G


script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/04_${1}_fitQTL.R --vanilla "${1}"
