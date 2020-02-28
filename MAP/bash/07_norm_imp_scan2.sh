#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/imp_scan2_%x_%a.out
#SBATCH  --error=/home/jmiller1/QTL_agri/MAP/bash/slurms/imp_scan2_%x_%a.err

perms="$HOME/QTL_agri/MAP/R/final/scan2"

Rscript $perms/07_norm_imp_scan2.R --vanilla "${1}" "${2}"
