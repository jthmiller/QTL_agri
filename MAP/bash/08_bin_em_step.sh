#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/em_step_%x_%a.out
#SBATCH  --error=/home/jmiller1/QTL_agri/MAP/bash/slurms/em_step_%x_%a.err

perms="$HOME/QTL_agri/MAP/R/final/step"

Rscript $perms/08_bin_em_step.R --vanilla "${1}" "${2}" "${3}"
