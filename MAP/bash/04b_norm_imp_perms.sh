#!/bin/bash -l
#SBATCH -t 3:00:00
#SBATCH -p high
#SBATCH --mem=8G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a.out

perms="$HOME/QTL_agri/MAP/R/permutations"

Rscript $perms/04b_norm_imp_perms.R --vanilla "${1}" "${2}"
