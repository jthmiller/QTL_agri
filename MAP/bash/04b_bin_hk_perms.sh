#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a.out

perms="$HOME/QTL_agri/MAP/R/permutations"

Rscript $perms/04b_bin_hk_perms.R --vanilla "${1}" "${2}" "${3}"
