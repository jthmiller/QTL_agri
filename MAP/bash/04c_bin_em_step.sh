#!/bin/bash -l
#SBATCH -t 3:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a.out

perms="$HOME/QTL_agri/MAP/R/permutations"

Rscript $perms/04c_bin_em_step.R --vanilla "${1}" "${2}"
