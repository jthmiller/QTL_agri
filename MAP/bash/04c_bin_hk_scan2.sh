#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=10G
#SBATCH --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/scan2%x_%a.out
#SBATCH --error=/home/jmiller1/QTL_agri/MAP/bash/slurms/error/perm%x_%a.error

perms="$HOME/QTL_agri/MAP/R/permutations"

Rscript $perms/04c_bin_hk_scan2.R --vanilla "${1}" "${2}" "${3}"
