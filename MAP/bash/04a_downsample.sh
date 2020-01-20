#!/bin/bash -l
#SBATCH -t 3:00:00
#SBATCH -p low
#SBATCH --mem=8G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a.out

perms='/QTL_agri/MAP/R/permutations'

Rscript $perms/04a_downsample.R --vanilla "${1}" "${2}"
