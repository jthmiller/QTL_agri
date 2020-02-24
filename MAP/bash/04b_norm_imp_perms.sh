#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/perms/perm%x_%a.out
#SBATCH  --error=/home/jmiller1/QTL_agri/MAP/bash/slurms/perms/perm%x_%a.err


perms="$HOME/QTL_agri/MAP/R/permutations"

Rscript $perms/04b_norm_imp_perms.R --vanilla "${1}" "${2}" "${3}" "${4}" "${SLURM_ARRAY_TASK_ID}"
