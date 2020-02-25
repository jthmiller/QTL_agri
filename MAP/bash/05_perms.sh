#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=30G
#SBATCH --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/perms/%x_%a.out
#SBATCH --error=/home/jmiller1/QTL_agri/MAP/bash/slurms/perms/%x_%a.error

perms="$HOME/QTL_agri/MAP/R/final"

## 04b_bin_em_perms.sh --vanilla pop perm_count cores arraynum
# sbatch -J "NBH_PBE" -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' 12 1

Rscript $perms/05_perms.R "${1}" "${2}" "${3}" "${4}" "${5}" "${6}" "${SLURM_ARRAY_TASK_ID}"

#vanilla <- as.numeric(commandArgs(TRUE)[1])
#pop <- as.numeric(commandArgs(TRUE)[2])
#perm_count <- as.numeric(commandArgs(TRUE)[3])
#cores <- as.numeric(commandArgs(TRUE)[4])
#model <- as.numeric(commandArgs(TRUE)[5])
#method <- as.numeric(commandArgs(TRUE)[6])
#arraynum <- as.numeric(commandArgs(TRUE)[7])
