#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=60G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/perm%x_%a.out

perms="$HOME/QTL_agri/MAP/R/permutations"

## 04b_bin_em_perms.sh --vanilla pop perm_count cores arraynum
# sbatch -J "NBH_PBE" -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' 12 1

Rscript $perms/04b_bin_em_perms.R "${1}" "${2}" "${3}" "${SLURM_ARRAY_TASK_ID}"

"--vanilla" 'NBH' '1' '12'
