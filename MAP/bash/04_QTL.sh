#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p med
#SBATCH --mem=60G
#SBATCH --array=1-3

script_dir='/home/jmiller1/QTL_agri/MAP'


if [ "$SLURM_ARRAY_TASK_ID" == 1 ]; then
 Rscript $script_dir/R/04_final_bin_hk_QTL.R --vanilla "${1}"
fi

##if [ "$SLURM_ARRAY_TASK_ID" == 2 ]; then
## Rscript $script_dir/R/04_final_bin_imp_QTL.R --vanilla "${1}"
##fi
##
##if [ "$SLURM_ARRAY_TASK_ID" == 3 ]; then
## Rscript $script_dir/R/04_final_norm_imp_QTL.R --vanilla "${1}"
##fi
