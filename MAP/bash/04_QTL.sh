#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=30G
#SBATCH --array=1-2
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a.out

script_dir='/home/jmiller1/QTL_agri/MAP'

if [ "$SLURM_ARRAY_TASK_ID" == 1 ]; then
 Rscript $script_dir/R/04_final_bin_hk_QTL.R --vanilla "${1}"
fi

if [ "$SLURM_ARRAY_TASK_ID" == 2 ]; then
 Rscript $script_dir/R/04_final_bin_em_QTL.R --vanilla "${1}"
fi

if [ "$SLURM_ARRAY_TASK_ID" == 3 ]; then
 Rscript $script_dir/R/04_final_norm_imp_QTL.R --vanilla "${1}"
fi
