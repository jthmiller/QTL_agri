#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p med
#SBATCH --mem=60G

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/04_final_bin_hk_QTL.R --vanilla "${1}"
Rscript $script_dir/R/04_final_bin_imp_QTL.R --vanilla "${1}"
Rscript $script_dir/R/04_final_norm_imp_QTL.R --vanilla "${1}"
