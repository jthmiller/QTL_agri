#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -p high
#SBATCH --mem=10G
#SBATCH  --output=/home/jmiller1/QTL_agri/MAP/bash/slurms/%x_%a_write.out

script_dir='/home/jmiller1/QTL_agri/MAP'

Rscript $script_dir/R/03_write_map_cross.R --vanilla "${1}"
