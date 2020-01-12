### MASTER
script_dir='/home/jmiller1/QTL_agri/MAP'

#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

## sbatch 01_filter.sh 'BRP'
## sbatch 01_filter.sh 'NBH'
## sbatch 01_filter.sh 'ELR'
## sbatch 01_filter.sh 'NEW'
## srun Rscript $script_dir/R/01b_ELR_add_AHR_genotypes.R

## sbatch 02_map.sh 'BRP'
## sbatch 02_map.sh 'NBH'
## sbatch 02_map.sh 'ELR'
## sbatch 02_map.sh 'NEW'
## sbatch 02_map_missing.sh 'ELR'

sbatch 03_write_map_cross.sh 'NBH'
sbatch 03_write_map_cross.sh 'ELR'
sbatch 03_write_map_cross.sh 'BRP'
sbatch 03_write_map_cross.sh 'NEW'

sbatch 04_QTL.sh 'NBH'
sbatch 04_QTL.sh 'ELR'
sbatch 04_QTL.sh 'BRP'
sbatch 04_QTL.sh 'NEW'
