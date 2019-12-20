### MASTER

#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

## sbatch 01_filter.sh 'BRP'
## sbatch 01_filter.sh 'NBH'
## sbatch 01_filter.sh 'ELR'
##sbatch 01_filter.sh 'NEW'

## sbatch 02_map.sh 'BRP'
## sbatch 02_map.sh 'NBH'
## sbatch 02_map.sh 'ELR'
 sbatch 02_map.sh 'NEW'
## sbatch 02_map_missing.sh 'ELR'

sbatch 03_fit.sh 'NBH'
sbatch 03_fit.sh 'ELR'
sbatch 03_fit.sh 'NEW'
sbatch 03_fit.sh 'BRP'
sbatch 03_fit_missing.sh 'ELR'
##
##sbatch 04_fit_short.sh 'NBH'
##sbatch 04_fit_short.sh 'BRP'
##sbatch 04_fit_short.sh 'NEW'
##sbatch 04_fit_short.sh 'ELR'
##
sbatch 04_fit_long.sh 'NBH'
sbatch 04_fit_long.sh 'BRP'
sbatch 04_fit_long.sh 'NEW'
sbatch 04_fit_long.sh 'ELR'
