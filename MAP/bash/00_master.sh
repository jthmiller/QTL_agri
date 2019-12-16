### MASTER

#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

sbatch 01_filter.sh 'BRP'
sbatch 01_filter.sh 'NBH'
sbatch 01_filter.sh 'ELR'
sbatch 01_filter.sh 'NEW'

##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

sbatch 02_map.sh 'BRP'
sbatch 02_map.sh 'NBH'
sbatch 02_map.sh 'ELR'
sbatch 02_map.sh 'NEW'


sbatch elr_fitQTL.sh
sbatch brp_fitQTL.sh
sbatch new_fitQTL.sh
sbatch nbh_fitQTL.sh


sbatch fitQTL_long.sh NBH
sbatch fitQTL_long.sh BRP
sbatch fitQTL_long.sh NEW
sbatch fitQTL_long.sh ELR
