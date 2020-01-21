### MASTER
script_dir='/home/jmiller1/QTL_agri/MAP'

#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

## sbatch 01_filter.sh 'BRP'
## sbatch -J "NBH" filter/01_filter.sh 'NBH'
## sbatch -J "ELR" filter/01_filter.sh 'ELR'
## sbatch -J "NEW" filter/01_filter.sh 'NEW'
## srun Rscript $script_dir/R/01b_ELR_add_AHR_genotypes.R

## sbatch 02_map.sh 'BRP'
## sbatch -J "NBH_map" mapping/02_map.sh 'NBH'
## sbatch -J "ELR_map" mapping/02_map.sh 'ELR'
## sbatch -J "NEW" mapping/02_map.sh 'NEW'
## sbatch -J "ELR_M" mapping/02_map_missing.sh 'ELR'

##sbatch -J "NBH" --depend=afterany:17464611 mapping/03_write_map_cross.sh 'NBH'
sbatch -J "NBH_wc" mapping/03_write_map_cross.sh 'NBH'
sbatch -J "ELR_wc" mapping/03_write_map_cross.sh 'ELR'
##sbatch -J "BRP" mapping/03_write_map_cross.sh 'BRP'
##sbatch -J "NEW"  mapping/03_write_map_cross.sh 'NEW'

sbatch -J "NBH_qtl" ../models/04_QTL.sh 'NBH'
sbatch -J "ELR_qtl" ../models/04_QTL.sh 'ELR'
sbatch -J "ELR_Mis_qtl" ../models/04_QTL.sh 'ELR.missing'
##sbatch -J "BRP" models/04_QTL.sh 'BRP'
##sbatch -J "NEW" models/04_QTL.sh 'NEW'
##

sbatch -J "NBH_perms" ../permutations/05_perms.sh 'NBH' 22
sbatch -J "ELR_perms" ../permutations/05_perms.sh 'ELR' 22
sbatch -J 'ELR_perms' ../permutations/05_perms.sh 'ELR.missing' 22

################################################################################

sbatch -J "NBH_perms" ../permutations/04a_downsample.sh 'NBH' 22
sbatch -J "ELR_perms" ../permutations/04a_downsample.sh 'ELR' 22
sbatch -J "ELRM_perms" ../permutations/04a_downsample.sh 'ELR.missing' 22

################################################################################
bashsc="$HOME/QTL_agri/MAP/bash"

#sbatch -J "NBH_P.N.I" $bashsc/04b_norm_imp_perms.sh 'NBH' 22 22
#sbatch -J "ELR_P.N.I" $bashsc/04b_norm_imp_perms.sh 'ELR' 22 22
#sbatch -J "ELRM_P.N.I" $bashsc/04b_norm_imp_perms.sh 'ELR.missing' 22 22

sbatch -J "NBH_P.B.E" $bashsc/04b_bin_em_perms.sh 'NBH' 12 22
sbatch -J "ELR_P.B.E" $bashsc/04b_bin_em_perms.sh 'ELR' 12 22
sbatch -J "ELRM_P.B.E" $bashsc/04b_bin_em_perms.sh 'ELR.missing' 12 22

sbatch -J "NBH_P.B.K" $bashsc/04b_bin_hk_perms.sh 'NBH' 22 22
sbatch -J "ELR_P.B.K" $bashsc/04b_bin_hk_perms.sh 'ELR' 22 22
sbatch -J "ELRM_P.B.K" $bashsc/04b_bin_hk_perms.sh 'ELR.missing' 22 22

################################################################################
bashsc="$HOME/QTL_agri/MAP/bash"

sbatch -J "NBH_N.I"  $bashsc/04c_norm_imp_step.sh 'NBH' 22 22
sbatch -J "ELR_N.I"  $bashsc/04c_norm_imp_step.sh 'ELR' 22 22
sbatch -J "ELRM_N.I" $bashsc/04c_norm_imp_step.sh 'ELR.missing' 22 22

sbatch -J "NBH_B.E" --depend=afterany:17491470 $bashsc/04c_bin_em_step.sh 'NBH' 22 22
sbatch -J "ELR_B.E" --depend=afterany:17491471 $bashsc/04c_bin_em_step.sh 'ELR' 22 22
sbatch -J "ELRM_B.E" --depend=afterany:17491472 $bashsc/04c_bin_em_step.sh 'ELR.missing' 22 22

sbatch -J "NBH_B.K" --depend=afterany:17491473 $bashsc/04c_bin_hk_step.sh 'NBH' 22 22
sbatch -J "ELR_B.K" --depend=afterany:17491474 $bashsc/04c_bin_hk_step.sh 'ELR' 22 22
sbatch -J "ELRM_B.K" --depend=afterany:17491475 $bashsc/04c_bin_hk_step.sh 'ELR.missing' 22 22

################################################################################
