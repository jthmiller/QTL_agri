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

sbatch -J "NBH_dwns" ../permutations/04a_downsample.sh 'NBH' 22
sbatch -J "ELR_dwns" ../permutations/04a_downsample.sh 'ELR' 22
sbatch -J "ELRM_dwns" ../permutations/04a_downsample.sh 'ELR.missing' 22

################################################################################
##SCANTWO PERMUTATIONS

## 04b_bin_em_perms.sh --vanilla pop perm_count cores arraynum


bashsc="$HOME/QTL_agri/MAP/bash"

#sbatch -J "NBH_P.N.I" $bashsc/04b_norm_imp_perms.sh 'NBH' 22 22
#sbatch -J "ELR_P.N.I" $bashsc/04b_norm_imp_perms.sh 'ELR' 22 22
#sbatch -J "ELRM_P.N.I" $bashsc/04b_norm_imp_perms.sh 'ELR.missing' 22 22
## test
sbatch -J "NBH_PBE" --array=2-10 -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' 1 24
sbatch -J "ELR_PBE" --array=1-10 -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR' 12 1
sbatch -J "ELRM_PBE" --array=1-10 -p low -t 13:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR.missing' 12 1

sbatch -J "NBH_PBK" $bashsc/04b_bin_hk_perms.sh 'NBH' 22 22
sbatch -J "ELR_PBK" $bashsc/04b_bin_hk_perms.sh 'ELR' 22 22
sbatch -J "ELRM_PBK" $bashsc/04b_bin_hk_perms.sh 'ELR.missing' 22 22

################################################################################
##SCANTWO
bashsc="$HOME/QTL_agri/MAP/bash"

sbatch -J "NBH_S2NI"  -p high -t 48:00:00 $bashsc/04c_norm_imp_scan2.sh 'NBH' 22
sbatch -J "ELR_S2NI"  -p high -t 48:00:00 $bashsc/04c_norm_imp_scan2.sh 'ELR' 22
sbatch -J "ELRM_S2NI" -p high -t 48:00:00 $bashsc/04c_norm_imp_scan2.sh 'ELR.missing' 22

sbatch -J "NBH_S2BE" -p high -t 48:00:00 $bashsc/04c_bin_em_scan2.sh 'NBH' 22
sbatch -J "ELR_S2BE"  -p high -t 48:00:00 $bashsc/04c_bin_em_scan2.sh 'ELR' 22
sbatch -J "ELRM_S2BE" -p high -t 48:00:00 $bashsc/04c_bin_em_scan2.sh 'ELR.missing' 22

sbatch -J "NBH_S2BH" -p med -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'NBH' 22
sbatch -J "ELR_S2BH" -p med -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'ELR' 22
sbatch -J "ELRM_S2BH" -p med -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'ELR.missing' 22

################################################################################

################################################################################
### STEPWISE QTL
bashsc="$HOME/QTL_agri/MAP/bash"

sbatch -J "NBH_SWNI" --depend=afterany: $bashsc/04c_norm_imp_step.sh 'NBH' 22 22
sbatch -J "ELR_SWNI" --depend=afterany: $bashsc/04c_norm_imp_step.sh 'ELR' 22 22
sbatch -J "ELRM_SWNI" --depend=afterany: $bashsc/04c_norm_imp_step.sh 'ELR.missing' 22 22

sbatch -J "NBH_SWBE" --depend=afterany: $bashsc/04c_bin_em_step.sh 'NBH' 22 22
sbatch -J "ELR_SWBE" --depend=afterany: $bashsc/04c_bin_em_step.sh 'ELR' 22 22
sbatch -J "ELRM_SWBE" --depend=afterany: $bashsc/04c_bin_em_step.sh 'ELR.missing' 22 22

sbatch -J "NBH_SWBK" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'NBH' 22 22
sbatch -J "ELR_SWBK" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'ELR' 22 22
sbatch -J "ELRM_SWBK" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'ELR.missing' 22 22

################################################################################
sbatch -J "power_calc" $bashsc/06_power.sh 'NBH'
