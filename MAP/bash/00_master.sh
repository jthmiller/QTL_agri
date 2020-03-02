### MASTER
script_dir='/home/jmiller1/QTL_agri/MAP'

#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

### FILTER AND ORDER MARKERS
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_map" --mem=12G -p high --array=1-24 $bashsc/01_NBH_filter_impute.sh 'NBH'
sbatch -J "ELR_map" --mem=12G -p high --array=1-24 $bashsc/01_ELR_filter_impute.sh 'NBH'
#################################################################################

### ESTIMATE MAP
bashsc="$HOME/QTL_agri/MAP/bash"
srun -J "NBH_map" --mem=12G -p high $bashsc/02c_map_estmap.sh 'NBH'
sbatch -J "ELR_map" --mem=12G -p high --array=1-24 $bashsc/02c_map_estmap.sh 'ELR'
#################################################################################

### WRITE MAP
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "ELR" $bashsc/03_write_map_cross.sh 'ELR'
#################################################################################

## RUN SCAN1 SCRIPT MANUAL TO BUILD QTL MODELS

################################################################################
### SCAN2 PERMS
bashsc="$HOME/QTL_agri/MAP/bash"

args=( '--vanilla' 'NBH' '1' '1' 'binary' )
#args=( '--vanilla' 'ELR' '1' '1' 'binary' )
varem=$(sbatch --mem=10G -p low --array=1-50 -t 6:00:00 -J "PERM.EM.NBH" $bashsc/05_perms.sh "${args[@]}" 'em' | cut -f4 -d' ')
varhk=$(sbatch --mem=10G -p low --array=1-50 -t 6:00:00 -J "PERM.HK.NBH" $bashsc/05_perms.sh "${args[@]}" 'hk' | cut -f4 -d' ')
varmr=$(sbatch --mem=10G -p low --array=1-50 -t 6:00:00 -J "PERM.MR.NBH" $bashsc/05_perms.sh "${args[@]}" 'mr' | cut -f4 -d' ')

args=( '--vanilla' 'NBH' '1' '1' 'normal' 'imp' )
#args=( '--vanilla' 'ELR' '1' '1' 'normal' 'imp' )
varim=$(sbatch --mem=10G -p low --array=1-50 -t 6:00:00 -J "PERM.IMP.NBH" $bashsc/05_perms.sh "${args[@]}" | cut -f4 -d' ')

Rscript $perms/06_combine_permutations.R --vanilla 'NBH' "_all_perms_bin_em.rsave"
Rscript $perms/06_combine_permutations.R --vanilla 'NBH' "_all_perms_bin_hk.rsave"
Rscript $perms/06_combine_permutations.R --vanilla 'NBH' "_all_perms_bin_mr.rsave"
Rscript $perms/06_combine_permutations.R --vanilla 'NBH' "_all_perms_bin_imp.rsave"
################################################################################

### TWO LOCUS SCAN $1=population $2=number of cores to use
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "EMB.NBH" --mem=60G -p low -t 13:00:00 $bashsc/07_bin_em_scan2.sh 'NBH' 22
sbatch -J "MRB.NBH" --mem=60G -p high -t 13:00:00 $bashsc/07_bin_mr_scan2.sh 'NBH' 22
sbatch -J "EMN.NBH" --mem=60G -p high -t 13:00:00 $bashsc/07_norm_em_scan2.sh 'NBH' 22
sbatch -J "NIMP.NBH" --mem=60G -p med -t 13:00:00 $bashsc/07_norm_imp_scan2.sh 'NBH' 22
sbatch -J "MRN.NBH" --mem=60G -p high -t 13:00:00 $bashsc/07_norm_mr_scan2.sh 'NBH' 22


perms="$HOME/QTL_agri/MAP/R/final"
sbatch -J "EMB.ELR" --mem=60G -p low -t 13:00:00 $bashsc/07_bin_em_scan2.sh 'ELR' 22
sbatch -J "MRB.ELR" --mem=60G -p low -t 13:00:00 $bashsc/07_bin_mr_scan2.sh 'ELR' 22
sbatch -J "EMN.ELR" --mem=60G -p low -t 13:00:00 $bashsc/07_norm_em_scan2.sh 'ELR' 22
sbatch -J "NIMP.ELR" --mem=60G -p low -t 13:00:00 $bashsc/07_norm_imp_scan2.sh 'ELR' 22
sbatch -J "MRN.ELR" --mem=60G -p high -t 13:00:00 $bashsc/07_norm_mr_scan2.sh 'ELR' 22
#################################################################################

### BIN HK STEPWISE QTL
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_SWBH" $bashsc/04c_bin_hk_step.sh 'NBH' 22 22


### BIN HK STEPWISE QTL
bashsc="$HOME/QTL_agri/MAP/bash/"
sbatch -J "NBH_link" $bashsc/10_linkage.sh 'NBH' 22
sbatch -J "ELR_link" $bashsc/10_linkage.sh 'ELR' 22
#################################################################################

#################################################################################
#################################################################################
#################################################################################







### MASTER
script_dir='/home/jmiller1/QTL_agri/MAP'

#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH" $bashsc/01_filter.sh 'NBH' '17' '0.075' '2' '0.00001'
sbatch -J "ELR" $bashsc/01_filter.sh 'ELR' '14' '0.1' '2' '0.00001'
## sbatch 01_filter.sh 'BRP'
## sbatch -J "NEW" $bashsc/01_filter.sh 'NEW'
## srun Rscript $script_dir/R/01b_ELR_add_AHR_genotypes.R

## ORDER MARKERS
#bashsc="$HOME/QTL_agri/MAP/bash"
#sbatch -J "NBH_map" -p low --array=1-24 $bashsc/02a_map.sh 'NBH'
#sbatch -J "ELR_map" -p low --array=1-24 $bashsc/02a_map.sh 'ELR'

## sbatch -J "BRP" $bashsc/02_map.sh 'BRP'
## sbatch -J "NEW"  $bashsc/02_map.sh 'NEW'
## sbatch -J "ELR_M"  $bashsc/02_map_missing.sh 'ELR'

### ESTIMATE MAP
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_map" --mem=12G -p high --array=1-24 $bashsc/02c_map_estmap.sh 'NBH'
sbatch -J "ELR_map" --mem=12G -p high --array=1-24 $bashsc/02c_map_estmap.sh 'ELR'

### WRITE MAP
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_wc" $bashsc/03_write_map_cross.sh 'NBH'
sbatch -J "ELR_wc" $bashsc/03_write_map_cross.sh 'ELR'
##sbatch -J "NBH" --depend=afterany:17464611 mapping/03_write_map_cross.sh 'NBH'
##sbatch -J "BRP" mapping/03_write_map_cross.sh 'BRP'
##sbatch -J "NEW"  mapping/03_write_map_cross.sh 'NEW'

################################################################################
## Downsample loci
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_dwns"  $bashsc/04a_downsample.sh 'NBH' 1
sbatch -J "ELR_dwns" $bashsc/04a_downsample.sh 'ELR' 1
sbatch -J "ELRM_dwns" $bashsc/04a_downsample.sh 'ELR.missing' 1
################################################################################

################################################################################
##SCANTWO BIN EM
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_S2BE" -p high -t 48:00:00 $bashsc/04c_bin_em_scan2.sh 'NBH' 22
sbatch -J "ELR_S2BE"  -p high -t 48:00:00 $bashsc/04c_bin_em_scan2.sh 'ELR' 22
sbatch -J "ELRM_S2BE" -p high -t 12:00:00 $bashsc/04c_bin_em_scan2.sh 'ELR.missing' 22
################################################################################

################################################################################
##SCANTWO BIN EM PERMUTATIONS
##--depend=afterok:"${var1}_80"
## 04b_bin_hk_perms.sh --vanilla pop perm_count cores
bashsc="$HOME/QTL_agri/MAP/bash"

var1=$(sbatch \
 --mem=5G -p low --array=1-200 -t 1:00:00 \
 -J "NBH_PBE" \
 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' '1' '1' \
 | cut -f4 -d' ')

var2=$(sbatch \
 --mem=5G -p low --array=1-200 -t 1:00:00 \
 -J "ELR_PBE" \
 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR' '1' '1' \
 | cut -f4 -d' ')

var3=$(sbatch \
  --mem=5G -p low --array=1-200 -t 1:00:00 \
  -J "ELRM_PBE" \
  $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR.missing' '1' '1' \
 | cut -f4 -d' ')
################################################################################

################################################################################
### RQTL does not support bin EM stepwise QTL
################################################################################

################################################################################
##SCANTWO BIN HK
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_S2BH" --mem=6G -p high -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'NBH' 2
sbatch -J "ELR_S2BH" --mem=6G -p high -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'ELR' 2
sbatch -J "ELRM_S2BH" --mem=6G -p high -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'ELR.missing' 2
################################################################################

################################################################################
##SCANTWO BIN HK PERMUTATIONS
##--depend=afterok:"${var1}_80"
## 04b_bin_hk_perms.sh --vanilla pop perm_count cores
bashsc="$HOME/QTL_agri/MAP/bash"

var1=$(sbatch \
 --mem=5G -p low --array=1-100 -t 3:00:00 \
 -J "NBH_PBH" \
 $bashsc/04b_bin_hk_perms.sh "--vanilla" 'NBH' '1' '1' \
 | cut -f4 -d' ')

var2=$(sbatch \
 --mem=5G -p low --array=1-100 -t 3:00:00 \
 -J "ELR_PBH" \
 $bashsc/04b_bin_hk_perms.sh "--vanilla" 'ELR' '1' '1' \
 | cut -f4 -d' ')

var3=$(sbatch \
  --mem=5G -p low --array=1-100 -t 3:00:00 \
  -J "ELRM_PBH" \
  $bashsc/04b_bin_hk_perms.sh "--vanilla" 'ELR.missing' '1' '1' \
 | cut -f4 -d' ')
################################################################################

################################################################################
### BIN HK STEPWISE QTL
bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "NBH_SWBH" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'NBH' 22 22
sbatch -J "ELR_SWBH" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'ELR' 22 22
sbatch -J "ELRM_SWBH" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'ELR.missing' 22 22
################################################################################


sbatch -J "NBH_scans" $bashsc/05_scans_NBH.sh 12













sbatch -J "NBH_S2NI"  -p high -t 48:00:00 $bashsc/04c_norm_imp_scan2.sh 'NBH' 22
sbatch -J "ELR_S2NI"  -p high -t 48:00:00 $bashsc/04c_norm_imp_scan2.sh 'ELR' 22
sbatch -J "ELRM_S2NI" -p high -t 48:00:00 $bashsc/04c_norm_imp_scan2.sh 'ELR.missing' 22


#sbatch -J "NBH_S2BH" -p med -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'NBH' 22
#sbatch -J "ELR_S2BH" -p med -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'ELR' 22
#sbatch -J "ELRM_S2BH" -p med -t 48:00:00 $bashsc/04c_bin_hk_scan2.sh 'ELR.missing' 22

################################################################################

################################################################################
##SCANTWO PERMUTATIONS
## 04b_bin_em_perms.sh --vanilla pop perm_count cores
## 04b_bin_em_perms.sh --vanilla pop perm_count cores arraynum
# sbatch -J "NBH_PBE" -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' 12 1
bashsc="$HOME/QTL_agri/MAP/bash"

###sbatch -J "NBH_PBE" --mem=3G -p high --array=1-1%25  -t 48:00:00 -o "%j" $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' '1' '1'
###sbatch -J "ELR_PBE" --depend=afterok:17700225_95 --mem=3G -p high --array=1-100%25 -t 48:00:00 $bashsc/04b_bin_em_perms.sh '--vanilla' 'ELR' '1' '1'
###sbatch -J "ELRM_PBE" --depend=afterok:17700280_95 --mem=3G -p high --array=1-100%25 -t 48:00:00 $bashsc/04b_bin_em_perms.sh '--vanilla' 'ELR.missing' '1' '1'

###sbatch -J "NBH_P.N.I" $bashsc/04b_norm_imp_perms.sh 'NBH' 22 22
###sbatch -J "ELR_P.N.I" $bashsc/04b_norm_imp_perms.sh 'ELR' 22 22
###sbatch -J "ELRM_P.N.I" $bashsc/04b_norm_imp_perms.sh 'ELR.missing' 22 22

var1=$(sbatch -J "NBH_PBE" --mem=3G -p high --array=1-100%25  -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' '1' '1' | cut -f4 -d' ')
var2=$(sbatch -J "ELR_PBE" --depend=afterok:"${var1}_80" --mem=3G -p high --array=1-100%25  -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR' '1' '1' | cut -f4 -d' ')
var3=$(sbatch -J "ELRM_PBE" --depend=afterok:"${var2}_80" --mem=3G -p high --array=1-100%25  -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR.missing' '1' '1' | cut -f4 -d' ')

var1=$(sbatch -J "NBH_PBE" --mem=6G -p high --array=1-1%25 -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' '1' '1' | cut -f4 -d' ')
var2=$(sbatch -J "ELR_PBE" --mem=6G -p high --array=1-100%20 -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR' '1' '1' | cut -f4 -d' ')
var3=$(sbatch -J "ELRM_PBE" --mem=6G -p high --array=1-100%25 -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR.missing' '1' '1' | cut -f4 -d' ')


################################################################################

### STEPWISE QTL
bashsc="$HOME/QTL_agri/MAP/bash"

sbatch -J "NBH_SWNI" --depend=afterany: $bashsc/04c_norm_imp_step.sh 'NBH' 22 22
sbatch -J "ELR_SWNI" --depend=afterany: $bashsc/04c_norm_imp_step.sh 'ELR' 22 22
sbatch -J "ELRM_SWNI" --depend=afterany: $bashsc/04c_norm_imp_step.sh 'ELR.missing' 22 22

sbatch -J "NBH_SWBE" --depend=afterany: $bashsc/04c_bin_em_step.sh 'NBH' 22 22
sbatch -J "ELR_SWBE" --depend=afterany: $bashsc/04c_bin_em_step.sh 'ELR' 22 22
sbatch -J "ELRM_SWBE" --depend=afterany: $bashsc/04c_bin_em_step.sh 'ELR.missing' 22 22

################################################################################

bashsc="$HOME/QTL_agri/MAP/bash"
sbatch -J "power_calc" $bashsc/06_power.sh 'NBH'

################################################################################


























################################################################################
## test
sbatch -J "NBH_PBE" --array=1-5 -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' '1' '12'
sbatch -J "ELR_PBE" --array=1-5 -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR' '1' '12'
sbatch -J "ELRM_PBE" --array=1-5 -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR.missing' '1' '12'

sbatch -J "NBH_PBK" $bashsc/04b_bin_hk_perms.sh 'NBH' 22 22
sbatch -J "ELR_PBK" $bashsc/04b_bin_hk_perms.sh 'ELR' 22 22
sbatch -J "ELRM_PBK" $bashsc/04b_bin_hk_perms.sh 'ELR.missing' 22 22

################################################################################

sbatch -J "NBH_SWBK" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'NBH' 22 22
sbatch -J "ELR_SWBK" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'ELR' 22 22
sbatch -J "ELRM_SWBK" --depend=afterany: $bashsc/04c_bin_hk_step.sh 'ELR.missing' 22 22

sbatch -J "NBH_PBE" --depend=afterok:17542819 --array=6-100 -p high -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'NBH' 1 12
sbatch -J "ELR_PBE" -p med --array=1-100 -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR' 1 12
sbatch -J "ELRM_PBE" -p med --array=1-100 -t 48:00:00 $bashsc/04b_bin_em_perms.sh "--vanilla" 'ELR.missing' 1 12
