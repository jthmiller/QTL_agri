#!/bin/R

pop <- 'ELR'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
library(scales)


## mapped cross
mapfile <- paste0(pop,'_reorder_imp_mapped')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)


## unmapped w/additional genotypes in AHR2a and AIP
fl <- file.path(mpath,'ELR_unmapped_added_markers.csv')

cross2 <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

### Test all offspring genotyped at the deletion
cross_whoi <- pull.markers(cross2, c("AHR2a_del","AIP_252","AIP_261") )
cross_whoi$pheno$bin <- ifelse(cross_whoi$pheno$Pheno > 2, 1 , 0)
sone_whoi <- scanone(cross_whoi, pheno.col=4, method="mr", model="bin")


cross2 <- subset(cross2, ind = cross2$pheno$ID %in% cross$pheno$ID)

AHR2a_del <- pull.geno(cross2)[,"AHR2a_del"]
AIP_252 <- pull.geno(cross2)[,"AIP_252"]
AIP_261 <- pull.geno(cross2)[,"AIP_261"]

cross <- addmarker(cross, AHR2a_del, "AHR2a_del", "1", 1)
cross <- addmarker(cross, AIP_252, "AIP_252", "2", 1)
cross <- addmarker(cross, AIP_261, "AIP_261", "2", 2)

cross <- switchAlleles(cross, "AIP_252")
cross <- switchAlleles(cross, "AIP_261")

## map with deletion
cross  <- est.rf(cross)
cross  <- tspOrder(cross = cross , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
newmap <- est.map(cross, error.prob = erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross <- replace.map(cross, newmap)

i <- 100 ; plotit(cross, 'cross_noimp_exact_nopar')




erprob <- 0.001 #determined from mapping step
cross <- sim.geno(cross, n.draws = 150, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
cross <- calc.genoprob(cross, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
################################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)

sone_bin <- scanone(cross, pheno.col=4, method="hk", model="bin")
sone_reg <- scanone(cross, pheno.col=4, method="mr", model="bin")


sone_bin["AIP_252",]
### "AIP_252" LOD = 0.833792

sone_bin["AIP_261",]
### "AIP_261" LOD = 0.9058545

sone_bin["AHR2a_del",]
### "AHR2a_del" LOD = 0.1616044

### 
table(pull.geno(cross, 1)[which(cross$pheno$bin == 1),"AHR2a_del" ])