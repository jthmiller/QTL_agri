#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
cores <- as.numeric(commandArgs(TRUE)[3])

print(commandArgs(TRUE))
print(paste(pop))

library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

################################################################################
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################

################################################################################
print(paste(cores,'cores'))
erp <- 0.001
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'
cross$pheno$pheno_norm <- nqrank(cross$pheno$Pheno)
################################################################################

################################################################################
if(pop == 'NBH'){
 mar <- '2:27373969'
 g <- pull.geno(fill.geno(cross))[,mar]
 g <- cbind(as.numeric(g==1), as.numeric(g==2))
} else {
 mar <- '18:20422142'
 g <- pull.geno(fill.geno(cross))[,mar]
 g <- cbind(as.numeric(g==1), as.numeric(g==2))
}
################################################################################

################################################################################
bin.em.2 <- scantwo(cross, pheno.col=4, model="binary", method="em",
 clean.output=T, clean.nmar=50, clean.distance=50, maxit=1000, incl.markers=T,
 assumeCondIndep=T, n.cluster=cores, use="complete.obs", addcovar=g)
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
################################################################################

### summary(bin.em.2, thresholds=c(0, Inf, 5, Inf, Inf), what="int")
