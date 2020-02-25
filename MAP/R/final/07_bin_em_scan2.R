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
bin.em.2.cov <- scantwo(cross, pheno.col=4, model="binary", method="em",
 clean.output=T, clean.nmar=10, clean.distance=50, maxit=1000, incl.markers=T,
 assumeCondIndep=T, n.cluster=cores, use="complete.obs", addcovar=g)
################################################################################

################################################################################
bin.imp.2 <- scantwo(cross, pheno.col=5, model="normal", method="em",
 clean.output=T, clean.nmar=10, clean.distance=50, maxit = 1000, incl.markers=T,
 assumeCondIndep=T, n.cluster=cores, use="complete.obs")
################################################################################


################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
################################################################################

### summary(bin.em.2, thresholds=c(0, Inf, 5, Inf, Inf), what="int")

plot_test('lod',width=2000,heigh=2000)
plot(bin.em.2,col.scheme = "redblue",contours=T, zlim = c(5,5))
dev.off()

plot_test('lod',width=2000,heigh=2000)
plot(bin.em.2,col.scheme = "redblue", contours=c(15,5), zlim = c(18,16))
dev.off()
