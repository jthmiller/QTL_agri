#!/bin/R
pop <- 'ELR'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- 16

#################################################################################
mapfile <- paste0(pop,'_reorder_imp_mapped')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#################################################################################

#load(file.path(mpath,paste0(pop,'_imputed.rsave')))

################################################################################
cross$pheno <- as.data.frame(cross$pheno)
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'
cross$pheno$pheno_norm <- nqrank(cross$pheno$Pheno)
################################################################################

#save.image(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))

################################################################################
## INITIAL SCAN
################################################################################
erprob <- 0.01 #determined from mapping step
cross <- sim.geno(cross, n.draws = 150, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
cross <- calc.genoprob(cross, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
################################################################################
sone_bin <- scanone(cross, pheno.col=4, method="hk", model="bin", n.cluster = cores)
perm_cross <- subset(cross, chr = c(1,2:4,6:12,14,16:24))
## Drop chromosome 2,X, and 18 for permutations
sone.perms <- scanone(perm_cross, pheno.col=4, method="hk",  model="bin", n.perm=1000, n.cluster = cores)

lod <- summary(sone.perms)[[2]]
qtl.add <- summary(sone_bin,lod)
## make QTL object for the model
hk.qtl.1_perms <- makeqtl(cross, chr = qtl.add[['chr']], pos = qtl.add[['pos']], what="prob")
## refine the position of QTL considering all 4 QTL
hk.qtl.1_perms <- refineqtl(cross, pheno.col = 4, qtl = hk.qtl.1_perms, method = "hk", model='binary', incl.markers=T)
################################################################################

################################################################################
erprob_noperm <- 0.05 #determined from mapping step
noperms <- sim.geno(noperms, n.draws = 150, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
noperms <- calc.genoprob(noperms, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
## SCAN
sone_bin_noperm <- scanone(noperms, pheno.col=4, method="hk", model="bin", n.cluster = cores)
perm_cross <- subset(cross, chr = c(1,2:4,6:12,14:24))
## Drop chromosome 2,X, and 18 for permutations
sone.perms_noperm <- scanone(perm_cross, pheno.col=4, method="hk",  model="bin", n.perm=1000, n.cluster = cores)
lod_noperm <- summary(sone.perms_noperm)[[2]]
qtl.add_noperm <- summary(sone_bin_noperm, lod_noperm)
## make QTL object for the model
hk.qtl.4_noperms <- makeqtl(noperms, chr = qtl.add_noperm[['chr']], pos = qtl.add_noperm[['pos']], what="prob")
## refine the position of QTL considering all 4 QTL
hk.qtl.4_noperms <- refineqtl(noperms, pheno.col = 4, qtl = hk.qtl.4_noperms, method = "hk", model='binary', incl.markers=T)
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################

################################################################################
fit_hk_bin_perm <- fitqtl(cross, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.1_perms,
                covar=NULL, formula = y~Q1, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

fit_hk_bin_noperm <- fitqtl(noperms, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.4_noperms,
                covar=NULL, formula = y~Q1+Q2, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=100000, forceXcovar=FALSE)

################################################################################

################################################################################
## CONFID INTERVALS
bayesint(sone_bin, chr = 13, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin, chr = 15, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin, chr = 18, prob=0.99, lodcolumn=1, expandtomarkers=TRUE)

lodint(hk.qtl.4_perms, qtl.index = 1, expandtomarkers=TRUE)
################################################################################

################################################################################
## CONFID INTERVALS sone_bin_noperm
bayesint(sone_bin_noperm, chr = 2, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)

bayesint(sone_bin_noperm, chr = 13, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin_noperm, chr = 15, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin_noperm, chr = 18, prob=0.99, lodcolumn=1, expandtomarkers=TRUE)

lodint(hk.qtl.4_noperms, qtl.index = 1, expandtomarkers=TRUE)
lodint(hk.qtl.4_noperms, qtl.index = 2, expandtomarkers=TRUE)
################################################################################

plot_test('hk_noperm')
plot(sone_bin_noperm, chr = c(13,15,18))
dev.off()

plot_test('hk_perm')
plot(sone_bin, chr = c(13,15,18))
dev.off()
################################################################################

ahr_genes <- get_AHR(cross)

################################################################################
bin.em.2.cov <- scantwo(cross, pheno.col = 4, model="binary", method="em",
 clean.output = T, clean.nmar = 10, clean.distance = 10, maxit=1000, incl.markers=T,
 assumeCondIndep = T, n.cluster=cores)
################################################################################

################################################################################
bin.em.2.cov <- scantwo(noperms, pheno.col = 4, model="binary", method="em",
 clean.output = T, clean.nmar = 10, clean.distance = 10, maxit=1000, incl.markers=T,
 assumeCondIndep = T, n.cluster=cores)
################################################################################

save.image(file.path(mpath,paste0(pop,'_scan1.rsave')))














(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))


plot_test('terrain', width = 1000, height = 1000)
 plot(bin.imp.2, zmax = c(10,8), col.scheme = "terrain", contours=T)
dev.off()

plot_test('terrain', width = 1000, height = 1000)
 plot(bin.imp.2, zlim = c(12,6), col.scheme = "terrain", contours=c(2,2))
dev.off()

plot_test('terrain', width = 1000, height = 1000)
 plot(norm.em.2, zlim = c(12,6), col.scheme = "terrain", contours=c(2,2))
dev.off()


save.image(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
