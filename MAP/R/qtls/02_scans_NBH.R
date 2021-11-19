#!/bin/R
pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- 16
################################################################################

#################################################################################
mapfile <- paste0(pop,'_reorder_imp_nopar')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#################################################################################

load(file.path(mpath,paste0(pop,'_imputed.rsave')))

################################################################################
cross$pheno <- as.data.frame(cross$pheno)
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'
cross$pheno$pheno_norm <- nqrank(cross$pheno$Pheno)
################################################################################

#save.image(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))

################################################################################
## INITIAL SCAN BINARY
################################################################################
erprob <- 0.001 #determined from mapping step
cross <- sim.geno(cross, n.draws = 150, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
cross <- calc.genoprob(cross, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
################################################################################
sone_bin <- scanone(cross, pheno.col=4, method="imp", model="bin", n.cluster = cores)

sone_bin <- scanone(cross, pheno.col=4, method="hk", model="bin", n.cluster = cores)
perm_cross <- subset(cross, chr = c(1,3:4,6:17,19:24))
## Drop chromosome 2,X, and 18 for permutations
sone.perms <- scanone(perm_cross, pheno.col=4, method="hk",  model="bin", n.perm=1000, n.cluster = cores)
lod <- summary(sone.perms)[[2]]
qtl.add <- summary(sone_bin,lod)
## make QTL object for the model
hk.qtl.4_perms <- makeqtl(cross, chr = qtl.add[['chr']], pos = qtl.add[['pos']], what="prob")
## refine the position of QTL considering all 4 QTL
hk.qtl.4_perms <- refineqtl(cross, pheno.col = 4, qtl = hk.qtl.4_perms, method = "hk", model='binary', incl.markers=T)

### without chr24 (drop one analysis indicates that chr24 QTL does not provide additional explanatory power when condidered in the full-QTL model)
lod <- summary(sone.perms)[[1]]
qtl.add <- summary(sone_bin,lod)
## make QTL object for the model
hk.qtl.3_perms <- makeqtl(cross, chr = qtl.add[['chr']], pos = qtl.add[['pos']], what="prob")
## refine the position of QTL considering all 4 QTL
hk.qtl.3_perms <- refineqtl(cross, pheno.col = 4, qtl = hk.qtl.3_perms, method = "hk", model='binary', incl.markers=T)


### 2 qtl refine
hk.qtl.2_perms <- makeqtl(cross, chr = qtl.add[['chr']][c(1,3)], pos = qtl.add[['pos']][c(1,3)], what="prob")
hk.qtl.2_perms <- refineqtl(cross, pheno.col = 4, qtl = hk.qtl.2_perms, method = "hk", model='binary', incl.markers=T)
fit2_hk_bin_perm <- fitqtl(cross, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.2_perms,
                covar=NULL, formula = y~Q1+Q2, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

ahr_genes_perm <- get_AHR(cross)
################################################################################

################################################################################
erprob_noperm <- 0.05 #determined from mapping step
noperms <- sim.geno(noperms, n.draws = 150, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
noperms <- calc.genoprob(noperms, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
## SCAN
sone_bin_noperm <- scanone(noperms, pheno.col=4, method="hk", model="bin", n.cluster = cores)
perm_noperms <- subset(noperms, chr = c(1,3:4,6:17,19:24))
## Drop chromosome 2,X, and 18 for permutations
sone.perms_noperm <- scanone(perm_noperms, pheno.col=4, method="hk",  model="bin", n.perm=1000, n.cluster = cores)
lod_noperm <- summary(sone.perms_noperm)[[2]]
qtl.add_noperm <- summary(sone_bin_noperm, lod_noperm)
## make QTL object for the model
hk.qtl.4_noperms <- makeqtl(noperms, chr = qtl.add_noperm[['chr']], pos = qtl.add_noperm[['pos']], what="prob")
## refine the position of QTL considering all 4 QTL
hk.qtl.4_noperms <- refineqtl(noperms, pheno.col = 4, qtl = hk.qtl.4_noperms, method = "hk", model='binary', incl.markers=T)





ahr_genes_noperm <- get_AHR(cross)
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################
fit3_hk_bin_perm <- fitqtl(cross, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.3_perms,
                covar=NULL, formula = y~Q1+Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

fit2_hk_bin_perm <- fitqtl(cross, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.3_perms,
                covar=NULL, formula = y~Q1+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)


fit_hk_bin_perm <- fitqtl(cross, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.4_perms,
                covar=NULL, formula = y~Q1+Q2+Q3+Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

fit_hk_bin_noperm <- fitqtl(noperms, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.4_noperms,
                covar=NULL, formula = y~Q1+Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

################################################################################


 if (any(ahr$chr==ch)) {
  ind <- which(ahr$chr==ch)
  text(ahr[ind,'pos'], 2 , labels = ahr[ind,'gene'])
 }



################# SCAN 2 #######################################################
################################################################################

################################################################################
bin.em.2.cov <- scantwo(cross, pheno.col = 4, model="binary", method="em",
 clean.output = T, clean.nmar = 10, clean.distance = 10 , maxit=1000, incl.markers=T,
 assumeCondIndep = T, n.cluster=cores)
################################################################################

################################################################################
bin.em.2.noperms <- scantwo(noperms, pheno.col = 4, model="binary", method="em",
 clean.output = T, clean.nmar = 10, clean.distance = 10 , maxit=1000, incl.markers=T,
 assumeCondIndep = T, n.cluster=cores)
################################################################################

save.image(file.path(mpath,paste0(pop,'_scan2.rsave')))

################# SCAN 2 #######################################################
################################################################################









################################################################################
################ SCAN 2 SUMMARY ################################################

## lod.full lod.fv1 lod.int lod.add lod.av1
thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6)

best <- summary(bin.em.2.cov, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="best")
best[order(best$lod.int),]
best[order(best$lod.full),]

full <- summary(bin.em.2.cov, thresholds = thresholds, what="full")
full[order(full$lod.full),]
full[order(full$lod.int),]

## Pull out positions of 2 locus interactions that support interactive loci
ints <- summary(bin.em.2.cov, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6), what="int")
ints[order(ints$lod.int),]








################################################################################
############### CONFIDENCE INTERVALS ###########################################

s1boot_q2_bin <- scanoneboot(cross, chr = 2, pheno.col = 4, model= 'binary', method = "hk")
summary(s1boot_q2_bin)
s1boot_q2_np <- scanoneboot(cross, chr = 2, pheno.col = 1, model= 'np', method = "hk")
summary(s1boot_q2_np)


s1boot_q8 <- scanoneboot(cross, chr = 2, pheno.col = 1, model= 'np', method = "hk")
s1boot_q18 <- scanoneboot(cross, chr = 2, pheno.col = 1, model= 'np', method = "hk")
s1boot_q24 <- scanoneboot(cross, chr = 2, pheno.col = 1, model= 'np', method = "hk")


plot_test('hk_perm_qtl2')
plot(s1boot_q2_bin)
dev.off()












################################################################################
## CONFID INTERVALS
bayesint(sone_bin, chr = 2, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin, chr = 8, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin, chr = 18, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin, chr = 24, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)

lodint(hk.qtl.4_perms, qtl.index = 1, expandtomarkers=TRUE)
lodint(hk.qtl.4_perms, qtl.index = 2, expandtomarkers=TRUE)
lodint(hk.qtl.4_perms, qtl.index = 3, expandtomarkers=TRUE)
lodint(hk.qtl.4_perms, qtl.index = 4, expandtomarkers=TRUE)

lodint(hk.qtl.4_perms, qtl.index = 1, expandtomarkers=TRUE)


################################################################################

bayesint(sone_bin_noperm, chr = 2, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin_noperm, chr = 8, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
bayesint(sone_bin_noperm, chr = 18, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)

lodint(hk.qtl.4_noperms, qtl.index = 1, expandtomarkers=TRUE)
lodint(hk.qtl.4_noperms, qtl.index = 2, expandtomarkers=TRUE)
lodint(hk.qtl.4_noperms, qtl.index = 3, expandtomarkers=TRUE)
################################################################################




################################################################################
############ PLOTS #############################################################

################################################################################
plot_test('hk_noperm')
plot(sone_bin_noperm, chr = c(2,8,18,24))
dev.off()

plot_test('hk_perm')
plot(sone_bin, chr = c(2,8,18,24))
dev.off()

################################################################################

















################################################################################
## COVARIATES

if(pop == 'NBH'){
 mar <- "2:37159317"
 g.add <- pull.geno(fill.geno(cross))[,mar]
 g.add <- data.frame(cbind(as.numeric(g.add == 1), as.numeric(g.add == 2)))

 mar <- find.marker(cross,13,33.9)
 g.int <- pull.geno(fill.geno(cross))[,mar]
 g.int <- data.frame(cbind(as.numeric(g.int == 1), as.numeric(g.int == 2)))

} else {
 mar <- '18:20422142'
 g <- pull.geno(fill.geno(cross))[,mar]
 g <- cbind(as.numeric(g==1), as.numeric(g==2))
}

 mar <- find.marker(cross,qtl$chr[2],qtl$pos[2])
 g.add <- pull.geno(fill.geno(cross))[,mar]
 g.add <- data.frame(cbind(as.numeric(g.add == 1), as.numeric(g.add == 2)))

 mar <- find.marker(cross,qtl$chr[1],qtl$pos[1])
 g.int <- pull.geno(fill.geno(cross))[,mar]
 g.int <- data.frame(cbind(as.numeric(g.int == 1), as.numeric(g.int == 2)))


sone.int1 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.add)
sone.int2 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.int)
sone.add1 <- scanone(cross, pheno.col=5, model="normal", method="hk", addcovar=g.add)
sone.add2 <- scanone(cross, pheno.col=5, model="normal", method="hk", addcovar=g.int)
sone.1 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.add, addcovar=g.int)
sone.2 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.int, addcovar=g.add)

sone.int1 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.add)
int.cov <- summary(sone.int1)[13,]

sone.int2 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.int)

mara <- find.marker(cross,2,85)
g.add <- pull.geno(fill.geno(cross))[,mara]
g.add <- data.frame(cbind(as.numeric(g.add == 1), as.numeric(g.add == 2)))

marb <- find.marker(cross,13,39.4)
g.int <- pull.geno(fill.geno(cross))[,marb]
g.int <- data.frame(cbind(as.numeric(g.int == 1), as.numeric(g.int == 2)))

sone.int.add <- scanone(cross, pheno.col=5, model="normal", method="hk", addcovar=g.add, intcovar=g.int)
summary(sone.int.add)

sone.add <- scanone(cross, pheno.col=5, model="normal", method="hk", addcovar=g.add)
summary(sone.add)

sone.add.perms <- scanone(cross, pheno.col=5, model="normal", method="hk", n.perm=1000, n.cluster=cores, addcovar=g.add)
lod <- summary(sone.add.perms)[[2]]
qtl.add <- summary(sone.add,lod)

sone <- scanone(cross, pheno.col=5, model="normal", method="hk")
sone.perms <- scanone(cross, pheno.col=5, model="normal", method="hk", n.perm=1000, n.cluster=cores)
lod <- summary(sone.perms)[[2]]
qtl <- summary(sone,lod)

qtl5 <- rbind(qtl.int,qtl)

qtl5 <- qtl
hk.qtl.5 <- makeqtl(cross, chr=qtl5[['chr']], pos=qtl5[['pos']], what="prob")
hk.qtl.5 <-  addtoqtl(cross, qtl = hk.qtl.5, chr = 13, pos =  30.5)


################################################################################
#### FIT
hk.qtl.5.n <- refineqtl(cross, pheno.col = 5, qtl=hk.qtl.5, method = "hk", model='normal',incl.markers=T)
fit_hk_2wINT_normal <- fitqtl(cross, pheno.col=5, method="hk", model="normal", qtl = hk.qtl.5.n,
                covar=NULL, formula = y~Q1+Q2+Q3+Q4+Q1:Q6, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_hk_2wINT_normal)
hk.qtl.5.b <- refineqtl(cross, pheno.col = 4, qtl=hk.qtl.5, method = "hk", model='binary',incl.markers=T)
fit_hk_2wINT_bin <- fitqtl(cross, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.5.b,
                covar=NULL, formula = y~Q1+Q2+Q3+Q4+Q5+Q1:Q6, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_hk_2wINT_bin)

hk.qtl.5.o <- refineqtl(cross, pheno.col = 1, qtl=hk.qtl.5, method = "hk", model='normal',incl.markers=T)
fit_hk_2wINT <- fitqtl(cross, pheno.col=1, method="hk", model="normal", qtl = hk.qtl.5.o,
                covar=NULL, formula = y~Q2+Q3+Q4+Q2:Q6, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_hk_2wINT)

imp.qtl.5 <- refineqtl(cross, pheno.col = 1, qtl=hk.qtl.5, method = "imp", model='normal',incl.markers=T)
fit_imp_2wINT <- fitqtl(cross, pheno.col=1, method="imp", model="normal", qtl = imp.qtl.5,
                covar=NULL, formula = y~Q2+Q3+Q4+Q2:Q6, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_2wINT)

fit_imp_2wINT.all <- fitqtl(cross, pheno.col=1, method="imp", model="normal", qtl = imp.qtl.5,
                covar=NULL, formula = y~Q2+Q3+Q4+Q5+Q2:Q6, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_2wINT.all)
################################################################################

load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))


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
