#!/bin/R
pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- 6
################################################################################
##fl <- "NBH_5755_imputed_NW_tsp.csv"
fl <- "NBH_7370_imputed_tsp_mapped.csv"
cross <- read.cross(file=fl , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross$pheno <- as.data.frame(cross$pheno)
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'
cross$pheno$pheno_norm <- nqrank(cross$pheno$Pheno)
################################################################################

################################################################################

error <- 0.001
cross <- sim.geno(cross,n.draws=160, error.prob=error, map.function="kosambi", stepwidth="fixed")
cross <- calc.genoprob(cross, error.prob=error, map.function="kosambi", stepwidth="fixed")
cross <- est.rf(cross, maxit=100, tol=1e-6)

save.image(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))

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
