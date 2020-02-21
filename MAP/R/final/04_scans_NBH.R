#!/bin/R
pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- as.numeric(commandArgs(TRUE)[2])

################################################################################
##load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################
fl <- paste0(pop,'_imp.mapped.tsp.csv')
fl <- file.path(mpath,fl)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross$pheno <- as.data.frame(cross$pheno)
################################################################################
## cross <- subset(cross, chr=c(1:4,6:24))

sone <- scanone(cross,pheno.col=4, model="binary", method="mr")
sone.perms <- scanone(cross,pheno.col=4, model="binary", method="mr", n.perm=10000, n.cluster=cores)
summary(sone, alpha=0.1, lodcolumn=1, pvalues=T, perms=sone.perms, ci.function="bayesint")
lod <- summary(sone.perms)[[2]]
qtl <- summary(sone,lod)

add.qtl1 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")

add.qtl1 <- refineqtl(cross, pheno.col = 4, qtl=add.qtl1, method = "hk", model='binary',
                      incl.markers=F)

int.em <- addint(cross, pheno.col = 4, qtl = add.qtl1, method='hk', model='binary',
                 covar=data.frame(cross$pheno$sex) ,formula=y~Q1+Q2+Q3, maxit=10000)

add_Q4_hk <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=10000, incl.markers=F,
                    formula = y~Q1+Q2+Q3+Q4)

add_Q4_wInts <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
                       incl.markers=T, verbose=FALSE, tol=1e-4, maxit=10000, incl.markers=F,
                       formula = y~Q1*Q3+Q2+Q4)



cross.grid.fixed <- calc.genoprob(cross, step=1, off.end=1, error.prob=0.01, map.function="kosambi", stepwidth="fixed")
cross.grid.fixed <- reduce2grid(cross.grid.fixed)
sone.grid <- scanone(cross.grid.fixed, pheno.col=4, model="binary", method="em")
sone.grid.perms <- scanone(cross.grid.fixed, pheno.col=4, model="binary", method="em", n.perm=1000, n.cluster=1)
lod <- summary(sone.grid.perms)[[2]]
qtl <- summary(sone.grid,lod)

add.qtl1 <- makeqtl(cross.grid.fixed, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")

add.qtl1 <- refineqtl(cross.grid.fixed, pheno.col = 4, qtl=add.qtl1, method = "hk", model='binary',
                      incl.markers=F)

int.em <- addint(cross.grid.fixed, pheno.col = 4, qtl = add.qtl1, method='hk', model='binary',
                 covar=data.frame(cross.grid.fixed$pheno$sex) ,formula=y~Q1+Q2+Q3, maxit=10000)

add_Q4_hk <- addqtl(cross.grid.fixed, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=10000, incl.markers=F,
                    formula = y~Q1+Q2+Q3+Q4)

add_Q4_wInts <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
                       incl.markers=T, verbose=FALSE, tol=1e-4, maxit=10000, incl.markers=F,
                       formula = y~Q1*Q3+Q2+Q4)
#
#int.em <- addint(cross, qtl = add.qtl1, formula=y~Q1+Q2+Q3+Q4, method='hk', model='binary')

################################################################################

fit_hk_3int <- fitqtl(cross, pheno.col=4, method="hk", model="binary", qtl = add.qtl1,
                covar=NULL, formula=y~Q1+Q2+Q3+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE, run.checks=F)

################################################################################


################################################################################
AHR.qtl <- makeqtl(cross, chr=1, pos=0, what="prob")

add_Q2_AHR <- addqtl(cross, pheno.col=4, qtl = AHR.qtl, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=10000,
            formula = y~Q1+Q2)

int_Q2_AHR <- addqtl(cross, pheno.col=4, qtl = AHR.qtl, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=10000,
            formula = y~Q1*Q2)

plot_test('only_AHR_qtl.png',width=1000)
 plot(no_qtl_hk, add_Q2_AHR, int_Q2_AHR)
dev.off()
################################################################################

#### EM ########################################################################
no_qtl_em <- scanone(cross, pheno.col=4, method="em", model="binary", maxit=10000)

add.perms <- scanone(cross, pheno.col=4, method="em", model="binary", n.perm = 10000, n.cluster=6, perm.Xsp=T)

lod <- summary(add.perms)[1]

qtl <- summary(no_qtl_em,lod)

Q3 <- makeqtl(cross, chr=qtl[['chr']][c(1,3)], pos=qtl[['pos']][c(1,3)], what="draws")

################################################################################

#### IMP #######################################################################
cross <- sim.geno(cross, stepwidth="fixed", step=1,  error.prob=erp, off.end=1, map.function="kosambi", n.draws=1000)

fit_3_em <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
  covar=NULL, formula=y~Q1+Q2, dropone=TRUE, get.ests=T,
  run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

fit_3_em_int <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
  covar=NULL, formula=y~Q1+Q2+Q3+Q1:Q3, dropone=TRUE, get.ests=T,
  run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
################################################################################


##### IMP ######################################################################
##### To fit imp, must use em to scan for first QTL  (not offered rQTL)
no_qtl_im <- scanone(cross, pheno.col=4, method="imp", model="binary")

imp_perms <- scanone(cross, pheno.col=4, method="imp", model="binary", n.perm = 10000, n.cluster=6)

lod <- summary(imp_perms)[1]

qtl <- summary(no_qtl_im, lod)

Q3 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")

Q3 <- refineqtl(cross, pheno.col = 4, qtl=Q3, method = "imp", model='binary',
                incl.markers=F)

int.imp <- addint(cross, pheno.col = 4, qtl = Q3, method='imp', model='binary',
                  covar=data.frame(cross$pheno$sex) ,formula=y~Q1+Q2+Q3, maxit=10000)

fit_3_imp <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
                    formula=y~Q1+Q2+Q3, dropone=TRUE, get.ests=T, covar=data.frame(cross$pheno$sex),
                    run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

scan_3_imp <- scanqtl(cross, pheno.col=4,method="imp", model="binary",
                      chr=as.numeric(qtl[['chr']]), pos=as.numeric(qtl[['pos']]),
                      covar=data.frame(cross$pheno$sex), formula=y~Q1*Q3+Q2,
                      incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=10000,
                      forceXcovar=FALSE)
################################################################################

##### MR #######################################################################
no_qtl_mr <- scanone(cross, pheno.col=4, method="mr", model="binary")
no_qtl_mram <- scanone(cross, pheno.col=5, method="mr-argmax", model="normal")
no_qtl_mrim <- scanone(cross, pheno.col=5, method="mr-imp")
cbind(summary(no_qtl_mr),am=summary(no_qtl_mram)$lod,im=summary(no_qtl_mrim)$lod)

mr_perms <- scanone(cross, pheno.col=4, method="mr", model="binary", n.perm = 10000, n.cluster=6)

lod <- summary(mr_perms)[1]

qtl <- summary(no_qtl_mr, lod)

Q3mr <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")

################################################################################

int.imp <- addint(cross, pheno.col = 4, qtl = Q3, method='imp', model='binary',
                  covar=data.frame(cross$pheno$sex) ,formula=y~Q1+Q2+Q3, maxit=10000)

fit_3_imp <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
                    formula=y~Q1+Q2+Q3, dropone=TRUE, get.ests=T, covar=data.frame(cross$pheno$sex),
                    run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

scan_3_imp <- scanqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
                      covar=data.frame(cross$pheno$sex), formula=y~Q1+Q2+Q3,
                      incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
                      forceXcovar=FALSE)

scan_3_imp_wo_sex <- scanqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
                      formula=y~Q1+Q2+Q3,
                      incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
                      forceXcovar=FALSE)

plot_test('qtl2_18')
effectplot(cross, pheno.col=4, mname1="2@89.0", mname2="18@47.1", var.flag="pooled")
dev.off()

plot_test('qtl2_18')
effectplot(cross, pheno.col=4, mname2="2@89.0", mname1="18@47.1", var.flag="group")
dev.off()


plot_test('qtl2_18')
effectplot(cross, pheno.col=1, mname2="2@89.0", mname1="18@47.1", var.flag="pooled")
dev.off()

geno.crosstab(cross, mname2=find.marker(cross,2,89.0), mname1=find.marker(cross,18,47.1), eliminate.zeros=F)

plot_test('qtl2_18')
effectplot(cross, pheno.col=1, mname2="2@89.0", mname1="18@47.1", var.flag="group", ylim=c(0,5))
dev.off()

find.psuedo
2 102.67


plot_test('qtl2_18')
effectplot(cross, pheno.col=1, mname2="c2.loc102", mname1="18@47.1", var.flag="group")
dev.off()

plot_test('qtl2_18')
effectplot(cross, pheno.col=4, mname2="c2.loc102", var.flag="group")
dev.off()

plot_test('qtl2_18')
effectplot(cross, pheno.col=1, mname2="2@89.0", var.flag="group")
dev.off()

geno.crosstab(cross, mname2="2@89.0", mname1="18@47.1", eliminate.zeros=F)

save.image('~/saveR.rsave')
