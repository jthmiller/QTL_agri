#!/bin/R
pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

cores <- as.numeric(commandArgs(TRUE)[2])

################################################################################
## perms.1
## perms.2
## pens
load(file.path(mpath,paste0(pop,'_all_perms_bin_hk.rsave')))
################################################################################
## bin.em.2
load(file.path(mpath,paste0(pop,'_scan2_bin_hk.rsave')))
################################################################################
#sone.o <- scanone(cross,pheno.col=4, model="binary", method="em")
#sone.a <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1])
#sone.i <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#sone.io <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#cbind(summary(sone.o),summary(sone.a)$lod,summary(sone.i)$lod,summary(sone.io)$lod)
cross$pheno <- as.data.frame(cross$pheno)
################################################################################
## cross <- subset(cross, chr=c(1:4,6:24))

no_qtl_hk <- scanone(cross, pheno.col=4, method="hk", model="binary", verbose=FALSE, tol=1e-4, maxit=10000)

add.perms <- scanone(cross, pheno.col=4, method="hk", model="binary", n.perm = 1000, n.cluster=6, perm.Xsp=T)

lod <- summary(add.perms)[[1]][1]

qtl <- summary(no_qtl_hk,lod)

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
                    run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

scan_3_imp <- scanqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
                      covar=data.frame(cross$pheno$sex), formula=y~Q1+Q2+Q3,
                      incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
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
