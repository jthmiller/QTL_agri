#!/bin/R
pop <- 'ELR.missing'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
## perms.1
## perms.2
## pens
load(file.path(mpath,paste0(pop,'_all_perms_bin_em.rsave')))
################################################################################
## bin.em.2
load(file.path(mpath,paste0(pop,'_scan2_bin_hk.rsave')))
################################################################################
#sone.o <- scanone(cross,pheno.col=4, model="binary", method="em")
#sone.a <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1])
#sone.i <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#sone.io <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#cbind(summary(sone.o),summary(sone.a)$lod,summary(sone.i)$lod,summary(sone.io)$lod)

################################################################################
no_qtl <- scanone(cross, pheno.col=4, method="hk", model="binary", verbose=FALSE, tol=1e-4, maxit=1000)

add_Q3 <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q3)


add_Q3 <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q3)

add_Q3Q1_int <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q1*Q3)

add_Q3Q1_int <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q2*Q3)

plot_test('int_Q3.png',width=1000)
 plot(add_Q3, add_Q3Q1_int, add_Q3Q1_int)
dev.off()

plot_test('add_Q3.png',width=1000)
 plot(no_qtl, add_Q3)
dev.off()
################################################################################

################################################################################
no_qtl <- scanone(cross, pheno.col=4, method="hk", model="binary", verbose=FALSE, tol=1e-4, maxit=1000)

AHR.qtl <- makeqtl(cross, chr=1, pos=0, what="prob")

add_Q2_AHR <- addqtl(cross, pheno.col=4, qtl = AHR.qtl, method="hk", model="binary",
            incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1+Q2)

int_Q2_AHR <- addqtl(cross, pheno.col=4, qtl = AHR.qtl, method="hk", model="binary",
            incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2)

plot_test('only_AHR_qtl.png',width=1000)
 plot(no_qtl, add_Q2_AHR, int_Q2_AHR)
dev.off()
################################################################################

int_Q2_AHR <- addqtl(cross, pheno.col=4, qtl = AHR.qtl, method="hk", model="binary",
            incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2)






int_AHR <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
            incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2)

plot_test('only_AHR_qtl.png',width=1000)
 plot(no_qtl, add_AHR, int_AHR)
dev.off()

###########################################################################
### TEST WHETHER ADDED AHR GENOTYPE AS QTL CHANGES GENOME WIDE LOD
add.perms <- scanone(cross, pheno.col=4, model='binary', method = "hk", n.perm = 1000, n.cluster=6)
lod <- summary(add.perms)[2]
add <- scanone(cross, pheno.col=4, model='binary', method = "hk")
qtl <- summary(add,lod)

add.qtl1 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
add.qtl1 <- refineqtl(cross, qtl=add.qtl1, pheno.col=4, model='binary', method = "hk", incl.markers=F)

int.em <- addint(cross, qtl=add.qtl1, formula=y~Q1+Q2, method='hk')

qtls2 <- refineqtl(cross, pheno.col=4, model='binary',
   qtl=add.qtl1, method='hk', incl.markers=F,
   formula=y~Q1*Q2)

qtls2_AHR <- addtoqtl(cross, qtls2, chr=1, pos=0)

int.em <- addint(cross, qtl=qtls2_AHR, formula=y~Q1*Q2+Q3, method='hk')

AHR_noint <- addqtl(cross, pheno.col=4, qtl = qtls2_AHR, method="hk", model="binary",
            incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q3)

AHR_int <- addqtl(cross, pheno.col=4, qtl = qtls2_AHR, method="hk", model="binary",
            incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q3*Q4)



###########################################################################
add_additive <- addqtl(cross, pheno.col=4, qtl = qtls2, method="hk", model="binary",
            incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q3)

AHR_int_Q1 <- addqtl(cross, pheno.col=4, qtl = qtls2, method="hk", model="binary",
            incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q3*Q1)

AHR_int_Q2 <- addqtl(cross, pheno.col=4, qtl = qtls2, method="hk", model="binary",
            incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q3*Q2)

plot_test('elr_missing_AHR.int_no.int.png', width=2000)
plot(add_additive, AHR_int_Q1, AHR_int_Q2, ylab="LOD score", incl.markers=F)
dev.off()
###########################################################################







############################################################
load(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
############################################################

full.norm.imp <- stepwiseqtl(cross, penalties=pens, incl.markers=T, qtl=norm.add.qtls, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8, chr=qtls_chr)


## manual add stepwise qtl
bin.add.hk.perms <- scanone(cross, pheno.col=4, model='binary', method = "hk", n.perm = 10000, n.cluster=10)
lod <- summary(bin.add.hk.perms)[2]
bin.add.hk <- scanone(cross, pheno.col=4, model='binary', method = "hk")
qtl <- summary(bin.add.hk,lod)
bin.add.hk.qtls1 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
bin.add.hk.qtls1 <- refineqtl(cross, qtl=bin.add.hk.qtls1, pheno.col=4, model='binary', method = "hk", incl.markers=F)

int.hk <- addint(cross, qtl=bin.add.hk.qtls1, formula=y~Q1+Q2+Q3, method='hk')
bin.add.hk.qtls2 <- refineqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls1, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
#int.hk <- addint(cross, qtl=bin.add.hk.qtls, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')
##scan for an additional additive QTL
add.hk.a <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls2, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
##scan for an additional interactive QTL
int.hk <- addint(cross, qtl=bin.add.hk.qtls2, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')



add.hk.i <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q1:Q4)

add.hk.2 <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q4)
add.hk.3 <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q3:Q4)

## Scan for interacting pair to add (long)
add.hk.c <- addpair(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)





out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)
