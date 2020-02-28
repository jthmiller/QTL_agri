#!/bin/R
pop <- 'ELR'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- as.numeric(commandArgs(TRUE)[2])
cores <- 20
################################################################################

pop <- 'NBH'
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))


#### NBH INTERACTIONS ###########################################################
## summary(norm.em.2, thresholds=c(0, Inf, 6, Inf, Inf), what="int")
        pos1  pos2 lod.full lod.fv1 lod.int lod.add lod.av1
c2:c13 80.59 39.42    27.32    4.68    9.38  17.937  -4.705
c3:c24 33.14 77.79     9.18    5.69    6.41   2.767  -0.727
c4:cX   6.14  8.93     7.19    5.88    6.95   0.242  -1.068
cX:c7  54.49 62.82     7.25    5.53    6.01   1.244  -0.475
cX:c24  1.66 45.50     9.70    6.15    6.91   2.784  -0.760

## summary(scanone(cross, pheno.col=5, model="normal", method="em"))
## 18:17874376            18 39.51  7.017
## 1:191503                1  0.00  2.297
## 2:37728362              2 89.41 22.642
## NW_012224670.1:105636   8 22.28  3.818
## 18:17874376            18 39.51  7.017
## 24:31080447            24 53.23  3.494

nbh_ahra <- ahr_genes[which(ahr_genes$chr == 1 & ahr_genes$gene = 'AHR1'),'close_marker']
nbh_aip <- ahr_genes[which(ahr_genes$chr == 2 & ahr_genes$gene = 'aip'),'close_marker']
nbh_arnt8 <- ahr_genes[which(ahr_genes$chr == 8 & ahr_genes$gene = 'ARNT'),'close_marker']
nbh_ahrb <- ahr_genes[which(ahr_genes$chr == 18 & ahr_genes$gene = 'AHR2b'),'close_marker']
nbh_arnt13 <- ahr_genes[which(ahr_genes$chr == 13 & ahr_genes$gene = 'ARNT'),'close_marker']
nbh_hsp13 <- ahr_genes[which(ahr_genes$chr == 13 & ahr_genes$gene = 'hspa14'),'close_marker']
nbh_hspb8  <- ahr_genes[which(ahr_genes$chr == 24 & ahr_genes$gene = 'hspb8'),'close_marker']


nbh_qtl2 <- find.marker(cross,2,89.41)
nbh_qtl8 <- find.marker(cross,8,22.28)
nbh_qtl13 <- find.marker(cross,13,39.42)
nbh_qtl18 <- find.marker(cross,18,39.51)
nbh_qtl24 <- find.marker(cross,24,52.23)


nbh_qtl4 <- find.marker(cross,4,6.14)
nbh_qtlx <- find.marker(cross,'X',8.93)


chr <- c(2,8,13,18,24)
pos <- c(89.41,22.28,39.42,39.51,52.23)


error <- 0.001
cross <- sim.geno(cross, n.draws=160, error.prob=error, map.function="kosambi", stepwidth="fixed")
cross <- calc.genoprob(cross, error.prob=error, map.function="kosambi", stepwidth="fixed")


###### NBH HK FIT ##################################################################
em.qtl <- makeqtl(cross,chr =chr, pos = pos, what="prob")
##em.qtl <- refineqtl(cross, pheno.col = 5, qtl = em.qtl, method = "hk", model='normal',incl.markers=T)

fit_hk_5wINT_normal <- fitqtl(cross, pheno.col=5, method="hk", model="normal", qtl = em.qtl,
                covar=NULL, formula = y~Q1+Q2+Q3+Q4+Q5+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=5000, forceXcovar=FALSE)

fit_hk_3wINT_binary <- fitqtl(cross, pheno.col=4, method="hk", model="binary", qtl = em.qtl,
                covar=NULL, formula = y~Q1+Q3+Q4+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=5000, forceXcovar=FALSE)

fit_hk_3wINT_normal <- fitqtl(cross, pheno.col=4, method="hk", model="normal", qtl = em.qtl,
                covar=NULL, formula = y~Q1+Q3+Q4+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_hk_3wINT_normal)
summary(fit_hk_3wINT_binary)
summary(fit_hk_5wINT_normal)
################################################################

###### NBH IMP FIT ##################################################################

imp.qtl <- makeqtl(cross,chr =chr, pos = pos, what="draws")
### DONT REFINE QTL (moves interaction at 13)

fit_imp_5wINT_normal <- fitqtl(cross, pheno.col=5, method="imp", model="normal", qtl = imp.qtl,
                covar=NULL, formula = y~Q1+Q2+Q3+Q4+Q5+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_5wINT_normal)

fit_imp_3wINT_normal <- fitqtl(cross, pheno.col=5, method="imp", model="normal", qtl = imp.qtl,
                covar=NULL, formula = y~Q1+Q3+Q4+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_3wINT_normal)

fit_imp_3wINT_bin <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = imp.qtl,
                covar=NULL, formula =y~Q1+Q3+Q4+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_3wINT_bin)

fit_imp_5wINT_bin <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = imp.qtl,
                covar=NULL, formula = y~Q1+Q2+Q3+Q4+Q5+Q1:Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_5wINT_bin)


### BEST NORMAL MODEL:
nbh_norm <- summary(fit_imp_3wINT_bin)

################################################################################

pop <- 'ELR'
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
load(file.path(mpath,paste0(pop,'_scan2_normal_em.rsave')))

ahr_genes <- get_AHR(cross)
gt <- geno.table(cross)
ahr_genes$segdist <- -log10(gt[ahr_genes$close_marker,'P.value'])
ahr_genes_sub <- ahr_genes[!is.na(ahr_genes$PATH),]

sone <- scanone(cross, pheno.col=4, model="binary", method="em")
summary(sone)
8:37635736    8 1.02e+02 3.056
18:20273448  18 7.47e+01 5.352
13:4070714   13 1.67e+01 4.696

### ELR
#### ELR INTERACTIONS ###########################################################
## summary(norm.em.2, thresholds=c(0, Inf, 6, Inf, Inf), what="int")

         pos1  pos2 lod.full lod.fv1 lod.int lod.add lod.av1
c1 :c1  90.31 123.4     7.97    5.47    6.25   1.718  -0.778
c4 :cX  57.16   6.9     6.69    5.21    6.04   0.658  -0.825
cX :c20 32.37  31.3     7.05    5.57    6.12   0.934  -0.549
c9 :c12 51.22  75.8     6.62    4.60    6.09   0.523  -1.493
c15:c18  8.82  62.0    10.82    5.16    7.11   3.709  -1.948
################################################################################

elr_ahra <- ahr_genes[which(ahr_genes$chr == 1 & ahr_genes$gene == 'AHR1'),'close_marker']
elr_aip <- ahr_genes[which(ahr_genes$chr == 2 & ahr_genes$gene == 'aip'),'close_marker']
elr_arnt8 <- ahr_genes[which(ahr_genes$chr == 8 & ahr_genes$gene == 'ARNT'),'close_marker']
elr_ahrb <- ahr_genes[which(ahr_genes$chr == 18 & ahr_genes$gene == 'AHR2b'),'close_marker']
elr_arnt13 <- ahr_genes[which(ahr_genes$chr == 13 & ahr_genes$gene == 'ARNT'),'close_marker']
elr_hsp13 <- ahr_genes[which(ahr_genes$chr == 13 & ahr_genes$gene == 'hspa14'),'close_marker']
elr_hspb8  <- ahr_genes[which(ahr_genes$chr == 24 & ahr_genes$gene == 'hspb8'),'close_marker']

elr_qtl8 <- find.marker(cross,8,10.2)
elr_qtl13 <- find.marker(cross,13,16.7)
elr_qtl18 <- find.marker(cross,18,74.7)
elr_qtl15 <- find.marker(cross,15,8.82)

elr.add.qtl <- c(elr_qtl8,elr_qtl13,elr_qtl18,elr_qtl15)

chr <- c(8,13,15,18)
pos <- c(10.2,16.7,74.7,8.82)


################################################################################
em.qtl <- makeqtl(cross,chr =chr, pos = pos, what="prob")
##em.qtl <- refineqtl(cross, pheno.col = 5, qtl = em.qtl, method = "hk", model='normal',incl.markers=T)

fit_hk_3wINT_normal <- fitqtl(cross, pheno.col=5, method="hk", model="normal", qtl = em.qtl,
                covar=NULL, formula = y~Q2+Q3+Q3:Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=5000, forceXcovar=FALSE)
summary(fit_hk_3wINT_normal)


fit_hk_3wINT_binary <- fitqtl(cross, pheno.col=4, method="hk", model="binary", qtl = em.qtl,
                covar=NULL, formula = y~Q2+Q3+Q3:Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=5000, forceXcovar=FALSE)
summary(fit_hk_3wINT_binary)
################################################################################


imp.qtl <- makeqtl(cross,chr = chr, pos = pos, what="draws")
### DONT REFINE QTL (moves interaction at 13)

fit_imp_5wINT_normal <- fitqtl(cross, pheno.col=5, method="imp", model="normal", qtl = imp.qtl,
                covar=NULL, formula = y~Q1+Q2+Q3+Q4+Q3:Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_5wINT_normal)

fit_imp_3wINT_normal <- fitqtl(cross, pheno.col = 5, method="imp", model="normal", qtl = imp.qtl,
                covar=NULL, formula = y~Q2+Q3+Q3:Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_3wINT_normal)

fit_imp_3wINT_normal <- fitqtl(cross, pheno.col = 4, method="imp", model="binary", qtl = imp.qtl,
                covar=NULL, formula = y~Q2+Q3+Q3:Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
summary(fit_imp_3wINT_normal)


################################################################################
elr.add.qtl <- c(elr_qtl8,elr_qtl13,elr_qtl18,elr_qtl15)
sapply(elr.add.qtl, single, pop = 'ELR')

intxs.bin(loc_b = elr_qtl15,loc_a = elr_qtl18, popchr = 'ELR',locbN = 'CYP1b1', main = 'CYP1b1 v AHRb')
intxs.bin(loc_b = elr_aip,loc_a = elr_arnt13, popchr = 'ELR',locbN = 'CYP1b1', main = 'CYP1b1 v AHRb')


elr_aip


nbh.add.qtl <- c(nbh_qtl2,nbh_qtl8,nbh_qtl13,nbh_qtl18,nbh_qtl24)



#### NBH INCOMPATABILITY AT 13 v 2
#### WHAT ELSE IS 13 incompatable with?
c2 :c13 80.59 39.42    27.32    4.68    9.38  17.937 -4.7046
################################################################################

(file.path(mpath,paste0(pop,'_scan2_normal_em.rsave')))
summary(bin.em.2, thresholds=c(0, Inf, 9.1, Inf, Inf), what="int")
###########

summary(bin.em.2, thresholds=c(0, Inf, 9.1, Inf, Inf), what="int")
