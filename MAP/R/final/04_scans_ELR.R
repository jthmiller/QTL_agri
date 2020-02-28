#!/bin/R
pop <- 'ELR'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- as.numeric(commandArgs(TRUE)[2])
cores <- 20
################################################################################
################################################################################
#fl <- paste0(pop,'_',sum(nmar(cross)),'_imputed_high_confidence_tsp_mapped.csv')
fl <- 'ELR_7457_imputed_high_confidence_tsp_mapped.csv'
cross <- read.cross(file=fl , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

cross$pheno <- as.data.frame(cross$pheno)
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'
cross$pheno$pheno_norm <- nqrank(cross$pheno$Pheno)
################################################################################

##error <- 1e-04
error <- 0.001
cross <- sim.geno(cross,n.draws=160, error.prob=error, map.function="kosambi", stepwidth="fixed")
cross <- calc.genoprob(cross, error.prob=error, map.function="kosambi", stepwidth="fixed")


















### HK
sone <- scanone(cross, pheno.col=4, model="binary", method="hk")
sone.perms <- scanone(cross, pheno.col=4, model="binary", method="hk", n.perm=1000, n.cluster=cores)
lod <- summary(sone.perms)[[2]]
qtl <- summary(sone,lod)

### MAKE COVARIATES

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

add1 <- summary(sone.int1)[which.max(summary(sone.int1)$lod),]
add2 <- summary(sone.int2)[which.max(summary(sone.int2)$lod),]
int1 <- summary(sone.1)[which.max(summary(sone.1)$lod),]
int2 <- summary(sone.2)[which.max(summary(sone.2)$lod),]

qtl <- rbind(add1,add2,int1,int2)

hk.qtl.4 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
hk.qtl.5 <-  addtoqtl(cross, qtl = hk.qtl.4, chr = 15, pos =  8.82)

int.hk <- addint(cross, pheno.col = 5, qtl = hk.qtl.5, method='hk', model='normal', formula=y~Q1+Q2+Q3+Q4+Q5, maxit=1000)
## 13*18
## 18*23

fit_hk_5int <- fitqtl(cross, pheno.col=1, method="hk", model="normal", qtl = hk.qtl.5,
                covar=NULL, formula = y~Q1+Q2*Q5+Q3+Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)


##hk.qtl.5 <- refineqtl(cross, pheno.col = 5, qtl=hk.qtl.5, method = "hk", model='normal',incl.markers=T)


add_Q5_hk <- addqtl(cross, pheno.col = 5, qtl = hk.qtl.4, method = "hk", model = "normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1+Q2+Q3:Q4+Q5)

add_Q5_hk <- addqtl(cross, pheno.col = 5, qtl = hk.qtl.4, method = "hk", model = "normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1+Q2*Q5)


mar <- find.marker(cross,qtl.int$chr[1],qtl.int$pos[1])
g.int2 <- pull.geno(fill.geno(cross))[,mar]
g.int2 <- data.frame(cbind(as.numeric(g.int2 == 1), as.numeric(g.int2 == 2)))




sone.int1 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.int)
sone.int1.perms <- scanone(cross, pheno.col=5, model="normal", method="hk", n.perm=1000, n.cluster=cores, intcovar=g.int)
lod <- summary(sone.int1.perms)[[2]]
qtl.int2 <- summary(sone.int1,lod)

sone.int2 <- scanone(cross, pheno.col=5, model="normal", method="hk", intcovar=g.int)
sone.int2.perms <- scanone(cross, pheno.col=5, model="normal", method="hk", n.perm=1000, n.cluster=cores, intcovar=g.int)
lod <- summary(sone.int2.perms)[[2]]
qtl.int2 <- summary(sone.int2,lod)

mara <- find.marker(cross,2,58.55)
g.add <- pull.geno(fill.geno(cross))[,mara]
g.add <- data.frame(cbind(as.numeric(g.add == 1), as.numeric(g.add == 2)))

marb <- find.marker(cross,13,30.5)
g.int <- pull.geno(fill.geno(cross))[,marb]
g.int <- data.frame(cbind(as.numeric(g.int == 1), as.numeric(g.int == 2)))

sone.int.add <- scanone(cross, pheno.col=5, model="normal", method="hk", addcovar=g.add, intcovar=g.int)
sone.int.add.perms <- scanone(cross, pheno.col=5, model="normal", method="hk", n.perm=1000, n.cluster=cores, addcovar=g.add, intcovar=g.int)
lod <- summary(sone.int.add.perms)[[2]]
qtl.int <- summary(sone.int.add,lod)

sone.add <- scanone(cross, pheno.col=5, model="normal", method="hk", addcovar=g.add)
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

#### HK ##########################################################################################
sone <- scanone(cross, pheno.col=4, model="binary", method="hk")
sone.perms <- scanone(cross, pheno.col=4, model="binary", method="hk", n.perm=1000, n.cluster=cores, addcovar=g.em)
summary(sone, alpha=0.1, lodcolumn=1, pvalues=T, perms=sone.perms, ci.function="bayesint")
lod <- summary(sone.perms)[[2]]
qtl <- summary(sone,lod)

hk.qtl.2 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")

hk.qtl.2 <- refineqtl(cross, pheno.col = 4, qtl=hk.qtl.2, method = "hk", model='binary',incl.markers=T)

int.hk <- addint(cross, pheno.col = 4, qtl = hk.qtl.2, method='hk', model='binary',
                 formula=y~Q1+Q2, maxit=1000)

add_Q3_hk <- addqtl(cross, pheno.col=4, qtl = hk.qtl.2, method="hk", model="binary",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1*Q2+Q3)

add <- summary(add_Q3_hk)
ind <- which.max(add$lod)

hk.qtl.3 <-  addtoqtl(cross, qtl = hk.qtl.2, chr = as.character(add[ind,'chr']), pos =  add[ind,'pos'])

hk.qtl.3 <- refineqtl(cross, pheno.col = 4, qtl=hk.qtl.3, method = "hk", model='normal',incl.markers=T)

int.hk <- addint(cross, pheno.col = 4, qtl = hk.qtl.3, method='hk', model='binary',
                 formula=y~Q1*Q2+Q3, maxit=1000)

fit_hk_3int <- fitqtl(cross, pheno.col=4, method="hk", model="binary", qtl = hk.qtl.3,
                covar=NULL, formula = y~Q1*Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

summary(int.hk)

add_Q4_hk <- addqtl(cross, pheno.col=4, qtl = hk.qtl.3, method="hk", model="binary",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1*Q2+Q3+Q4)

add <- summary(add_Q4_hk)
ind <- which.max(add$lod)

hk.qtl.4 <-  addtoqtl(cross, qtl = hk.qtl.3, chr = as.character(add[ind,'chr']), pos =  add[ind,'pos'])
hk.qtl.4 <- refineqtl(cross, pheno.col = 4, qtl=hk.qtl.4, method = "hk", model='normal',incl.markers=T)



summary(fit_hk_3int)

add_Q4_hk <- addqtl(cross, pheno.col=4, qtl = hk.qtl.3, method="hk", model="binary",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1*Q2+Q3+Q4, covar = as.data.frame(cross$pheno$sex))

add <- summary(add_Q4_hk)

hk.qtl.3 <-  addtoqtl(cross, qtl = hk.qtl.3, chr = as.character(add[16,'chr']), pos =  add[16,'pos'])

hk.qtl.3 <- refineqtl(cross, pheno.col = 4, qtl=hk.qtl.3, method = "hk", model='normal',incl.markers=T)

fit_hk_3int.sex <- fitqtl(cross, pheno.col=4, method="hk", model="binary", qtl = hk.qtl.3,
                formula = y~Q1*Q2+Q3+Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE, covar = as.data.frame(cross$pheno$sex))

summary(fit_hk_3int.sex)
################################################################################

################################################################################
################################################################################
#### HK NORMAL #################################################################
sone <- scanone(cross, pheno.col = 5, model="normal", method="hk")
sone.perms <- scanone(cross, pheno.col = 5, model="normal", method="hk", n.perm=1000, n.cluster=cores, addcovar=g.em)
summary(sone, alpha=0.1, lodcolumn = 1, pvalues=T, perms=sone.perms, ci.function="bayesint")
lod <- summary(sone.perms)[[2]]
qtl <- summary(sone,lod)

hk.qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
hk.qtl <- refineqtl(cross, pheno.col = 5, qtl=hk.qtl, method = "hk", model='normal',incl.markers=T)

int.hk <- addint(cross, pheno.col = 5, qtl = hk.qtl, method='hk', model='normal',
                 formula=y~Q1+Q2, maxit=1000)

summary(int.hk)

add_Q3_hk <- addqtl(cross, pheno.col = 5, qtl = hk.qtl, method="hk", model="normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1*Q2+Q3)

summary(add_Q3_hk)

fit_hk_2int <- fitqtl(cross, pheno.col = 5, method="hk", model="normal", qtl = hk.qtl,
                covar=NULL, formula = y~Q1*Q2, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

summary(fit_hk_2int)

add_Q3_hk.sex <- addqtl(cross, pheno.col = 5, qtl = hk.qtl, method="hk", model="normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1*Q2+Q3, covar = as.data.frame(cross$pheno$sex))

summary(add_Q3_hk.sex)
##########################################################################################

##########################################################################################
##########################################################################################
##########################################################################################
#### IMP ##########################################################################################
imp <- scanone(cross, pheno.col = 5, model="normal", method="imp")
imp.perms <- scanone(cross,pheno.col = 5, model="normal", method="imp", n.perm=100, n.cluster=cores, addcovar=g.imp)
summary(imp, alpha=0.1, lodcolumn = 1, pvalues=T, perms=imp.perms, ci.function="bayesint")

lod <- summary(imp.perms)[[2]]
qtl <- summary(imp)[c(13,18),]
qtl

imp.qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")

imp.qtl <- refineqtl(cross, pheno.col = 5, qtl=imp.qtl, method = "imp", model='normal',incl.markers=T)

int.em.sex <- addint(cross, pheno.col = 5, qtl = imp.qtl, method='imp', model='normal',
                 covar=data.frame(cross$pheno$sex) ,formula=y~Q1+Q2, maxit=1000)
summary(int.em.sex)

int.em <- addint(cross, pheno.col = 5, qtl = imp.qtl, method='imp', model='normal',
                 formula=y~Q1+Q2, maxit=1000)
summary(int.em)

add_Q3_imp <- addqtl(cross, pheno.col = 5, qtl = imp.qtl, method="imp", model="normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1+Q2+Q3)

add <- summary(imp)[23,]

imp.qtl.3 <-  addtoqtl(cross, qtl = imp.qtl, chr = as.character(add[23,'chr']), pos =  add[23,'pos'])

imp.qtl.3 <- refineqtl(cross, pheno.col = 5, qtl=imp.qtl.3, method = "imp", model='normal',incl.markers=T)

fit_imp_3int <- fitqtl(cross, pheno.col = 5, method="imp", model="normal", qtl = imp.qtl.3,
                covar=NULL, formula=y~Q1+Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

summary(fit_imp_3int)

add_Q4_imp <- addqtl(cross, pheno.col = 5, qtl = imp.qtl.3, method="imp", model="normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1+Q2+Q3+Q4)

add <- summary(add_Q4_imp)

imp.qtl.3 <-  addtoqtl(cross, qtl = imp.qtl.3, chr = as.character(add[2,'chr']), pos =  add[2,'pos'])
imp.qtl.3 <- refineqtl(cross, pheno.col = 5, qtl=imp.qtl.3, method = "imp", model='normal',incl.markers=T)

fit_imp_3int <- fitqtl(cross, pheno.col = 5, method="imp", model="normal", qtl = imp.qtl.3,
                covar=NULL, formula=y~Q1+Q2+Q3+Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

summary(fit_imp_3int)

int.imp <- addint(cross, pheno.col = 5, qtl = imp.qtl.3, method='imp', model='normal',
                 formula=y~Q1+Q2+Q3+Q4, maxit=1000)

summary(int.imp)

fit_imp_3int <- fitqtl(cross, pheno.col = 5, method="imp", model="normal", qtl = imp.qtl.3,
                covar=NULL, formula=y~Q1+Q2*Q3+Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

summary(fit_imp_3int)

add_Q4_imp <- addqtl(cross, pheno.col = 5, qtl = imp.qtl.3, method="imp", model="normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula=y~Q1+Q2*Q3+Q4)

summary(add_Q4_imp)

add_Q4_imp.sex <- addqtl(cross, pheno.col = 5, qtl = imp.qtl.3, method="imp", model="normal",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula=y~Q1+Q2*Q3+Q4, covar = as.data.frame(cross$pheno$sex))

summary(add_Q4_imp.sex)

add <- summary(add_Q4_imp)
imp.qtl.4 <-  addtoqtl(cross, qtl = imp.qtl.3, chr = as.character(add[24,'chr']), pos =  add[24,'pos'])
imp.qtl.4 <- refineqtl(cross, pheno.col = 5, qtl=imp.qtl.4, method = "imp", model='normal',incl.markers=T)

fit_imp_5int <- fitqtl(cross, pheno.col = 5, method="imp", model="normal", qtl = imp.qtl.4,
                covar=NULL, formula=y~Q1+Q2*Q3+Q4+Q5, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)
summary(fit_imp_5int)

int.imp <- addint(cross, pheno.col = 5, qtl = imp.qtl.4, method='imp', model='normal',
                 formula=y~Q1+Q2*Q3+Q4+Q5, maxit=1000)
summary(int.imp)
##########################################################################################
##########################################################################################
##########################################################################################

hk.norm.out <- summary(fit_hk_3int.sex)
imp.norm.out <- summary(fit_imp_5int)

save.image(file.path(mpath,paste0(pop,'_imputed.rsave')))

full.bin.em.step <- stepwiseqtl(cross, model='normal', method = "imp", pheno.col = 5,
 incl.markers=T, qtl=imp.qtl.4, additive.only = T,  scan.pairs = T, max.qtl=8)


save.image(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
