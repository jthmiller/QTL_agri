#!/bin/R
pop <- 'ELR'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- as.numeric(commandArgs(TRUE)[2])

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

################################################################################
## COVARIATES
cov <- ifelse(pop == 'ELR',18,2)

imp <- summary(scanone(cross, pheno.col=4, model="normal", method="imp"))[cov,]
mar.imp <- find.marker(cross, imp$chr, imp$pos)
g <- pull.geno(fill.geno(cross))[,mar.imp]
g.imp <- cbind(as.numeric(g==1), as.numeric(g==2))
summary(scanone(cross,pheno.col=4, model="normal", method="imp",addcovar=g.imp))

em <- summary(scanone(cross, pheno.col=4, model="bin", method="em"))[cov,]
mar.em <- find.marker(cross, em$chr, em$pos)
g <- pull.geno(fill.geno(cross))[,mar.em]
g.em <- cbind(as.numeric(g==1), as.numeric(g==2))
summary(scanone(cross,pheno.col=4, model="bin", method="em",addcovar=g.em))

mar.dist <- '17:14629450'
g.dist <- pull.geno(fill.geno(cross))[,'17:14629450']
g.dist <- cbind(as.numeric(g.dist==1), as.numeric(g.dist==2))
summary(scanone(cross,pheno.col=4, model="bin", method="em",addcovar=g.dist))
################################################################################

#### HK ##########################################################################################
sone <- scanone(cross, pheno.col=4, model="binary", method="hk")
sone.perms <- scanone(cross, pheno.col=4, model="binary", method="hk", n.perm=1000, n.cluster=cores, addcovar=g.em)
summary(sone, alpha=0.1, lodcolumn=1, pvalues=T, perms=sone.perms, ci.function="bayesint")
lod <- summary(sone.perms)[[2]]
qtl <- summary(sone,lod)

hk.qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")

hk.qtl <- refineqtl(cross, pheno.col = 4, qtl=hk.qtl, method = "hk", model='binary',incl.markers=T)

int.hk <- addint(cross, pheno.col = 4, qtl = hk.qtl, method='hk', model='binary',
                 formula=y~Q1+Q2, maxit=1000)

add_Q3_hk <- addqtl(cross, pheno.col=4, qtl = hk.qtl, method="hk", model="binary",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1*Q2+Q3)

add <- summary(add_Q3_hk)
ind <- which.max(add$lod)

hk.qtl.3 <-  addtoqtl(cross, qtl = hk.qtl, chr = as.character(add[ind,'chr']), pos =  add[ind,'pos'])

hk.qtl.3 <- refineqtl(cross, pheno.col = 4, qtl=hk.qtl.3, method = "hk", model='normal',incl.markers=T)

int.hk <- addint(cross, pheno.col = 4, qtl = hk.qtl.3, method='hk', model='binary',
                 formula=y~Q1*Q2+Q3, maxit=1000)

summary(int.hk)

add_Q4_hk <- addqtl(cross, pheno.col=4, qtl = hk.qtl.3, method="hk", model="binary",
                    incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
                    formula = y~Q1*Q2+Q3+Q4)

summary(add_Q4_hk)

fit_hk_3int <- fitqtl(cross, pheno.col=4, method="hk", model="binary", qtl = hk.qtl.3,
                covar=NULL, formula = y~Q1*Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

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
