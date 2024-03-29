#!/bin/R

pop <- 'ELR.missing'
library('qtl')
library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
erp <- 0.0025
################################################################################
## Read cross
cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
cross$pheno <- as.data.frame(cross$pheno)
cross <- jittermap(cross, amount=1e-6)
strat <- c("ELR_10977","ELR_10988","ELR_10983","ELR_10980","ELR_10974","ELR_10987","ELR_10971","ELR_10972")
cross$pheno$strat <- ifelse(cross$pheno$ID %in% strat, 1,0)

cross <- sim.geno(cross,step=0,off.end=5, error.prob=erp,map.function="kosambi")
cross <- calc.genoprob(cross,step=1,error.prob=erp,off.end=5)

gg_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross)[[X]], 2)} ))
#gg_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross)[[X]], 0.50)} ))

if(pop == 'ELR.missing') gg_marks <- c(gg_marks,"AHR2a_del")
gg <- pull.markers(cross,gg_marks)
ggmap <- est.map(gg,error.prob=erp,map.function="kosambi",sex.sp=F,n.cluster=6)
gg <- replace.map(gg,ggmap)
gg <- jittermap(gg)
gg <- sim.geno(gg, step=1, error.prob=erp, off.end=5, map.function="kosambi", n.draws=100)
gg <- calc.genoprob(gg, step=1, error.prob=erp, off.end=5, map.function="kosambi")
gg_step2 <- gg
##gg_step2 <- reduce2grid(gg)

norm.add <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = F, max.qtl=5)
norm.add.qtls <- summary(norm.add)
norm.add.qtls <- makeqtl(gg_step2, chr=as.character(norm.add.qtls[['chr']]), pos=as.numeric(norm.add.qtls[['pos']]), what="draws")
qtls_chr <- unique(c(norm.add.qtls[['chr']],1,2,5,8,13,18,23,24))

################################################################################
save.image(file.path(mpath,paste0(pop,'_AHR_cov.rsave')))
################################################################################

perms <- scantwo(gg_step2, chr = c(1:4,6:24), pheno.col=5, model="normal", method="imp",incl.markers=F, clean.output=T, clean.nmar=10,clean.distance=10,n.perm=1000,assumeCondIndep=T,n.cluster=22,perm.strata=gg_step2$pheno$strat)
pens <- calc.penalties(perms, alpha=0.10)

################################################################################
save.image(file.path(mpath,paste0(pop,'_AHR_cov.rsave')))
################################################################################

full.norm.imp <- stepwiseqtl(gg_step2, penalties=pens, incl.markers=T, qtl=norm.add.qtls, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8, chr=qtls_chr)

################################################################################
save.image(file.path(mpath,paste0(pop,'_AHR_cov.rsave')))
################################################################################

################################################################################
## normal
perms.norm.imp.dwnsmpl <- scanone(gg_step2, method = "imp", model = "normal", maxit = 1000, n.perm = 10000, pheno.col = 5, n.cluster = 22)
perms.norm.imp <- scanone(cross, method = "imp", model = "normal", maxit = 1000, n.perm = 10000, pheno.col = 5, n.cluster = 22)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
################################################################################

####################################################################################
print(summary(perms.norm.imp))
################################################################################
save.image(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
################################################################################
