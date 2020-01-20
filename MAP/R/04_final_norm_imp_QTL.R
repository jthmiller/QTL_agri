#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cores <- detectCores() - 2

################################################################################
################################################################################
print(pop)
print('norm imp')
## Error prob = 0.025
if(pop == 'NBH') erp <- 0.0025
if(pop == 'ELR') erp <- 0.0025
if(pop == 'ELR.missing') erp <- 0.0025
################################################################################
################################################################################
## Read cross
cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
cross$pheno <- as.data.frame(cross$pheno)
cross <- jittermap(cross)

cross <- sim.geno(cross, stepwidth="fixed", step=1,off.end=5, error.prob=erp ,map.function="kosambi")
cross <- calc.genoprob(cross, stepwidth="fixed", step=1, error.prob=erp, off.end=5)

gg_step2 <- reduce2grid(cross)

################################################################################

norm.add.imp.perms <- scanone(gg_step2, pheno.col=4, model='normal', method = "imp", n.perm = 2000, n.cluster=cores)
lod <- summary(norm.add.imp.perms)[1]
norm.add.imp <- scanone(gg_step2, pheno.col=4, model='normal', method = "imp")

qtl <- summary(norm.add.imp,lod)
norm.add.imp.qtls <- makeqtl(gg_step2, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
qtls_chr <- unique(c(norm.add.imp.qtls[['chr']],1,2,5,8,13,18,24))

################################################################################
save.image(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
################################################################################

perms <- scantwo(gg_step2, incl.markers=F, chr = c(1:4,6:24), pheno.col=5, model="normal", method="imp", clean.output=T, clean.nmar=10,clean.distance=10,n.perm=1,assumeCondIndep=T,n.cluster=22)
summary(perms)
pens <- calc.penalties(perms, alpha=0.10)
summary(pens)

################################################################################
save.image(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
################################################################################

full.norm.imp <- stepwiseqtl(gg_step2, penalties=pens, incl.markers=F, qtl=norm.add.qtls, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)

################################################################################
save.image(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
################################################################################

################################################################################
full.norm.hk <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = F, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=8)
save.image(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
