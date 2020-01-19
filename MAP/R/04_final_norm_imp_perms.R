#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
perm_count <- as.numeric(commandArgs(TRUE)[!is.na(as.numeric(commandArgs(TRUE)))])

print(commandArgs(TRUE))
print(paste(pop,perm_count))

library('qtl')
library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cores <- detectCores() - 2

################################################################################
if(pop == 'NBH') erp <- 0.0025
if(pop == 'ELR') erp <- 0.0025
if(pop == 'ELR.missing') erp <- 0.0025
################################################################################

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
cross$pheno <- as.data.frame(cross$pheno)

################################################################################

dups <- findDupMarkers(cross, exact.only=F, adjacent.only=F)
cross <- pull.markers(cross, names(dups))
cross <- sim.geno(cross, stepwidth="fixed", step=1, error.prob=erp, off.end=5, map.function="kosambi", n.draws=100)
cross <- calc.genoprob(cross, stepwidth="fixed", step=1, error.prob=erp, off.end=5, map.function="kosambi")

################################################################################

perms <- scantwo(cross, incl.markers=F, chr = c(1:4,6:24), pheno.col=5, model="normal", method="imp", clean.output=T, clean.nmar=10,clean.distance=10,n.perm=perm_count, assumeCondIndep=T, n.cluster=cores)
pens <- calc.penalties(perms, alpha=0.05)
print(paste('done with', perm_count, 'permutations'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################

norm.add <- stepwiseqtl(cross, incl.markers=F, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = F, max.qtl=8)
norm.add.qtls <- summary(norm.add)
norm.add.qtls <- makeqtl(cross, chr=as.character(norm.add.qtls[['chr']]), pos=as.numeric(norm.add.qtls[['pos']]), what="draws")
qtls_chr <- unique(c(norm.add.qtls[['chr']],1,2,5,8,13,18,23,24))
full.norm.imp <- stepwiseqtl(cross, penalties=pens, incl.markers=F, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8, chr=qtls_chr, qtl=norm.add.qtls)
print(paste('done with', perm_count, 'full.norm.imp'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################

full.norm.hk <- stepwiseqtl(cross, penalties=pens, incl.markers=F, additive.only = F, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=8, chr=qtls_chr, qtl=norm.add.qtls)
print(paste('done with', perm_count, 'full.norm.hk'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################
