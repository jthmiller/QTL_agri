#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
perm_count <- as.numeric(commandArgs(TRUE)[!is.na(as.numeric(commandArgs(TRUE)))])

print(commandArgs(TRUE))
print(paste(pop,perm_count))

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cores <- 22
print(paste(cores,'cores'))

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
cross <- jittermap(cross)
################################################################################

dups <- findDupMarkers(cross, exact.only=F, adjacent.only=F)
if(pop == 'ELR.missing') dups <- c(dups,"AHR2a_del")
cross <- pull.markers(cross, names(dups))
cross <- sim.geno(cross, stepwidth="fixed", step=1, error.prob=erp, off.end=5, map.function="kosambi", n.draws=100)
#cross <- calc.genoprob(cross, stepwidth="fixed", step=1, error.prob=erp, off.end=5, map.function="kosambi")

################################################################################

perms <- scantwo(cross, incl.markers=F, chr = c(1:4,6:24), pheno.col=5, model="normal", method="imp", clean.output=T, clean.nmar=10,clean.distance=10,n.perm=perm_count, assumeCondIndep=T, n.cluster=cores)
summary(perms)
pens <- calc.penalties(perms, alpha=0.05)
summary(pens)
print(pens)
print(paste('done with', perm_count, 'permutations'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################

norm.add.imp.perms <- scanone(cross, pheno.col=5, model='normal', method = "imp", n.perm = 2000, n.cluster=cores)
lod <- summary(norm.add.imp.perms)[2]
norm.add.imp <- scanone(cross, pheno.col=5, model='normal', method = "imp")

qtl <- summary(norm.add.imp,lod)
norm.add.imp.qtls <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################

full.norm.hk <- stepwiseqtl(cross, penalties=pens, incl.markers=F, additive.only = F, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=8, qtl=norm.add.imp.qtls)
print(paste('done with', perm_count, 'full.norm.hk'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################
