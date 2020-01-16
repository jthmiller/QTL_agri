#!/bin/R
### binary trait with imputation

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
################################################################################
## Read cross
cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross,step=0,off.end=5, error.prob=0.025,map.function="kosambi")
cross <- calc.genoprob(cross,step=1,error.prob=0.025,off.end=5)

gg_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross)[[X]], 0.50)} ))
gg <- pull.markers(cross,gg_marks)
gg <- sim.geno(gg, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=100)
gg <- calc.genoprob(gg, step=1, error.prob=0.025, off.end=5, map.function="kosambi")
gg_step2 <- reduce2grid(gg)
################################################################################
################################################################################
################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_imp.rsave')))
################################################################################
################################################################################
bin.add.imp <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = F, max.qtl=5)
bin.add.imp.qtls <- summary(bin.add.imp)
bin.add.imp.qtls <- makeqtl(gg_step2, chr=as.character(bin.add.imp.qtls$chr), pos=as.numeric(bin.add.imp.qtls$pos), what="draws")
qtls_chr <- unique(c(bin.add.imp.qtls$chr,1,2,5,8,13,18,24))
full.bin.imp <- stepwiseqtl(gg_step2, incl.markers=F, qtl=bin.add.imp.qtls, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=5, chr=qtls_chr)

## binary imp is not supported in rQTL
grid.perms.bin.em <- scanone(gg_step2, method = "em", model = "binary", maxit = 1000, n.perm = 10000, pheno.col = 4, n.cluster = 10)
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_imp.rsave')))
################################################################################

################################################################################
cross <- sim.geno(cross,step=0,off.end=5, error.prob=0.025,map.function="kosambi")
cross <- calc.genoprob(cross,step=1,error.prob=0.025,off.end=5)

## binary imp is not supported in rQTL
## binary
scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)
##############################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_imp.rsave')))
################################################################################
