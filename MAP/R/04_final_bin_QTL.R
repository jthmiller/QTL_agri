#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

plotpub <- function(X) { png(paste0('~/public_html/',X,'.png')) }

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

gg <- sim.geno(cross, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=160)
gg <- calc.genoprob(gg, step=1, error.prob=0.025, off.end=5, map.function="kosambi")
gg_step2 <- reduce2grid(gg)
################################################################################

norm.add <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = F, max.qtl=5)
norm.add.qtls <- summary(norm.add)
norm.add.qtls <- makeqtl(gg_step2, chr=as.character(norm.add.qtls$chr), pos=as.numeric(norm.add.qtls$pos), what="draws")
full.norm.imp <- stepwiseqtl(gg_step2, incl.markers=F, qtl=norm.add.qtls, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=5, chr=c(1,2,5,8,13,18,24))

bin.add.imp <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = F, max.qtl=5)
bin.add.imp.qtls <- summary(bin.add.imp)
bin.add.imp.qtls <- makeqtl(gg_step2, chr=as.character(bin.add.imp.qtls$chr), pos=as.numeric(bin.add.imp.qtls$pos), what="draws")
full.bin.imp <- stepwiseqtl(gg_step2, incl.markers=F, qtl=bin.add.imp.qtls, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=5, chr=c(1,2,5,8,13,18,24))

bin.add.hk <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = T, model='binary', method = "hk", pheno.col = 4, scan.pairs = F, max.qtl=5)
bin.add.hk.qtls <- summary(bin.add.hk)
bin.add.hk.qtls <- makeqtl(gg_step2, chr=as.character(bin.add.hk.qtls$chr), pos=as.numeric(bin.add.hk.qtls$pos), what="draws")
full.bin.hk <- stepwiseqtl(gg_step2, incl.markers=F, qtl=bin.add.hk.qtls, additive.only = F, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=5, chr=c(1,2,5,8,13,18,24)

################################################################################
## binary
scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)
bins <- data.frame(em=summary(scan.bin.em),mr=summary(scan.bin.mr))

## normal
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
norm.imp.qtl <- makeqtl(cross, chr=c(2,13,18), pos=c(30,33,28),  what="draws")
################################################################################


####################################################################################
## PERMS WITH ALL LOCI
perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 1000,
  n.perm = 10000, pheno.col = 4, n.cluster = 10)

perms.norm.imp <- scanone(cross, method = "imp", model = "normal", maxit = 1000,
  n.perm = 10000, pheno.col = 5, n.cluster = 10)
####################################################################################
#LOD thresholds (1000 permutations)
#     lod
#5%  4.20
#10% 3.84
# LOD thresholds (1000 permutations)
#      lod
# 5%  4.17
# 10% 3.83
################################################################################
save.image(file.path(mpath,paste0(pop,'.rsave'))
################################################################################
