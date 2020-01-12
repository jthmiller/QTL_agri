#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

plotpub <- function(X) { png(paste0('~/public_html/',X,'.png')) }

################################################################################
print(pop)
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

## normal
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
norm.imp.qtl <- makeqtl(cross, chr=c(2,13,18), pos=c(30,33,28),  what="draws")
################################################################################

perms.norm.imp <- scanone(cross, method = "imp", model = "normal", maxit = 1000,
  n.perm = 10000, pheno.col = 5, n.cluster = 10)
####################################################################################
print(summary(perms.norm.imp))
################################################################################
save.image(file.path(mpath,paste0(pop,'.rsave'))
################################################################################
