#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
##cores <- detectCores() - 2
cores <- 22
################################################################################
################################################################################
print(pop)
print('binary imp')
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
cross <- jittermap(cross)
dups <- findDupMarkers(cross, exact.only=F, adjacent.only=F)
if(pop == 'ELR.missing') dups <- c(dups,"AHR2a_del")
cross <- pull.markers(cross, names(dups))
cross$pheno <- as.data.frame(cross$pheno)
cross <- jittermap(cross)
cross <- sim.geno(cross, stepwidth="fixed", step=1,off.end=5, error.prob=erp ,map.function="kosambi", n.draws=100)
cross <- calc.genoprob(cross, stepwidth="fixed", step=1, error.prob=erp, off.end=5)

################################################################################
################################################################################
#gg_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross)[[X]], 2)} ))
#if(pop == 'ELR.missing') gg_marks <- c(gg_marks,"AHR2a_del")
#gg <- pull.markers(cross,gg_marks)
#ggmap <- est.map(gg,error.prob=erp,map.function="kosambi",sex.sp=F,n.cluster=6)
#gg <- replace.map(gg,ggmap)
#gg <- jittermap(gg)
#gg <- sim.geno(gg, step=1, error.prob=erp, off.end=5, map.function="kosambi", n.draws=100)
#gg <- calc.genoprob(gg, step=1, error.prob=erp, off.end=5, map.function="kosambi")
#gg_step2 <- gg
################################################################################
gg_step2 <- reduce2grid(cross)
################################################################################

bin.add.hk.perms <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk", n.perm = 2000,n.cluster=cores)
lod <- summary(bin.add.hk.perms)[1]
bin.add.hk <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk")

qtl <- summary(bin.add.hk,lod)
bin.add.hk.qtls <- makeqtl(gg_step2, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
qtls_chr <- unique(c(bin.add.hk.qtls[['chr']],1,2,5,8,13,18,24))

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_hk.rsave')))
################################################################################

perms <- scantwo(gg_step2, incl.markers=F, chr = c(1:4,6:24), pheno.col=4, model="binary", method="hk", clean.output=T, clean.nmar=10,clean.distance=10,n.perm=22, assumeCondIndep=T, n.cluster=cores)
pens <- calc.penalties(perms, alpha=0.05)

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_hk.rsave')))
################################################################################

full.bin.hk <- stepwiseqtl(gg_step2,  penalties=pens, incl.markers=F, qtl=bin.add.hk.qtls, additive.only = F, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=8)

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_hk.rsave')))
################################################################################
