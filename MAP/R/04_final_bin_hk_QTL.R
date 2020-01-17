#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
################################################################################
################################################################################
print(pop)
print('binary imp')
## Error prob = 0.025
if(pop == 'NBH') erp <- 0.0025
if(pop == 'ELR') erp <- 0.0025
if(pop == 'ELR_Missing') erp <- 0.0025
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
cross <- jittermap(cross, amount=1e-6)

cross <- sim.geno(cross,step=1, stepwidth="fixed", off.end=5, error.prob=erp,map.function="kosambi")
cross <- calc.genoprob(cross,step=1,stepwidth="fixed", error.prob=erp, off.end=5)
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
bin.add.hk <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = T, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=5)
bin.add.hk.qtls <- summary(bin.add.hk)
bin.add.hk.qtls <- makeqtl(gg_step2, chr=bin.add.hk.qtls[['chr']], pos=bin.add.hk.qtls[['pos']], what="prob")
qtls_chr <- unique(c(bin.add.hk.qtls[['chr']],1,2,5,8,13,18,24))
################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_hk.rsave')))
################################################################################

perms <- scantwo(gg_step2, incl.markers=F, chr = c(1:4,6:24), pheno.col=4, model="binary", method="imp", clean.output=T, clean.nmar=10,clean.distance=10,n.perm=1000,assumeCondIndep=T,n.cluster=22)
pens <- calc.penalties(perms, alpha=0.10)

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_hk.rsave')))
################################################################################

full.bin.imp <- stepwiseqtl(gg_step2,  penalties=pens, incl.markers=F, qtl=bin.add.imp.qtls, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=5, chr=qtls_chr)

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_hk.rsave')))
################################################################################

full.bin.hk <- stepwiseqtl(gg_step2, penalties=pens, incl.markers=F, qtl=bin.add.hk.qtls, additive.only = F, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=5, chr=qtls_chr)

################################################################################
## binary
scan.bin.hk <- scanone(cross, method = "hk", model = "binary", pheno.col = 4)
scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)
scan.bin.em <- scanone(gg_step2, method = "em", model = "binary", pheno.col = 4)
perms.bin.em.dwnsmpl <- scanone(gg_step2, method = "em", model = "binary", maxit = 1000, n.perm = 10000, pheno.col = 4, n.cluster = 22)
bins <- data.frame(em=summary(scan.bin.em),mr=summary(scan.bin.mr),hk=summary(scan.bin.hk))

####################################################################################
## PERMS WITH ALL LOCI  s
perms.bin.em <- scanone(cross, method = "hk", model = "binary", maxit = 1000,
  n.perm = 10000, pheno.col = 4, n.cluster = 10)
####################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_bin_hk.rsave')))
################################################################################
