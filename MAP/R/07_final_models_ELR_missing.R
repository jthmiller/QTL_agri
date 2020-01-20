#!/bin/R

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
pop <- 'ELR.missing'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

############################################################
load(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
############################################################

full.norm.imp <- stepwiseqtl(cross, penalties=pens, incl.markers=T, qtl=norm.add.qtls, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8, chr=qtls_chr)


## manual add stepwise qtl
bin.add.hk.perms <- scanone(cross, pheno.col=4, model='binary', method = "hk", n.perm = 10000, n.cluster=10)
lod <- summary(bin.add.hk.perms)[2]
bin.add.hk <- scanone(cross, pheno.col=4, model='binary', method = "hk")
qtl <- summary(bin.add.hk,lod)
bin.add.hk.qtls1 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
bin.add.hk.qtls1 <- refineqtl(cross, qtl=bin.add.hk.qtls1, pheno.col=4, model='binary', method = "hk", incl.markers=F)

int.hk <- addint(cross, qtl=bin.add.hk.qtls1, formula=y~Q1+Q2+Q3, method='hk')
bin.add.hk.qtls2 <- refineqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls1, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
#int.hk <- addint(cross, qtl=bin.add.hk.qtls, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')
##scan for an additional additive QTL
add.hk.a <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls2, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
##scan for an additional interactive QTL
int.hk <- addint(cross, qtl=bin.add.hk.qtls2, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')



add.hk.i <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q1:Q4)

add.hk.2 <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q4)
add.hk.3 <- addqtl(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q3:Q4)

## Scan for interacting pair to add (long)
add.hk.c <- addpair(cross, pheno.col=4, model='binary', qtl=bin.add.hk.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)





out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)
