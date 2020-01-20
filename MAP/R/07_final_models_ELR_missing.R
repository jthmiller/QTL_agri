#!/bin/R

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

############################################################
pop <- 'ELR.missing'
load(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
############################################################

## manual add stepwise qtl
bin.add.em.perms <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk", n.perm = 10000, n.cluster=10)
lod <- summary(bin.add.em.perms)[2]
bin.add.em <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk")
qtl <- summary(bin.add.em,lod)
bin.add.em.qtls1 <- makeqtl(gg_step2, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
bin.add.em.qtls1 <- refineqtl(gg_step2, qtl=bin.add.em.qtls1, pheno.col=4, model='binary', method = "hk", incl.markers=F)

int.em <- addint(gg_step2, qtl=bin.add.em.qtls1, formula=y~Q1+Q2+Q3, method='hk')
bin.add.em.qtls2 <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls1, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
#int.em <- addint(gg_step2, qtl=bin.add.em.qtls, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')
##scan for an additional additive QTL
add.em.a <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls2, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
##scan for an additional interactive QTL
int.em <- addint(gg_step2, qtl=bin.add.em.qtls2, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')



add.em.i <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q1:Q4)

add.em.2 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q4)
add.em.3 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q3:Q4)

## Scan for interacting pair to add (long)
add.em.c <- addpair(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)





out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)
