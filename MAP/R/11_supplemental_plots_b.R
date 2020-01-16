#!/bin/R
### first run combine pops for multi-pop cross objects

pop <- 'NBH'

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')

mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

#load(file.path(mpath,'08_phys_plots_pos.rsave'))
###############################################################################
## Effect plot (uses sim.geno)
cross_NBH <- sim.geno(cross_NBH, step=1, error.prob=0.05, off.end=5, map.function="kosambi", n.draws=160)
cross_ELR <- sim.geno(cross_ELR, step=1, error.prob=0.05, off.end=5, map.function="kosambi", n.draws=160)


png("/home/jmiller1/public_html/NBH_effect_scan.png", width=1500, height=500)
par(mfrow=c(2,1))
effectscan(cross_NBH, pheno.col=5, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-2,2), main= 'normal')
effectscan(cross_NBH, pheno.col=4, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'binary')
dev.off()


png("/home/jmiller1/public_html/ELR_effect_scan.png", width=1500, height=500)
par(mfrow=c(2,1))
effectscan(cross_ELR, pheno.col=5, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-2,2), main= 'normal')
effectscan(cross_ELR, pheno.col=4, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'binary')
dev.off()



## CHR8 ###################################
################################################################################
elr_ab <- ahr_elr[which(ahr_elr$chr==8),'pos1']
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==8),'pos1']

png("/home/jmiller1/public_html/NBH_ELR_effect_scan_8.png", width=750, height=250)
par(mfrow=c(1,2))
effectscan(cross_NBH, pheno.col=4, chr=8, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'NBH',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
effectscan(cross_ELR, pheno.col=4, chr=8, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'ELR',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=elr_ab,col='red')
dev.off()

################################################################################

png("/home/jmiller1/public_html/NBH_ELR_8_effectplot.png", width=750, height=250)
par(mfrow=c(1,2))
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank, chs=8)
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=8)
dev.off()



################################################################################
################################################################################
col <- c("slateblue", "violetred", "green3")

nbh <- convert2cross2(cross_NBH)
nbh_map <- insert_pseudomarkers(nbh$gmap, step=1)
nbh_pr <- calc_genoprob(nbh, nbh_map, error_prob=0.025, cores=4)

elr <- convert2cross2(cross_ELR)
elr_map <- insert_pseudomarkers(elr$gmap, step=1)
elr_pr <- calc_genoprob(elr, elr_map, error_prob=0.025, cores=4)

cands <- c("AHR1","aip","ARNT","ARNT2","ahrr","ahr1b","AHR2b")

ahr_nbh <- nbh.gens[which(nbh.gens$gene %in% cands),]
ahr_elr <- nbh.gens[which(elr.gens$gene %in% cands),]








plot_ef <- function(crs,map,pr,ahr,popgen){

 for (chr in 1:24){

  c2eff <- scan1coef(pr[,as.character(chr)], crs$pheno[,"pheno_norm"])

  plot(c2eff, map[as.character(chr)], columns=1:3, col=col, ylim=c(0,5), cex.axis = 2)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      abline(v=as.numeric(ahr[indx,'pos1']), col='red')
    }

    if(any( chr %in% popgen$chr )) {
      indx <- which(popgen$chr %in% chr)
      abline(v=as.numeric(popgen[indx,'pos1']), col='red')
    }


  last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients

  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
  }


 dev.off()
}

png("/home/jmiller1/public_html/NBH_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr, ahr = ahr_nbh, popgen = nbh.rank)

png("/home/jmiller1/public_html/ELR_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank)
################################################################################
################################################################################
## Correlate lod and segregation distortion

elr_c2eff <- lapply(1:24,function(X) {
 scan1coef(elr_pr[,as.character(X)], elr$pheno[,"bin"])
})
elr_c2eff <- do.call(rbind,elr_c2eff)
elr_seg <- geno.table(cross_ELR)[rownames(elr_c2eff),'P.value']

nbh_c2eff <- lapply(1:24,function(X) {
 scan1coef(nbh_pr[,as.character(X)], nbh$pheno[,"bin"])
})
nbh_c2eff <- do.call(rbind,nbh_c2eff)
nbh_seg <- geno.table(cross_NBH)[rownames(nbh_c2eff),'P.value']



png("/home/jmiller1/public_html/lodxdist.png", width=500, height=500)
plot(-log10(nbh_seg),nbh_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(nbh_seg),nbh_c2eff[,'BB'], col = 'blue',pch=19)
points(-log10(nbh_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)

dev.off()

png("/home/jmiller1/public_html/lodxdist_nbh_AABB.png", width=500, height=500)
plot(-log10(nbh_seg),nbh_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(nbh_seg),nbh_c2eff[,'BB'], col = 'red',pch=19)
##points(-log10(nbh_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)
dev.off()


png("/home/jmiller1/public_html/lodxdist_elr_AABB.png", width=500, height=500)
plot(-log10(elr_seg),elr_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(elr_seg),elr_c2eff[,'BB'], col = 'red',pch=19)
##points(-log10(elr_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)
dev.off()


################################################################################
################################################################################

plot_pgen <- function(crs,stat, map, ahr, ahr_clm, colnm, popgen, ylimo,rank_clm){

 for (chr in 1:24){

  xl <- summary(pull.map(crs))[chr,'length']
  ind <- which(stat$chr == chr)

  Y <- stat[ind,colnm]
  X <- stat[ind,map]
##  plot(X, Y, col='blue', cex.axis = 2, ylim = ylimo, xlim = c(0,xl), main=paste('CHR',chr), cex.main=2)

  plot(X, Y, col='blue', cex.axis = 2, ylim = ylimo, main=paste('CHR',chr), cex.main=2)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      abline(v=as.numeric(ahr[indx,ahr_clm]), col='red')
    }

    if(any( chr %in% popgen$chr )) {
      indx <- which(popgen$chr %in% chr)
      abline(v=as.numeric(popgen[indx,rank_clm]), col='red')
    }
 }
dev.off()
}


##PHYS##############################
png("/home/jmiller1/public_html/nbh_pbs_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pbs, map = 'mid' , ahr = ahr_nbh, ahr_clm= 'stp',  colnm = 'NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(-0.25,1.1) )

png("/home/jmiller1/public_html/nbh_pfst_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pfst, map = 'mid', ahr = ahr_nbh, ahr_clm= 'stp', colnm = 'BI.NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(0,1) )

png("/home/jmiller1/public_html/elr_pbs_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pbs, map = 'mid' , ahr = ahr_elr,ahr_clm= 'stp',  colnm = 'ER', popgen = elr.rank, rank_clm='end', ylimo=c(-1,2) )

png("/home/jmiller1/public_html/elr.sh_pfst_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'mid', ahr = ahr_elr, ahr_clm= 'stp', colnm = 'ER.SH', popgen = elr.rank, rank_clm='end', ylimo=c(0,1.5) )

png("/home/jmiller1/public_html/elr.kc_pfst_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'mid', ahr = ahr_elr, ahr_clm= 'stp', colnm = 'ER.KC', popgen = elr.rank, rank_clm='end', ylimo=c(0,1.5) )


####################################

##CM POS############################

png("/home/jmiller1/public_html/nbh_pfst.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pfst, map = 'nbh_cm' , ahr = ahr_nbh, ahr_clm= 'pos1',  colnm = 'BI.NBH' , popgen = nbh.rank, ylimo=c(0,0.75) )

png("/home/jmiller1/public_html/nbh_pbs.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pbs, map = 'nbh_cm' , ahr = ahr_nbh,  colnm = 'NBH' , popgen = nbh.rank, ylimo=c(-1,2) )

png("/home/jmiller1/public_html/elr_pfst.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'elr_cm', ahr = ahr_elr, ahr_clm= 'stp',  colnm = 'ER.KC' , popgen = elr.rank,rank_clm='end', ylimo=c(0,1) )

png("/home/jmiller1/public_html/elr_pfst.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'elr_cm', ahr = ahr_elr, ahr_clm= 'stp',  colnm = 'ER.KC' , popgen = elr.rank,rank_clm='end', ylimo=c(0,1) )


get_genes_cm(chr=2, start = 32809365,stop = 32962365, models = nbh_gene_models, colm = 'start')

######

a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_ELR)[[X]], 0.25)} ))
a <- pull.markers(cross_ELR,a)
a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

png("/home/jmiller1/public_html/elr_rf.png", width=1000, height=1000)
 plotRF(rf, zmax=8, col.scheme="redblue",what='lod')
dev.off()



a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_NBH)[[X]], 0.25)} ))
a <- pull.markers(cross_NBH,a)
#a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

png("/home/jmiller1/public_html/nbh_rf.png", width=1000, height=1000)
 plotRF(rf, zmax=10, col.scheme="redblue",what='lod')
dev.off()
