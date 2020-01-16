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

pdf("/home/jmiller1/public_html/NBH_effect_scan.pdf", width=1500, height=500)
par(mfrow=c(2,1))
effectscan(cross_NBH, pheno.col=5, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-2,2), main= 'normal')
effectscan(cross_NBH, pheno.col=4, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'binary')
dev.off()


pdf("/home/jmiller1/public_html/ELR_effect_scan.pdf", width=1500, height=500)
par(mfrow=c(2,1))
effectscan(cross_ELR, pheno.col=5, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-2,2), main= 'normal')
effectscan(cross_ELR, pheno.col=4, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'binary')
dev.off()

pdf("/home/jmiller1/public_html/NBH_effectplot.pdf", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr, ahr = ahr_nbh, popgen = nbh.rank)

pdf("/home/jmiller1/public_html/ELR_effectplot.pdf", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank)
## CHR8 ###################################
################################################################################
elr_ab <- ahr_elr[which(ahr_elr$chr==8),'pos1']
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==8),'pos1']

pdf("/home/jmiller1/public_html/NBH_ELR_effect_scan_8.pdf", width=750, height=250)
par(mfrow=c(1,2))
effectscan(cross_NBH, pheno.col=4, chr=8, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'NBH',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
effectscan(cross_ELR, pheno.col=4, chr=8, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'ELR',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=elr_ab,col='red')
dev.off()

pdf("/home/jmiller1/public_html/NBH_ELR_8_effectplot.pdf", width=750, height=250)
par(mfrow=c(1,2))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=8, main='NBH',model="bin")
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank, chs=8, main='ELR',model="bin")
dev.off()
################################################################################

## CHR13 ###################################
################################################################################
elr_ab <- ahr_elr[which(ahr_elr$chr==13),'pos1']
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==13),'pos1']

pdf("/home/jmiller1/public_html/NBH_ELR_effect_scan_13.pdf", width=750, height=250)
par(mfrow=c(1,2))
effectscan(cross_NBH, pheno.col=4, chr=13, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'NBH',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
effectscan(cross_ELR, pheno.col=4, chr=13, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'ELR',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=elr_ab,col='red')
dev.off()

pdf("/home/jmiller1/public_html/NBH_ELR_13_effectplot.pdf", width=750, height=250)
par(mfrow=c(1,2))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=13, main='NBH',model="bin")
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank, chs=13, main='ELR',model="bin")
dev.off()
################################################################################
################################################################################

## CHR2 ###################################
################################################################################
elr_ab <- ahr_elr[which(ahr_elr$chr==2),'pos1']
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==2),'pos1']

pdf("/home/jmiller1/public_html/NBH_ELR_effect_scan_2.pdf", width=750, height=250)
par(mfrow=c(1,2))
effectscan(cross_NBH, pheno.col=4, chr=2, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'NBH',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
effectscan(cross_ELR, pheno.col=4, chr=2, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'ELR',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=elr_ab,col='red')
dev.off()

pdf("/home/jmiller1/public_html/NBH_ELR_2_effectplot.pdf", width=750, height=250)
par(mfrow=c(1,2))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=2, main='NBH',model="bin")
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank, chs=2, main='ELR',model="bin")
dev.off()
################################################################################
################################################################################

## CHR18 ###################################
################################################################################
elr_ab <- ahr_elr[which(ahr_elr$chr==18),'pos1']
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==18),'pos1']

pdf("/home/jmiller1/public_html/NBH_ELR_effect_scan_18.pdf", width=750, height=250)
par(mfrow=c(1,2))
effectscan(cross_NBH, pheno.col=4, chr=18, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'NBH',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
effectscan(cross_ELR, pheno.col=4, chr=18, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'ELR',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=elr_ab,col='red')
dev.off()

pdf("/home/jmiller1/public_html/NBH_ELR_18_effectplot.pdf", width=750, height=250)
par(mfrow=c(1,2))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=18, main='NBH',model="bin")
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank, chs=18, main='ELR',model="bin")
dev.off()
################################################################################
################################################################################

## CHR24 ###################################
################################################################################
elr_ab <- ahr_elr[which(ahr_elr$chr==24),'pos1']
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==24),'pos1']

pdf("/home/jmiller1/public_html/NBH_ELR_effect_scan_24.pdf", width=750, height=250)
par(mfrow=c(1,2))
effectscan(cross_NBH, pheno.col=4, chr=24, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'NBH',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
effectscan(cross_ELR, pheno.col=4, chr=24, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'ELR',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=elr_ab,col='red')
dev.off()

pdf("/home/jmiller1/public_html/NBH_ELR_24_effectplot.pdf", width=750, height=250)
par(mfrow=c(1,2))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=24, main='NBH',model="bin")
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank, chs=24, main='ELR',model="bin")
dev.off()
################################################################################


## CHR1 ###################################
################################################################################
elr_ab <- ahr_elr[which(ahr_elr$chr==1),'pos1']
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==1),'pos1']

pdf("/home/jmiller1/public_html/NBH_ELR_effect_scan_1.pdf", width=750, height=250)
par(mfrow=c(1,2))
effectscan(cross_NBH, pheno.col=4, chr=1, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'NBH',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
effectscan(cross_ELR, pheno.col=4, chr=1, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'ELR',ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=elr_ab,col='red')
dev.off()

pdf("/home/jmiller1/public_html/NBH_ELR_1_effectplot.pdf", width=750, height=250)
par(mfrow=c(1,2))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=1, main='NBH',model="bin")
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank, chs=1, main='ELR',model="bin")
dev.off()
################################################################################





################################################################################
################################################################################
## Correlate lod and segregation distortion

pdf("/home/jmiller1/public_html/lodxdist.pdf", width=500, height=500)
plot(-log10(nbh_seg),nbh_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(nbh_seg),nbh_c2eff[,'BB'], col = 'blue',pch=19)
points(-log10(nbh_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)

dev.off()

pdf("/home/jmiller1/public_html/lodxdist_nbh_AABB.pdf", width=500, height=500)
plot(-log10(nbh_seg),nbh_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(nbh_seg),nbh_c2eff[,'BB'], col = 'red',pch=19)
##points(-log10(nbh_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)
dev.off()


pdf("/home/jmiller1/public_html/lodxdist_elr_AABB.pdf", width=500, height=500)
plot(-log10(elr_seg),elr_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(elr_seg),elr_c2eff[,'BB'], col = 'red',pch=19)
##points(-log10(elr_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)
dev.off()


################################################################################
################################################################################



##PHYS##########################################################################
png("/home/jmiller1/public_html/nbh_pbs_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pbs, map = 'mid' , ahr = ahr_nbh, ahr_clm= 'stp',  colnm = 'NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(-0.25,1.1),stat_name='pbs' )

png("/home/jmiller1/public_html/nbh_pfst_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pfst, map = 'mid', ahr = ahr_nbh, ahr_clm= 'stp', colnm = 'BI.NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(0,1),stat_name='pfst' )

png("/home/jmiller1/public_html/elr_pbs_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pbs, map = 'mid' , ahr = ahr_elr,ahr_clm= 'stp',  colnm = 'ER', popgen = elr.rank, rank_clm='end', ylimo=c(-1,2),stat_name='pbs' )

png("/home/jmiller1/public_html/elr.sh_pfst_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'mid', ahr = ahr_elr, ahr_clm= 'stp', colnm = 'ER.SH', popgen = elr.rank, rank_clm='end', ylimo=c(0,1.5),stat_name='pfst' )

png("/home/jmiller1/public_html/elr.kc_pfst_phys.png", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'mid', ahr = ahr_elr, ahr_clm= 'stp', colnm = 'ER.KC', popgen = elr.rank, rank_clm='end', ylimo=c(0,1.5),stat_name='pfst' )
#####################################################################################


##CM POS################################################################

pdf("/home/jmiller1/public_html/nbh_pfst.pdf", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pfst, map = 'nbh_cm' , ahr = ahr_nbh, ahr_clm= 'pos1',  colnm = 'BI.NBH' , popgen = nbh.rank, ylimo=c(0,0.75) )

pdf("/home/jmiller1/public_html/nbh_pbs.pdf", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pbs, map = 'nbh_cm' , ahr = ahr_nbh,  colnm = 'NBH' , popgen = nbh.rank, ylimo=c(-1,2) )

pdf("/home/jmiller1/public_html/elr_pfst.pdf", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'elr_cm', ahr = ahr_elr, ahr_clm= 'stp',  colnm = 'ER.KC' , popgen = elr.rank,rank_clm='end', ylimo=c(0,1) )

pdf("/home/jmiller1/public_html/elr_pfst.pdf", width=2500, height=750)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'elr_cm', ahr = ahr_elr, ahr_clm= 'stp',  colnm = 'ER.KC' , popgen = elr.rank,rank_clm='end', ylimo=c(0,1) )

##get_genes_cm(chr=2, start = 32809365,stop = 32962365, models = nbh_gene_models, colm = 'start')

################################################################################
## rf interaction plot #########################################################
################################################################################
a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_ELR)[[X]], 0.25)} ))
a <- pull.markers(cross_ELR,a)
a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

png("/home/jmiller1/public_html/elr_rf.png", width=1000, height=1000)
 plotRF(rf, zmax=8, col.scheme="redblue",what='lod')
dev.off()

################################################################################
a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_NBH)[[X]], 0.25)} ))
a <- pull.markers(cross_NBH,a)
#a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

png("/home/jmiller1/public_html/nbh_rf.png", width=1000, height=1000)
 plotRF(rf, zmax=10, col.scheme="redblue",what='lod')
dev.off()
################################################################################
################################################################################
