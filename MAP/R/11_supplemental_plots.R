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
21610008 21865786
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
for(ch in 1:24){
svg(paste0("/home/jmiller1/public_html/nbh_pbs_phys",ch,".svg"), width=6, height=3.5)
plot_pgen(crs = cross_NBH, chrs=ch, stat = pbs, map = 'mid' , ahr = ahr_nbh, ahr_clm= 'stp',  colnm = 'NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(-0.25,1.1),stat_name='pbs' )
dev.off()

svg(paste0("/home/jmiller1/public_html/nbh_pfst_phys",ch,".svg"), width=6, height=3.5)
plot_pgen(crs = cross_NBH,  chrs=ch, stat = pfst, map = 'mid', ahr = ahr_nbh, ahr_clm= 'stp', colnm = 'BI.NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(0,1),stat_name='pfst' )
dev.off()
}

##PHYS##########################################################################
pdf("/home/jmiller1/public_html/nbh_pbs_phys.pdf", width=35, height=10)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pbs, map = 'mid' , ahr = ahr_nbh, ahr_clm= 'stp',  colnm = 'NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(-0.25,1.1),stat_name='pbs' )
dev.off()

pdf("/home/jmiller1/public_html/nbh_pfst_phys.pdf", width=35, height=10)
par(mfrow=c(4,6))
plot_pgen(crs = cross_NBH, stat = pfst, map = 'mid', ahr = ahr_nbh, ahr_clm= 'stp', colnm = 'BI.NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(0,1),stat_name='pfst' )

pdf("/home/jmiller1/public_html/elr_pbs_phys.pdf", width=35, height=10)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pbs, map = 'mid' , ahr = ahr_elr,ahr_clm= 'stp',  colnm = 'ER', popgen = elr.rank, rank_clm='end', ylimo=c(-1,2),stat_name='pbs' )

pdf("/home/jmiller1/public_html/elr.sh_pfst_phys.pdf", width=35, height=10)
par(mfrow=c(4,6))
plot_pgen(crs = cross_ELR, stat = pfst, map = 'mid', ahr = ahr_elr, ahr_clm= 'stp', colnm = 'ER.SH', popgen = elr.rank, rank_clm='end', ylimo=c(0,1.5),stat_name='pfst' )

pdf("/home/jmiller1/public_html/elr.kc_pfst_phys.pdf", width=35, height=10)
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

pdf("/home/jmiller1/public_html/elr_rf.pdf", width=1000, height=1000)
 plotRF(rf, zmax=8, col.scheme="redblue",what='lod')
dev.off()

################################################################################
a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_NBH)[[X]], 0.25)} ))
a <- pull.markers(cross_NBH,a)
#a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

pdf("/home/jmiller1/public_html/nbh_rf.pdf", width=1000, height=1000)
 plotRF(rf, zmax=10, col.scheme="redblue",what='lod')
dev.off()
################################################################################
################################################################################

################################################################################
################################################################################
### CANDIDATE GENES ######
################################################################################
################################################################################

## get_genes_cm(chr=2, start = 32809365,stop = 32962365, models = nbh_gene_models, colm = 'start')
## get_genes_cm(chr=1, start = 20,stop = 30,models = nbh_gene_models, colm = 'cm_mid')
chr24_195 <- get_genes_cm(chr=24, start = 37087357,stop = 37192357,models = nbh_gene_models, colm = 'start')
chr24_196 <- get_genes_cm(chr=24, start = 38000000,stop = 39000000,models = nbh_gene_models, colm = 'start')
chr24_198 <- get_genes_cm(chr=24, start = 34684831,stop = 34754831,models = nbh_gene_models, colm = 'start')

chr13_293 <- get_genes_cm(chr=13, start = 7039628,stop = 7101628,models = nbh_gene_models, colm = 'start')

chr13_fst_out <- get_genes_cm(chr=13, start = 23896922, stop = 24223422,models = nbh_gene_models, colm = 'start')

chr13_AB_QTL <- get_genes_cm(chr=13, start = 35, stop = 45,models = nbh_gene_models, colm = 'cm_mid')
chr13_AB_QTL <- chr13_AB_QTL[grep('LOC',chr13_AB_QTL$name,invert=T),]

thrt <- pfst[which(pfst$chr == 13),]
thrt <- thrt[!is.na(thrt$ER.KC),]
tail(thrt[order(thrt$ER.KC),])
37087357 37192357
t[grep('LOC',t$name,invert=T),]


sb_pbs <- pbs[which(pbs$chr == 13),]
sb_pbs <- sb_pbs[!is.na(sb_pbs$ER),]
sb_pbs <- tail(sb_pbs[order(sb_pbs$ER),],30)
sb_pbs <- sb_pbs[order(sb_pbs$start),]
###       chr    start      end        NBH          BP          NYC        ER
### 790629  13  8173628  8178628 0.14751252 0.016933220 -0.022930839 1.6930752
chr13_elr_pbs_peak <- get_genes_cm(chr=13, start = 7973628, stop = 8178628,models = nbh_gene_models, colm = 'start')



sb_pfst <- pfst[which(pfst$chr == 13),]
sb_pfst <- sb_pfst[!is.na(sb_pfst$ER.S),]
tail(sb_pfst[order(sb_pfst$ER.SH),],30)

chr13_elr_pfst_peak <- get_genes_cm(chr=13, start = 7873628, stop = 8378628,models = nbh_gene_models, colm = 'start')

c('16419','939044','939045')


chr13_elr_pfst_peak2 <- get_genes_cm(chr=13, start = 23891922, stop = 25728077, models = nbh_gene_models, colm = 'start')

pfst['790629',]

thrt <- pbs[which(pbs$chr == 13),]
thrt <- thrt[!is.na(thrt$ER),]
tail(thrt[order(thrt$ER),])



sb_pbs <- pbs[which(pbs$chr == 8),]
sb_pbs <- sb_pbs[!is.na(sb_pbs$ER),]
sb_pbs <- tail(sb_pbs[order(sb_pbs$ER),],60)
sb_pbs <- sb_pbs[order(sb_pbs$start),]

chr8_elr_pbs_peak <- get_genes_cm(chr=8, start = 16449074, stop = 18239074,models = nbh_gene_models, colm = 'start')

################################################################################

chr1_nbh_rank_incompat_peak <- get_genes_cm(chr=1, start = 21600000, stop =  21800000,models = nbh_gene_models, colm = 'start')
## kcnj3 regulates heartbeat

################################################################################
################################################################################
