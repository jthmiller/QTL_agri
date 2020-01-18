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

load(file.path(mpath,'supplemental_plot_env.rsave'))

erp <- 0.0025

cross_NBH <- sim.geno(cross_NBH, step=1, error.prob=erp, off.end=5, map.function="kosambi", n.draws=100, stepwidth="fixed")
cross_NBH <- calc.genoprob(cross_NBH, step=1, error.prob=erp, off.end=5, map.function="kosambi", stepwidth="fixed")
cross_grid <- reduce2grid(cross_NBH)
scan_nbh <- scanone(cross_grid, method = "em", model = "binary", pheno.col = 4)

cross_NBH$pheno$pheno_norm <- round(nqrank(cross_NBH$pheno$Pheno),5)
norm_nbh <- scanone(cross_grid, method = "imp", model = "normal", pheno.col = 5)
bin_nbh <- scanone(cross_grid, method = "em", model = "binary", pheno.col = 4)

col <- c("slateblue", "violetred", "green3")

nbh <- convert2cross2(cross_NBH)
nbh_map <- insert_pseudomarkers(nbh$gmap, step=1)
nbh_pr <- calc_genoprob(nbh, nbh_map, error_prob=0.0025, cores=4)
nbh_bin <- scan1(nbh_pr, pheno=nbh$pheno[,'bin'] , model = "binary", cores = 8)

nbh_ab <- ahr_nbh[which(ahr_nbh$chr==ch),'pos1']

gt <- geno.table(cross_NBH)
ind <- which(gt$chr == ch)
segX <- c(as.numeric(gsub(".*:","",rownames(gt)[ind])))
segY <- c(-log10(gt[ind,'P.value']))


##plot.scanone, use incl.markers=FALSE

## plot:
PBS/fst
LOD profile
QTL effect
Model support
seg distortion


ch <- 2

pdf(paste0("/home/jmiller1/public_html/NBH_",ch,"segdist.pdf"), width=6, height=3.5)
plot(segX,segY,ylim=c(0,4),xlim=c(0,max(segX)), xlab=paste('chr',ch), ylab='-log10 pval', cex.axis=2,pch=18,cex.main=2)
dev.off()

pdf(paste0("/home/jmiller1/public_html/nbh_pbs_phys",ch,".pdf"), width=6, height=3.5)
plot_pgen(crs = cross_NBH, chrs=ch, stat = pbs, map = 'mid' , ahr = ahr_nbh, ahr_clm= 'stp',  colnm = 'NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(-0.25,1.1),stat_name='pbs' )
dev.off()

pdf(paste0("/home/jmiller1/public_html/nbh_pfst_phys",ch,".pdf"), width=6, height=3.5)
plot_pgen(crs = cross_NBH, chrs=ch, stat = pfst, map = 'mid' , ahr = ahr_nbh, ahr_clm= 'stp',  colnm = 'BI.NBH', popgen = nbh.rank, rank_clm='end', ylimo=c(-0.25,1.1),stat_name='pbs' )
dev.off()

##pdf(paste0("/home/jmiller1/public_html/nbh_scanone",ch,".pdf"), width=6, height=3.5)
##plot(bin_nbh,chr=2,incl.markers=FALSE)
##dev.off()

pdf(paste0("/home/jmiller1/public_html/nbh_scan1",ch,".pdf"), width=6, height=3.5)
plot(nbh_bin, nbh_map,chr=2, lodcolumn=1)
dev.off()

pdf("/home/jmiller1/public_html/NBH_2_effectplot.pdf", width=6, height=3.5)
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh, popgen = nbh.rank, chs=ch, main=NULL,model="bin")
dev.off()

pdf("/home/jmiller1/public_html/NBH_effect_scan_2.pdf", width=6, height=3.5)
effectscan(cross_grid, pheno.col=4, chr=ch, get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main=NULL,ylab='Effect Est. (+/- 1 SE)',cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=nbh_ab,col='red')
dev.off()
