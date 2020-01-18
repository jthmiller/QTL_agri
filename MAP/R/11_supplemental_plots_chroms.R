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
cross_NBH <- jittermap(cross_NBH)
gt <- geno.table(cross_NBH)
gt.map <- pull.map(cross_NBH)

cross_grid <- reduce2grid(cross_NBH)
scan_nbh <- scanone(cross_grid, method = "em", model = "binary", pheno.col = 4)
mpl <- summary(pull.map(cross_grid))

gg_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_NBH)[[X]], 1)} ))
gg <- pull.markers(cross_NBH,gg_marks)

cross_NBH$pheno$pheno_norm <- round(nqrank(cross_NBH$pheno$Pheno),5)
norm_nbh <- scanone(cross_grid, method = "imp", model = "normal", pheno.col = 5)
bin_nbh <- scanone(cross_grid, method = "em", model = "binary", pheno.col = 4)

col <- c("slateblue", "violetred", "green3")

nbh <- convert2cross2(cross_NBH)
nbh_map <- insert_pseudomarkers(nbh$gmap, step=1)
nbh_pr <- calc_genoprob(nbh, nbh_map, error_prob=0.0025, cores=4)
nbh_bin <- scan1(nbh_pr, pheno=nbh$pheno[,'bin'] , model = "binary", cores = 8)

################################################################################

#bottom, left, top and right
#location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks

ch <- 8
len <- mpl[ch,'length']
ind <- which(gt$chr == ch)
segX <- gt.map[[ch]][rownames(gt)[ind]]
#segX <- c(as.numeric(gsub(".*:","",rownames(gt)[ind])))/1000000
segY <- c(-log10(gt[ind,'P.value']))
nbh_ab <- ahr_nbh[which(ahr_nbh$chr==ch),'pos1']

pdf(paste0("/home/jmiller1/public_html/NBH_all",ch,"segdist.pdf"), width=4.5,height=7)
mat<-matrix(c(1:7),7,1, byrow=T)

layout(mat, widths=1, heights= c(0.1, 0.1, 0.025, 0.15, 0.15, 0.1, 0.1))
par(mar=c(0.5,2,0,1)+0.1, oma = c(1, 1, 0, 1))

plot_pgen(crs = cross_NBH, chrs=ch, stat = pbs, map = 'mid' ,
 ahr = ahr_nbh, ahr_clm= 'stp',  colnm = 'NBH', popgen = nbh.rank,
 rank_clm='end', ylimo=c(-0.08,0.5), stat_name='pbs', pch=16, mgp = c(1.25, 0.5, 0),
 cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n")

plot_pgen(crs = cross_NBH, chrs=ch, stat = taj, map = 'mid',
 ahr = ahr_nbh, ahr_clm= 'stp', colnm = 'NBH', popgen = nbh.rank,
 rank_clm = 'end', ylimo=c(-4,4), stat_name='taj', pch=16, mgp = c(1.25, 0.5, 0),
 cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75)

#plot(0,type='n',cex=0.25, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n")
plot(0,type='n',axes=FALSE,ann=FALSE)
## CM
plot(nbh_bin, nbh_map,chr=ch, lodcolumn=1, ylim=c(0,10), xaxt="n", mgp = c(1.25, 0.5, 0),)
abline(v=nbh_ab, col='red')

plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr_nbh,
 popgen = nbh.rank, chs=ch, main=NULL, model="bin", xaxt="n",
 cex.axis = 0.5, cex.lab= 0.75, cex.main = 0.75)

effectscan(cross_grid, pheno.col=4, chr=ch, get.se=T, draw=TRUE,
 gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5),
 xlim=c(0,len), main=NULL,ylab='Effect Est. (+/- 1 SE)', xaxs="i",, xaxt="n",
 cex.axis = 0.75, cex.lab= 0.75, cex.main = 0.75, mgp = c(1.5, 0.5, 0),)
abline(v=nbh_ab,col='red',lwd=0.5)

plot(segX,segY,ylim=c(0,4),xlim=c(0,max(segX)), xlab="",
 ylab='-log10 pval', pch=16, xaxs="i", mgp = c(1.25, 0.5, 0),
 cex=0.25, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75)
abline(v=nbh_ab,col='red',lwd=0.5)
dev.off()

################################################################################
