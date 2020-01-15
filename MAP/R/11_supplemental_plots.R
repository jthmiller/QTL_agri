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

load(file.path(mpath,'08_phys_plots_pos.rsave'))
######## Plot phys pos x marker order ##########################################

png("/home/jmiller1/public_html/ELR_NBH_physpo_filt.png", width=1500, height=1500)
par(mfrow=c(4,6))

for (i in 1:24){

 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross_NBH,i))))
 X <- 1:length(Y)

 A <- c(0, as.numeric(gsub(".*:","",markernames(cross_ELR,i))))
 B <- 1:length(A)

 ymax <- max(c(max(Y),max(A))
 xmax <- max(c(length(Y),length(A))

 plot(c(0,xmax),c(0,ymax), type="n", xlab=paste('chr',i), ylab='physical position')
 points(X,Y, col='blue')
 points(A,B, col='yellow')

 }
dev.off()
################################################################################

######## Plot phys pos x marker order ##########################################

png("/home/jmiller1/public_html/NBH_physpo_filt.png", width=1500, height=1500)
par(mfrow=c(4,6))

for (i in 1:24){

 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross,i))))
 X <- 1:length(Y)

 plot(X,Y, xlab=paste('chr',i), ylab='physical position',cex.axis=2,pch=21,cex.main=2)


 }
dev.off()
################################################################################

#### SEG DISTORTION X PHYSICAL POSITION
gt <- geno.table(cross_nbh)

png("/home/jmiller1/public_html/NBH_segdist.png", width=1500, height=1500)
par(mfrow=c(4,6))

for (i in 1:24){

 ind <- which(gt$chr == i)
 X <- c(as.numeric(gsub(".*:","",rownames(gt)[ind])))
 Y <- c(-log10(gt[ind,'P.value']))

 plot(X,Y,ylim=c(0,4),xlim=c(0,max(X)), xlab=paste('chr',i), ylab='-log10 pval', cex.axis=2,pch=21,cex.main=2)


 }
dev.off()
################################################################################

fl <- paste0('ELR.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cross_ELR <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

gt.elr <- geno.table(cross_ELR)

fl <- paste0('NBH.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cross_NBH <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

gt.nbh <- geno.table(cross_NBH)



#### ELR and NBH seg dist
png("/home/jmiller1/public_html/NBH_ELR_segdist.png", width=1750, height=1750)
par(mfrow=c(4,6))

for (i in 1:24){


 ind.nbh <- which(gt.nbh$chr == i)
 X <- c(as.numeric(gsub(".*:","",rownames(gt.nbh)[ind.nbh])))
 Y <- c(-log10(gt.nbh[ind.nbh,'P.value']))


 ind.elr <- which(gt.elr$chr == i)
 A <- c(as.numeric(gsub(".*:","",rownames(gt.elr)[ind.elr])))
 B <- c(-log10(gt.elr[ind.elr,'P.value']))

 ymax <- max(c(max(Y),max(B)))
 xmax <- max(c(max(X),max(A)))

 plot(c(0,xmax),c(0,ymax),ylim=c(0,4), type="n", xlab=paste('chr',i,'phys position'), ylab='-log10 p.value',cex.lab=2,cex.sub=2, cex.axis=2,cex.main=2)
 points(X,Y, col='blue',pch=19)
 points(A,B, col='red',pch=19)

 }
dev.off()


################################################################################
cross_NBH$pheno$pheno_norm <- round(nqrank(cross_NBH$pheno$Pheno),5)
cross_ELR$pheno$pheno_norm <- round(nqrank(cross_ELR$pheno$Pheno),5)

cross_NBH <- sim.geno(cross_NBH, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=160)
cross_ELR <- sim.geno(cross_ELR, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=160)
cross_NBH <- calc.genoprob(cross_NBH, step=1, error.prob=0.025, off.end=5, map.function="kosambi")
cross_ELR <- calc.genoprob(cross_ELR, step=1, error.prob=0.025, off.end=5, map.function="kosambi")
scan_nbh <- scanone(cross_NBH, method = "mr", model = "binary", pheno.col = 4)
scan_elr <- scanone(cross_ELR, method = "mr", model = "binary", pheno.col = 4)
scan_nbh <- scanone(cross_NBH, method = "imp", model = "binary", pheno.col = 4)
scan_elr <- scanone(cross_ELR, method = "imp", model = "binary", pheno.col = 4)



png("/home/jmiller1/public_html/NBH_ELR_segdist_binary.png", width=1750, height=1750)
par(mfrow=c(4,6))

for (i in 1:24){

 ind.nbh <- which(gt.nbh$chr == i)
 Ylod <- scan_nbh$lod[ind.nbh]
 X <- c(as.numeric(gsub(".*:","",rownames(gt.nbh)[ind.nbh])))
 Y <- c(-log10(gt.nbh[ind.nbh,'P.value']))

 ind.elr <- which(gt.elr$chr == i)
 Blod <- scan_elr$lod[ind.elr]
 A <- c(as.numeric(gsub(".*:","",rownames(gt.elr)[ind.elr])))
 B <- c(-log10(gt.elr[ind.elr,'P.value']))

 ymax <- max(c(max(Y),max(B)))
 xmax <- max(c(max(X),max(A)))


 if(max(Ylod) > 5 ) Ylod <- rescale(scan_nbh$lod[ind.nbh],c(0,5))
 if(max(Blod) > 5 ) Blod <- rescale(scan_elr$lod[ind.elr],c(0,5))

 plot(c(0,xmax),c(0,ymax),ylim=c(0,5), type="n", xlab=paste('chr',i,'phys position'), ylab='-log10 p.value',cex.lab=2,cex.sub=2, cex.axis=2,cex.main=2)
 points(X,Ylod, col='cornflowerblue',pch=21)
 points(A,Blod, col='pink',pch=21)
 points(X,Y, col='blue',pch=19)
 points(A,B, col='red',pch=19)

 }
dev.off()





################################################################################
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

################################################################################

iron <- convert2cross2(cross_ELR)
map <- insert_pseudomarkers(iron$gmap, step=1)
pr <- calc_genoprob(iron, map, error_prob=0.025, cores=4)
apr <- genoprob_to_alleleprob(pr)

chr <- '1'

c2eff <- scan1coef(pr[,as.character(chr)], iron$pheno[,"bin"])

png("/home/jmiller1/public_html/ELR_effect_1.png", width=500, height=500)

par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map[as.character(chr)], columns=1:3, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

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


plot_ef <- function(crs,map,pr){

 for (chr in 1:24){

  c2eff <- scan1coef(pr[,as.character(chr)], crs$pheno[,"bin"])

  plot(c2eff, map[as.character(chr)], columns=1:3, col=col, ylim=c(0,1), cex.axis = 2)
  last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients

  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

  }
 dev.off()
}

png("/home/jmiller1/public_html/NBH_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr )

png("/home/jmiller1/public_html/ELR_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = elr, map = elr_map, pr = elr_pr )
################################################################################
################################################################################
