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




plot_ef <- function(crs,map,pr,ahr,popgen,chs,main,model=c("bin","pheno_norm")){

 for (chr in chs){

  c2eff <- scan1coef(pr[,as.character(chr)], crs$pheno[,model])

  plot(c2eff, map[as.character(chr)], columns=1:3, col=col, ylim=c(0,1), cex.axis = 2,main=main)

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

}

png("/home/jmiller1/public_html/NBH_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr, ahr = ahr_nbh, popgen = nbh.rank,chs=1:24)
dev.off()

png("/home/jmiller1/public_html/ELR_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank,chs=1:24)
dev.off()
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

plot_pgen <- function(crs,stat, map, ahr, ahr_clm, colnm, popgen, ylimo,rank_clm,stat_name){

 for (chr in 1:24){

  xl <- summary(pull.map(crs))[chr,'length']
  ind <- which(stat$chr == chr)

  Y <- stat[ind,colnm]
  X <- stat[ind,map]
##  plot(X, Y, col='blue', cex.axis = 2, ylim = ylimo, xlim = c(0,xl), main=paste('CHR',chr), cex.main=2)

  plot(X, Y, col='blue', cex.axis = 2, ylim = ylimo, main=paste('CHR',chr), cex.lab=2, cex.main=2, xlab='physical position', ylab=stat_name)

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
