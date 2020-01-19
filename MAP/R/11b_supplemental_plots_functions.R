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

######## Plot phys pos x marker order ##########################################

pdf("/home/jmiller1/public_html/ELR_NBH_physpo_filt.pdf", width=21, height=21)
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

pdf("/home/jmiller1/public_html/NBH_physpo_filt.pdf", width=21, height=21)
par(mfrow=c(4,6))

for (i in 1:24){

 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross,i))))
 X <- 1:length(Y)

 plot(X,Y, xlab=paste('chr',i), ylab='physical position',cex.axis=2,pch=21,cex.main=2)


 }
dev.off()
################################################################################

#### SEG DISTORTION X PHYSICAL POSITION#########################################
gt <- geno.table(cross_nbh)

pdf("/home/jmiller1/public_html/NBH_segdist.pdf", width=21, height=21)
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
pdf("/home/jmiller1/public_html/NBH_ELR_segdist.pdf", width=1750, height=1750)
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
## gt.nbh, gt.elr, scan_elr,

pdf("/home/jmiller1/public_html/NBH_ELR_segdist_binary.pdf", width=1750, height=1750)
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

 plot(c(0,xmax),c(0,ymax),ylim=c(0,4), type="n", xlab=paste('chr',i,'phys position'), ylab='-log10 p.value',cex.lab=2,cex.sub=2, cex.axis=2,cex.main=2)
 points(X,Ylod, col='cornflowerblue',pch=21)
 points(A,Blod, col='pink',pch=21)
 points(X,Y, col='blue',pch=19)
 points(A,B, col='red',pch=19)

 }
dev.off()


################################################################################

plot_ef <- function(crs,map,pr,ahr,popgen,chs,main,model=c("bin","pheno_norm"),...){

 for (chr in chs){

  c2eff <- scan1coef(pr[,as.character(chr)], crs$pheno[,model])

  plot(c2eff, map[as.character(chr)], columns=1:3, col=col, ylim=c(0,1), cex.axis = 2,main=main,...)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      abline(v=as.numeric(ahr[indx,'pos1']), col='red',lwd=0.5)
      #xleft, ybottom, xright, ytop,

    }
    #if(any( chr %in% popgen$chr )) {
    #  indx <- which(popgen$chr %in% chr)
    #  abline(v=as.numeric(popgen[indx,'pos1']), col='red')
    #}


  last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients

  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
  }

}


################################################################################
################################################################################

plot_pgen <- function(crs,chrs,stat, map, ahr, ahr_clm, colnm, popgen, ylimo,rank_clm,stat_name,...){

 for (chr in chrs){

  xl <- summary(pull.map(crs))[chr,'length']
  ind <- which(stat$chr == chr)

  Y <- stat[ind,colnm]
  X <- stat[ind,map]/1000000
##  plot(X, Y, col='blue', cex.axis = 2, ylim = ylimo, xlim = c(0,xl), main=paste('CHR',chr), cex.main=2)

  plot(X, Y, col='black',type="n",xlim=c(0,max(X)), ylim = ylimo, main=NULL,
   xlab='physical position', ylab=stat_name, xaxs="i",yaxs="i", mgp = c(1, 0.5, 0),...)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      #rect(ahr[indx,ahr_clm]/1000000,ylimo[1],ahr[indx,'stp']/1000000,ylimo[2],lwd=0.5,col=alpha('lightgrey',.5))
      abline(v=as.numeric(ahr[indx,ahr_clm])/1000000,
       col='red',lwd=0.5)
    }

    if(any( chr %in% popgen$chr )) {
      indx <- which(popgen$chr %in% chr)
      rect(popgen[indx,'start']/1000000,ylimo[1],popgen[indx,'end']/1000000,ylimo[2],
       border = NA,lwd=0,col=alpha('lightgrey',.5))

      #abline(v=as.numeric(popgen[indx,rank_clm])/1000000, col='grey',lwd=2)
    }
  points(X, Y, col='black',...)
 }
}



#
######

a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_ELR)[[X]], 0.25)} ))
a <- pull.markers(cross_ELR,a)
a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

pdf("/home/jmiller1/public_html/elr_rf.pdf", width=1000, height=14)
 plotRF(rf, zmax=8, col.scheme="redblue",what='lod')
dev.off()



a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_NBH)[[X]], 0.25)} ))
a <- pull.markers(cross_NBH,a)
#a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

pdf("/home/jmiller1/public_html/nbh_rf.pdf", width=1000, height=14)
 plotRF(rf, zmax=10, col.scheme="redblue",what='lod')
dev.off()

################################################################################
save.image(file.path(mpath,'supplemental_plot_env.rsave'))
################################################################################
