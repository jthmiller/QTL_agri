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

 plot(X,Y, xlab=paste('chr',i), ylab='physical position')


 }
dev.off()
################################################################################
