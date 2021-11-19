#!/bin/R

## Make a figure of unfiltered and unmapped markers to ensure that the process does not throw away a QTL 
pop <- 'ELR'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

mpath <- '/home/jmiller1/QTL_agri/data'


#### NBH #########################################################################
### IND FILTERED CROSS
mapfile <- paste0(pop,'_reorder_imp_nopar')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
keep_inds <- cross$pheno$ID 

### UNFILTERED CROSS TO SUBSET AND PLOT
mapfile <- paste0(pop,'_nopar_nofilt_refined')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath,
  genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################





################################################################################
### ELR 
mapfile <- paste0(pop,'_reorder_imp_mapped')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
keep_inds <- cross$pheno$ID 
################################################################################

### Large cross pbject
mapfile <- paste0(pop,'_samples_checked')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross <- subset(cross, ind=cross$pheno$ID %in% keep_inds)
################################################################################



################################################################################
po2 <- function(cross,nme){
 sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
 X <- 1:length(Y)
 gt <- geno.table(cross)
 plot_test(nme, width = 5500, height = 250)
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 dev.off()
}

po2(cross,'elr_test')

### PLOTS ######################################################################
po <- function(cross,nme){
 sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
 X <- 1:length(Y)
 gt <- geno.table(cross)
 plot_test(nme, width = 5500, height = 750)
 par(mfrow=c(3,1))
  plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
  abline(h=5)
  plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
  abline(h=3)
  plot(c(1,length(X)),c(0,max(Y)),type="n", ylab='physical position')
   points(X,Y)
 dev.off()
}

