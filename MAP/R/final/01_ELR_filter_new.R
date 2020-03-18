#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

################################################################################
## read in the QTL cross
cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
### Pull names from plinkfile
path <- file.path(mpath, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
################################################################################

#### SEX #######################################################################
sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec
################################################################################

################################################################################
### The parents are wrong in this population
pars <- c('BLI_BI1124M','ELR_ER1124F')
cross.par <- subset(cross,ind=cross$pheno$ID %in% pars)

### Toss individuals that have high missing data
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.50))
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,pars))
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.35)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 1.00e-4),])
cross <- drop.markers(cross,toss)
################################################################################

################################################################################
toss.badata <- c("ELR_10869","ELR_10987","ELR_11580")
toss.missing <- c('ELR_10967','ELR_11103','ELR_11587','ELR_11593','ELR_11115')
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,toss.badata,'BLI_BI1124M','ELR_ER1124F'))

#### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
#misg <- function(X,perc) { nind(cross) * perc }
#mis <- misg(cross,0.10)
#drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
#cross <- drop.markers(cross,drop)
#################################################################################

### PLOTS ######################################################################
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross)
plot_test('elr_mar_regression', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,8), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,8), cex = 0.25)
 abline(h=4)
 plot(c(1,length(X)),c(0,max(Y)),type="n", ylab='physical position', pch = 19, cex = 0.5)
  points(X,Y)
dev.off()
################################################################################

#### THIS PVALUE APPEARS TO RETAIN EVEN DISTRIBUTION
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 1e-4),])
cross <- drop.markers(cross,toss)
################################################################################

## RETAIN BEST (LEAST DISTORTED) MARKER PER RADTAG #############################
cross <- thin_by_distortion(cross,1)
################################################################################

### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_filtered_unphased')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

### READ IN THE MARKER TABLE
#mfh <- file.path(mpath,'NBH_high_confid.tsv')
#mfl <- file.path(mpath,'NBH_low_confid.tsv')
#write.table(high_confid,mfh)
#write.table(low_confid,mfl)

pop <- 'ELR'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

### READ IN THE CROSS
fl <- paste0(pop,'_filtered_unphased.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

sd <- 1

mapit <- function(i){

 erprob <- 0.05
 Z <- i

 cross1 <- subset(cross,chr=Z)
 cross1 <- use_phys_map(cross1)
 cross1 <- est.rf(cross1)

 ### MAKE LOW/NO RECOMB GROUPS FOR CLEANUP #####################################
 ## high LOD initial check phase ###############################################
 RF <- 3/nind(cross1)
 LOD <- 15
 cross2 <- formLinkageGroups(cross1, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 #### REMOVE NON AB AB ###################################################
 dist <- sapply(chrnames(cross2), function(X) { mean(-log10(geno.table(cross2, chr=X)$P.value)) })
 drops <- names(which(dist > 3))
 drops <- markernames(cross2,chr=drops)
 cross3 <- drop.markers(cross1,drops)
 RF <- 4/nind(cross3)
 LOD <- 10
 cross4 <- formLinkageGroups(cross3, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 #####################################################################

 ## Switch phase of groups that are out of phase with LG1
 rf <- pull.rf(est.rf(cross4))
 chr <- chrnames(cross4)[-1]
 rf.mean <- sapply(chr, function(X) { mean(rf[markernames(cross4,chr=1),markernames(cross4,chr=X)],na.rm=T) })
 flips <- names(which(rf.mean > 0.5))
 cross5 <- switchAlleles(cross4, markers = markernames(cross4,chr=flips))
 RF <- 4/nind(cross5)
 LOD <- 10
 cross5 <- formLinkageGroups(cross5, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ## ALL MARKERS SHOULD BE IN SINGLE LG. THIN AND ORDER
 cross6 <- subset(cross5,chr = names(which.max(nmar(cross5))))
 cross6 <- thin_by_distortion(cross6,100)
 cross6 <- tspOrder(cross = cross6, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 cross6 <- removeDoubleXO(cross6)
 ###############################################################################

 ###############################################################################
 ## REMOVE MARKERS WITH HIGH MISSING DATA
 mis <- misg(cross6,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross6))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross7 <- drop.markers(cross6, bfixA)
 #############################

 ### WRITE MAP
 mapfile <- paste0(pop,'_',sd,'_noimpute_tsp_',i)
 filename <- file.path(mpath,mapfile)
 write.cross(cross6,filestem=filename,format="csv")
 ###############################################################################

 cross8 <- fill.geno(cross7, method="no_dbl_XO", error.prob = 0.05, min.prob=0.98)
 cross8 <- fill.geno(cross8, method="maxmarginal", error.prob = 0.05, min.prob=0.98)
 #############################
 ## REMOVE MARKERS WITH HIGH MISSING DATA
 mis <- misg(cross8,0.15)
 bfixA <- names(which(colSums(is.na(pull.geno(cross8))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross8 <- drop.markers(cross8, bfixA)
 cross8 <- tspOrder(cross = cross8, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 #############################

 ### WRITE IMPUTED AND CORRECTED MAP
 mapfile <- paste0(pop,'_',sd,'_impute_tsp_',i)
 filename <- file.path(mpath,mapfile)
 write.cross(cross8,filestem=filename,format="csv")

 plotit(cross8)
 ###############################################################################
 ###############################################################################
 ###############################################################################
}


library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% mapit(i)
#a <- c(11,10,9,13,12,16)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

png(paste0('~/public_html/',pop,'all_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
 fl  <- paste0(pop,'_',sd,'_impute_tsp_',B,'.csv')
 cross_plot <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
 plotRF(cross_plot)
}
dev.off()

png(paste0('~/public_html/',pop,'all_phys.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
 fl  <- paste0(pop,'_',sd,'_impute_tsp_',B,'.csv')
 cross_plot <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross_plot))))/1000000
 X <- 1:length(Y)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',B), ylab='physical position')
 points(X,Y)
}
dev.off()
