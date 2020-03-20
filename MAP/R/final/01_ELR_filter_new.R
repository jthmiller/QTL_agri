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
################################################################################

################################################################################
### Toss individuals that have high missing data
toss.parents <- c('BLI_BI1124M','ELR_ER1124F')
toss.badata <- c("ELR_10869","ELR_10987","ELR_11580")
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.50))
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.badata,toss.missing,pars))
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.125)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers on all LGs)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 5.0e-3),])
cross <- drop.markers(cross,toss)
################################################################################

################################################################################
### ALL BUT 17 can be filtered down to 1e-2. Truncates these LGS (see plots)
gt.sub <- geno.table(cross,chr=c(1:16,18:24))
toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
cross <- drop.markers(cross,toss.sub)
################################################################################

## RETAIN BEST (LEAST DISTORTED) MARKER PER RADTAG #############################
cross <- thin_by_distortion(cross,5)
################################################################################

### WRITE THE ABOVE CROSS OBJECT ###############################################
mapfile <- paste0(pop,'_filtered_unphased')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_filter.rsave')))
################################################################################

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
################################################################################

################################################################################
linked_marks <- function(X, LOD = 12, RF = 1){

 crossX <- est.rf(subset(cross,chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}
linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(X)
cross <- pull.markers(cross,unlist(linked))
################################################################################

################################################################################
unphased_marks <- function(X, crossX, LOD = 12, RF = 0.15){
 crossX <- est.rf(subset(crossX,chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 mk1 <- markernames(crossX, chr=1)
 crossX <- switchAlleles(crossX, mk1)
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX,1)[!markernames(crossX,1) %in% mk1]
}
switch <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% unphased_marks(X, crossX = cross)
cross <- switchAlleles(cross, unlist(switch))
################################################################################

################################################################################
switched_marks <- function(X){
 checkAlleles(subset(cross,chr=X), threshold = 3)
}
switched <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% switched_marks(X)
switched <- switched[!sapply(switched,is.null)]
switched <- do.call(rbind,switched)
cross <- switchAlleles(cross, as.character(switched$marker))
################################################################################

################################################################################
cross <- est.rf(cross)
################################################################################

################################################################################
cross.sub <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross.sub <- fill.geno(cross.sub, method="maxmarginal", error.prob = 0.05, min.prob=0.95)
plotit(cross.sub)


################################################################################
save.image(file.path(mpath,paste0(pop,'_filter.rsave')))
################################################################################






 cross1 <- subset(cross,chr=i)
 cross2 <- formLinkageGroups(cross1, max.rf = 0.2, min.lod = 10, reorgMarkers = TRUE)

 cross1 <- use_phys_map(cross1)
 cross1 <- est.rf(cross1)

 lod <- pull.rf(est.rf(cross1), what = 'lod')
 rf <- pull.rf(est.rf(cross1))

 cross2 <- formLinkageGroups(cross1, max.rf = 1, min.lod = 10, reorgMarkers = TRUE)




 ### MAKE LOW/NO RECOMB GROUPS FOR CLEANUP #####################################
 ## high LOD initial check phase ###############################################
 highlod <- markernames(formLinkageGroups(cross1, max.rf = 1, min.lod = LOD, reorgMarkers = TRUE),1)

 RF <- 5/nind(cross1)
 LOD <- 10
 cross2 <- formLinkageGroups(cross1, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 switch1 <- highlod[highlod %in% markernames(cross2,chrnames(cross2)[-1])]
 cross2 <- switchAlleles(cross2,switch1)
 cross3 <- formLinkageGroups(cross2, max.rf = 1, min.lod = LOD, reorgMarkers = TRUE)

 chk <- as.character(checkAlleles(cross2)$marker)
 cross2 <- switchAlleles(cross2,checkAlleles(cross2)$marker)

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
 lod <- pull.rf(est.rf(cross4), what = 'lod')

 chr <- chrnames(cross4)[-1]
 rf.mean <- sapply(chr, function(X) { mean(rf[markernames(cross4,chr=1), markernames(cross4,chr=X)],na.rm=T) })
 lod.mean <- sapply(chr, function(X) { mean(lod[markernames(cross4,chr=1), markernames(cross4,chr=X)],na.rm=T) })

 flips <- names(which(rf.mean > 0.5))
 switch <- markernames(cross4,flips)
 cross5 <- switchAlleles(cross4, markers = switch)

 RF <- 8/nind(cross5)
 LOD <- 10
 cross5 <- formLinkageGroups(cross5, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 keep <- markernames(cross5,1)

 list(switch,keep)
}

marks <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% mapit(i)
flips <- unlist(sapply(marks,"[[",1))
keep <- unlist(sapply(marks,"[[",2))


cross <- switchAlleles(cross,flips)
cross <- pull.markers(cross,keep)
################################################################################
################################################################################
################################################################################











################################################################################
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
sd <- 1
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
