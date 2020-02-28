#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

mpath <- '/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops'
fl <- 'NBH.um.unmapped.f2.csvr'
################################################################################
## read in the QTL cross
cross <- read.cross(file = file.path(mpath, fl),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
mpath <- '/home/jmiller1/QTL_agri/data'

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

#################################################################################
### read in the QTL cross
#cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
#format = "csvr", geno = c(1:3), estimate.map = FALSE)
#################################################################################

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
## drop invariant and ABxAB cross in grand parents
m <- which(cross$pheno$ID=='NBH_NBH1M')
f <- which(cross$pheno$ID=='NBH_NBH1F')
bfixbk <- pull.geno(cross)
drop <- names(which(bfixbk[m,] == bfixbk[f,]))
pars <- which(cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
table(bfixbk[pars,drop])
cross <- drop.markers(cross,drop)
################################################################################

################################################################################
### Switch phase and keep only parent conf markers #############################
### ENRICH FOR AAxBB ##########################################################

## DROP DANGEROUS ABxAB cross ##################################################
## ALREADY DROPPED IN INVARIANT FILTER
################################################################################

#### SWITCH ALLELES THAT ARE PROB AA x BB ######################################
bfixm <- pull.geno(cross)[m,]
bfix_swit1 <- names(bfixm)[which(as.numeric(bfixm)==1)]
bfixf <- pull.geno(cross)[f,]
bfix_swit2 <- names(bfixf)[which(as.numeric(bfixf)==3)]
bfix_swit12 <- unique(c(bfix_swit1 ,bfix_swit2))
cross <- switchAlleles(cross, markers = bfix_swit12)
################################################################################

## Parent markers
parc <- subset(cross,ind=c('NBH_NBH1M','NBH_NBH1F'))
plot_test('par_nbh', width=4000, height = 5000)
par(mfrow = c(24,1))
for(i in 1:24){ geno.image(parc, chr=i)} ; dev.off()

#################################################################################
#### Get highly likely (Parent) AB x AB markers ##########################################
bfix <- pull.geno(cross)
m <- which(cross$pheno$ID=='NBH_NBH1M')
f <- which(cross$pheno$ID=='NBH_NBH1F')
bfixm <- names(which(bfix[m,]==3))
bfixf <- names(which(bfix[f,]==1))

## Higher confidence markers
high_parABxAB <- intersect(bfixm,bfixf)
## Lower confidence markers
low_parABxAB <- unique(bfixf,bfixm)

crossh <- pull.markers(cross, high_parABxAB)
crossl <- pull.markers(cross, low_parABxAB)
plot_test('nbh_remove_AB_h', width=1000, height = 500); geno.image(crossh, chr=2);dev.off()
plot_test('nbh_remove_AB_l', width=1000, height = 500); geno.image(crossl, chr=2);dev.off()

parABxAB <- low_parABxAB
#parABxAB <- high_parABxAB
cross <- pull.markers(cross, parABxAB)

high_confid <- intersect(bfixm,bfixf)
low_confid <- low_parABxAB[!low_parABxAB %in% high_parABxAB]
################################################################################
## TOSS PARENTS AND HIGH MISSING DATA SAMPLES
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.20))
##toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137","NBH_5646")
## is "NBH_5646" another grandparent sample??
toss.missing <- c(toss.missing,"NBH_5646")
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################

### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

### PLOTS ######################################################################
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross)
plot_test('nbh_mar_regression_hi_confid', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 abline(h=6)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
  points(X,Y)
dev.off()
################################################################################

 nw_marks <- grep('NW_',markernames(cross), value = T)
 toss <- rownames(gt[nw_marks,][which(gt[nw_marks,'P.value'] < 1.0e-5),])
 cross <- drop.markers(cross,toss)

 nw_marks <- grep('NW_',chrnames(cross), value = T)
 cross_NW <- subset(cross, chr=nw_marks)

 RF <- 0.05
 LOD <- 20
 cross_NW <- formLinkageGroups(cross_NW, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 chr <- names(which(nmar(cross_NW) > 3))
 mar <- markernames(cross_NW,chr=chr)
 mar <- c(mar,markernames(cross,chr=1:24))
 cross <- pull.markers(cross, mar)
 nw_marks <- grep('NW_',markernames(cross), value = T)


################################################################################
### WRITE THE ABOVE CROSS OBJECT
#mapfile <- paste0(pop,'_filtered_unphased_NW')
#filename <- file.path(mpath,mapfile)
#write.cross(cross,filestem=filename,format="csv")

### READ IN THE MARKER TABLE
#mfh <- file.path(mpath,'NBH_high_confid.tsv')
#mfl <- file.path(mpath,'NBH_low_confid.tsv')
#write.table(high_confid,mfh)
#write.table(low_confid,mfl)


pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

### READ IN THE CROSS
fl <- paste0(pop,'_filtered_unphased_NW.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

### READ IN THE MARKER TABLE
mfh <- file.path(mpath,'NBH_high_confid.tsv')
mfl <- file.path(mpath,'NBH_low_confid.tsv')
low_confid <- as.character(read.table(mfl)[,1])
high_confid <- as.character(read.table(mfh)[,1])
################################################################################
## THIN TO ONE MARKER PER KB
cross <- thin_by_distortion(cross,1)
################################################################################

nw_marks <- grep('NW_',markernames(cross), value = T)
nw_chr <- grep('NW_',chrnames(cross), value = T)

mapit_noimpute <- function(i){

 erprob <- 0.05
 Z <- i

 cross1 <- subset(cross,chr=c(i,nw_chr))
 #cross1 <- use_phys_map(cross1)
 cross1 <- est.rf(cross1)
 premap <- high_confid[high_confid %in% markernames(cross1,i)]

 ### MAKE LOW/NO RECOMB GROUPS FOR CLEANUP #####################################
 ## high LOD initial check phase ###############################################
 LOD <- quantile(pull.rf(cross1,what='lod',chr=i)[premap,premap],probs = 0.90, na.rm=T)
 RF <- quantile(pull.rf(cross1,chr=i)[premap,premap],probs = 0.10, na.rm=T)

 cross2 <- formLinkageGroups(cross1, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 ####################################################################

 ### PVAL filt #################################################################
 ### Allow more distortion in markers that were confrimed in the g.parents (1.0e-5)
 ### Filter more strictly for the markers that might be ABxAA
 gt <- geno.table(cross2)
 pval.h <- 1.0e-6
 pval.l <- 1.0e-4
 mis <- misg(cross2,0.15)
 bfixH <- rownames(gt[which(gt$P.value > pval.h & gt$missing < mis),])
 bfixH <- high_confid[high_confid %in% bfixH]

 bfixL <- rownames(gt[which(gt$P.value > pval.l & gt$missing < mis),])
 bfixL <- low_confid[low_confid %in% bfixL]

 bfixA <- unique(c(as.character(bfixH),as.character(bfixL)))
 cross3 <- pull.markers(cross2,bfixA)
 ###############################################################################

 ### REMOVE NON AB AB ###################################################
 dist <- sapply(chrnames(cross3), function(X) { mean(-log10(geno.table(cross3, chr=X)$P.value)) })
 keep <- names(which(dist < 4))
 cross3 <- subset(cross3,chr=keep)
 ####################################################################
 ###############################################################################
 # Filter singletons by more strict threshold
 gt <- geno.table(cross3)
 keep2 <- names(which(nmar(cross3) < 3))
 toss <- which(gt[markernames(cross3,keep2),'P.value'] < 1.0e-3)
 cross4 <- drop.markers(cross3,markernames(cross3,keep2)[toss])
 ###############################################################################

 ###############################################################################
 chr <- names(which(table(gt[premap,'chr']) > 5))
 test <- chrnames(cross4)[!chrnames(cross4) %in% chr]
 lm <- sapply(test,function(z){ mean(pull.rf(cross4, what='lod')[markernames(cross4,chr),markernames(cross4,z)]) })
 drop <- names(which(lm < 2))
 cross5 <- drop.markers(cross4,markernames(cross4,drop))
 ###############################################################################

 ## RETURN ALL MARKERS TO THE SAME LG ##########################################
 gt <- geno.table(cross5)
 chri <- names(which(table(gt[premap,'chr']) > 5))

 ## use 50th percentile lod and rf to regroup all markers
 LOD <- quantile(unlist(pull.rf(cross5,what='lod',chr=chri)),probs = 0.50, na.rm=T)
 RF <- quantile(unlist(pull.rf(cross5,chr=chri)),probs = 0.50, na.rm=T)

 cross6 <- formLinkageGroups(cross5, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross6 <- subset(cross6, chr=names(which.max(nmar(cross6))))
 ###############################################################################

 ### USE PHYSICAL COORDS TO INITIALIZE
 #cross2 <- use_phys_map(cross2)
 #ord <- order(as.numeric(unlist(pull.map(cross2))))
 #cross2 <- switch.order(cross2, chr = 1, ord, error.prob = 0.01, map.function = "kosambi", maxit = 1, tol = 0.1, sex.sp = F)

 ### REMOVE MARKERS WITH FREQUENTLY HIGH LOD ERROR
 #erl <- calc.errorlod(cross2,error.prob=0.05, map.function="kosambi")
 #terl <- top.errorlod(erl, cutoff=4)

 ### REMOVE GENOTYPES WITH HIGH LOD ERROR
 #print(paste('dropped',length(terl[,1]),'genotypes due to error lod'))

 #if(dim(terl)[1] > 0)
 #for(m in 1:nrow(terl)) {
  #chr <- terl$chr[m]
  #id <- terl$id[m]
  #mar <- terl$marker[m]
  #cross2$geno[[chr]]$data[cross2$pheno$ID==id, mar] <- NA
 #}

 ###############################################################################
 cross7 <- tspOrder(cross = cross6, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

 mapfile <- paste0(pop,'_reorder_noimput_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross7,filestem=filename,format="csv")
###############################################################################
###############################################################################

###############################################################################
###############################################################################
 cross8 <- removeDoubleXO(cross7)
 cross8 <- fill.geno(cross8, method="no_dbl_XO", error.prob = 0.05, min.prob=0.99)

 ## REMOVE MARKERS WITH HIGH RATE OF MISSING DATA
 mis <- misg(cross8,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross8))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross9 <- drop.markers(cross8, bfixA)
 ###############################################################################

 cross10 <- tspOrder(cross = cross9, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

 map <- as.numeric(pull.map(cross10)[[1]])
 if(cor(1:length(map),map, use="complete.obs") < 0) cross10 <- flip.order(cross10, chrnames(cross10))

 mapfile <- paste0(pop,'_reorder2_noimput_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross10,filestem=filename,format="csv")
###############################################################################
###############################################################################

###############################################################################
###############################################################################
 cross11 <- fill.geno(cross10, method="maxmarginal", error.prob = 0.05, min.prob=0.95)
 cross11 <- fill.geno(cross11, method="no_dbl_XO", error.prob = 0.05, min.prob=0.95)

 mapfile <- paste0(pop,'_order_impute_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross11,filestem=filename,format="csv")
###############################################################################
###############################################################################
 ## TEST WHETHER UNMAPPED MARKERS ARE LINKED
 #goodmarks <- markernames(cross2)
 #cross_NW <- pull.markers(cross,c(goodmarks,nw_marks))
 #cross_NW <- subset(cross_NW,chr= names(which(nmar(cross_NW)>3)))
 #cross_NW <- removeDoubleXO(cross_NW)

 #RF <- 0.1
 #LOD <- 16
 #cross_NW <- formLinkageGroups(cross_NW, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 #cross_NW <- subset(cross_NW, chr=1)
 #cross_NW <- tspOrder(cross = cross_NW, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

 #cross_NW <- fill.geno(cross_NW, method="maxmarginal", error.prob = 0.05, min.prob=0.99)
 #cross_NW <- fill.geno(cross_NW, method="no_dbl_XO", error.prob = 0.05, min.prob=0.99)

 #mapfile <- paste0(pop,'_order_impute_NW_',i,'_tsp')
 #filename <- file.path(mpath,mapfile)
 #write.cross(cross_NW,filestem=filename,format="csv")

 plotit(cross7,nme='no_imp')
 plotit(cross10,nme='reorder')
 plotit(cross11,nme='reorder_imp')
 #plotit(cross_NW,nme='reorder_imp_NW')
 print(paste('done with chr',i))
}

library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% mapit_noimpute(i)
#a <- c(11,10,9,13,12,16)

png(paste0('~/public_html/',pop,'all_phys.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
 fl <- paste0(pop,'_reorder_noimput_',B,'_tsp.csv')
 cross_plot <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross_plot))))/1000000
 X <- 1:length(Y)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',B), ylab='physical position')
 points(X,Y)
}
dev.off()


png(paste0('~/public_html/',pop,'all_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
 fl <- paste0(pop,'_reorder_noimput_',B,'_tsp.csv')
 cross_plot <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
 plotRF(cross_plot)
}
dev.off()


plotit(cross5)

#######################################################################################
#######################################################################################

 #RF <- 7/nind(cross2)
 #LOD <- 14
 #cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ###############################################################################
 #keep <- chrnames(cross2)
 #lm <- sapply(keep[-1],function(z){ max(pull.rf(cross2, what='lod')[markernames(cross2,keep[1]),markernames(cross2,z)]) })
 #drop <- names(which(lm < 0.1))
 #cross2 <- drop.markers(cross2,markernames(cross2,drop))
 ###############################################################################

 ###############################################################################
 #rf <- pull.rf(est.rf(cross2))
 #keep <- chrnames(cross2)
 #minrf <- sapply(keep[-1],function(z){ min(rf[markernames(cross2,keep[1]),markernames(cross2,z)]) })
 #chr <- names(which(minrf > 0.25))
 #cross2 <- drop.markers(cross2,markernames(cross2,chr=chr))
 ###############################################################################
