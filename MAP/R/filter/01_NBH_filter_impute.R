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

##2 27500454 27504907      aip
# 2:27374265   2       0  1  0  1      0      0 0.36787944
# 2:27601321   2       0  1  0  1      0      0 0.36787944

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

 #toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137")
 #cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))

### NEEDS to be redone 8,10 (add scaffs?), 16(minor), 18 (still losing first 5 MB)

mapit_noimpute <- function(i){

 erprob <- 0.05
 Z <- i

 cross1 <- subset(cross,chr=Z)
 cross1 <- use_phys_map(cross1)
 cross1 <- est.rf(cross1)

 ### MAKE LOW/NO RECOMB GROUPS FOR CLEANUP #####################################
 ## high LOD initial check phase ###############################################
 RF <- 3/nind(cross1)
 LOD <- 20
 cross2 <- formLinkageGroups(cross1, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ### REMOVE NON AB AB ###################################################
 dist <- sapply(chrnames(cross2), function(X) { mean(-log10(geno.table(cross2, chr=X)$P.value)) })
 keep <- names(which(dist < 4))
 cross2 <- subset(cross2,chr=keep)
 ####################################################################

 ### PVAL filt #################################################################
 gt <- geno.table(cross2)
 pval <- 1.0e-6
 mis <- misg(cross2,0.15)
 bfixA <- rownames(gt[which(gt$P.value > pval & gt$missing < mis),])
 cross2 <- pull.markers(cross2,bfixA)
 ###############################################################################

 ## Switch phase of groups that are out of phase with LG1
 rf <- pull.rf(est.rf(cross2))
 chr <- chrnames(cross2)[-1]
 rf.mean <- sapply(chr, function(X) { mean(rf[markernames(cross2,chr=1),markernames(cross2,chr=X)],na.rm=T) })
 flips <- names(which(rf.mean > 0.5))
 cross2 <- switchAlleles(cross2, markers = markernames(cross2,chr=flips))

 ## RETURN ALL MARKERS TO THE SAME LG ##########################################
 RF <- 10/nind(cross2)
 LOD <- 8
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross2 <- subset(cross2, chr=chrnames(cross2)[1])
 ###############################################################################

 cross2 <- use_phys_map(cross2)
 ord <- order(as.numeric(unlist(pull.map(cross2))))
 cross2 <- switch.order(cross2, chr = 1, ord, error.prob = 0.01, map.function = "kosambi", maxit = 1, tol = 0.1, sex.sp = F)

 cross2 <- removeDoubleXO(cross2)
 cross2 <- fill.geno(cross2, method="no_dbl_XO", error.prob = 0.05, min.prob=0.995)

 ## Take the best marker every 1kb #############################################
 cross2 <- thin_by_distortion(cross2,10)
 ###############################################################################

 ## REMOVE MARKERS WITH HIGH MISSING DATA
 mis <- misg(cross2,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross2))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross2 <- drop.markers(cross2, bfixA)
 ###############################################################################

 #############################
 pos <- as.numeric(gsub(".*:","",markernames(cross2)))
 map <- as.numeric(pull.map(cross2)[[1]])

 if(cor(pos,map, use="complete.obs") < 0) cross2 <- flip.order(cross2, 1)
 cross2 <- shiftmap(cross2, offset=0)

 mapfile <- paste0(pop,'_unmapped_noimput_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross2,filestem=filename,format="csv")

 cross3 <- fill.geno(cross2, method="maxmarginal", error.prob = 0.05, min.prob=0.98)
 cross3 <- fill.geno(cross3, method="no_dbl_XO", error.prob = 0.05, min.prob=0.98)
 cross3 <- tspOrder(cross = cross3, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 pos <- as.numeric(gsub(".*:","",markernames(cross3)))
 map <- as.numeric(pull.map(cross3)[[1]])
 if(cor(pos,map, use="complete.obs") < 0) cross3 <- flip.order(cross3, 1)
 cross3 <- shiftmap(cross3, offset=0)
 mapfile <- paste0(pop,'_imputed_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross3,filestem=filename,format="csv")

 cross4 <- tspOrder(cross = cross2, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 mapfile <- paste0(pop,'_reordered_noImp_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross4,filestem=filename,format="csv")

 plotit(cross2,nme='no_imp')
 plotit(cross3,nme='imputed')
 plotit(cross4,nme='reorder_noimp')
}


library(qtl)
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% mapit_noimpute(i)


#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
