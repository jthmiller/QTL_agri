#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
LOD <- as.numeric(commandArgs(TRUE)[3])
RF <- as.numeric(commandArgs(TRUE)[4])
mis <- as.numeric(commandArgs(TRUE)[5])
pval <- as.numeric(commandArgs(TRUE)[6])

source("/home/jmiller1/QTL_agri/MAP/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'

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
high_parABxAB <- intersect(bfix1,bfix2)
## Lower confidence markers
low_parABxAB <- unique(bfix1,bfix2)

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

i <- 8

mapit <- function(i){

 erprob <- 0.05
 Z <- i

 cross2 <- subset(cross,chr=Z)
 toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137")
 cross2 <- subset(cross2, ind=!cross2$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
 cross2 <- est.rf(cross2)

 ### MAKE LOW/NO RECOMB GROUPS FOR CLEANUP
 ## high LOD initial check phase
 RF <- 2/nind(cross2)
 LOD <- 15
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ### REMOVE NON AB AB
 dist <- sapply(chrnames(cross2), function(X) { mean(-log10(geno.table(cross2, chr=X)$P.value)) })
 keep <- names(which(dist < 5))
 cross2 <- subset(cross2,chr=keep)

 cross2 <- subset(cross2, chr=names(which(nmar(cross2) > 1)))
 cross2 <- tspOrder(cross = cross2, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 cross2 <- removeDoubleXO(cross2)

 ### PUT BACK INTO 1 GROUP
 RF <- 20/nind(cross2)
 LOD <- 15
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross2 <- subset(cross2, chr=chrnames(cross2)[1])
 cross2 <- tspOrder(cross = cross2, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 cross2 <- removeDoubleXO(cross2)

 ca <- checkAlleles(cross2, threshold=5)
 cross2 <- switchAlleles(cross2, markers = ca$marker)

 RF <- 2/nind(cross2)
 LOD <- 15
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ## drop distorted lg

 cross2 <- fill.geno(cross2, method="no_dbl_XO", error.prob = 0.08, min.prob=0.995)

 ### PVAL filt #################
 gt <- geno.table(cross2)
 pval <- 1.0e-6
 mis <- misg(cross2,0.15)
 bfixA <- rownames(gt[which(gt$P.value > pval & gt$missing < mis),])
 cross2 <- pull.markers(cross2,bfixA)
 #############################

 cross2 <- thin_by_radtag(cross2, dist = 0.5)

 #############################
 swits <- markernames(cross2, chr= chrnames(cross2)[1])
 ## switch singletons before getting rid of them
 swits <- c(swits, markernames(cross2, names(which(nmar(cross2) < 2))))
 cross2 <- switchAlleles(cross2, markers = swits)
 RF <- 3/nind(cross2)
 LOD <- 15
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 swits <- markernames(cross2, chr=chrnames(cross2)[1])
 cross2 <- switchAlleles(cross2, markers = swits)
 #############################

 RF <- 2/nind(cross2)
 LOD <- 15
 cross30 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ca <- checkAlleles(cross30, threshold=6)
 cross30 <- switchAlleles(cross30, markers = ca$marker)
 cross30 <- formLinkageGroups(cross30, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ### REMOVES SINGLETONS
 cross30 <- subset(cross30, chr=names(which(nmar(cross30) > 1)))
 cross30 <- tspOrder(cross = cross30, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
plotit(cross30)
 ## HIGH LOD REMOVE CROSSOVERS
 cross30 <- removeDoubleXO(cross30)
 cross30 <- fill.geno(cross30, method="no_dbl_XO", error.prob = 0.05, min.prob=0.995)

 ## REMOVE MARKERS WITH HIGH MISSING DATA
 mis <- misg(cross30,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross30))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross30 <- drop.markers(cross30, bfixA)

 ## CALCULATE AVERAGE PVAL FOR EACH GROUP
 dist <- sapply(chrnames(cross30), function(X) { mean(-log10(geno.table(cross30, chr=X)$P.value)) })

 ### PLOT DISTORTION DISTRIB
 plot_test(paste0(pop,i,'dist')); hist(dist[which(nmar(cross30) > 10)], breaks=10); dev.off()

 ### Drop lod distortion greater than log10 5
 keep <- names(which(dist < 6))
 cross30 <- subset(cross30,chr=keep)

 ### Last linkage test before imputation
 RF <- 0.05
 LOD <- 20
 cross30 <- formLinkageGroups(cross30, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross30.1 <- subset(cross30, chr=names(which(nmar(cross30) > 1)))
 cross30.1 <- tspOrder(cross = cross30.1, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

 cross30.2 <- fill.geno(cross30.1, method="maxmarginal", error.prob = 0.05, min.prob=0.995)
 mis <- misg(cross30.2,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross30.2))) > mis))
 print(paste('dropping',length(bfixA),'markers due to missing data'))
 cross30.2 <- drop.markers(cross30.2,bfixA)

 swits <- markernames(cross30.2, chr=chrnames(cross30.2)[1])
 cross30.2 <- switchAlleles(cross30.2, markers = swits)
 cross30.2 <- formLinkageGroups(cross30.2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 swits <- markernames(cross30.2, chr=chrnames(cross30.2)[1])
 cross30.2 <- switchAlleles(cross30.2, markers = swits)
 cross30.2 <- formLinkageGroups(cross30.2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 cross30.2 <- subset(cross30.2, chr=names(which(nmar(cross30.2) > 1)))
 cross30.2 <- fill.geno(cross30.2, method="maxmarginal", error.prob = 0.05, min.prob=0.995)

 RF <- 0.1
 LOD <- 20
 cross30.3 <- formLinkageGroups(cross30.2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross30.3 <- subset(cross30.3, chr=names(which(nmar(cross30.3) > 1)))
 cross30.3 <- tspOrder(cross = cross30.3, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 cross30.3 <- fill.geno(cross30.3, method="no_dbl_XO", error.prob = 0.05, min.prob=0.995)

 ## REMOVE MARKERS WITH MISSING DATA
 mis <- misg(cross30.3,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross30.3))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross30.3 <- drop.markers(cross30.3, bfixA)

 ### FINAL ORDER
 cross30.4 <- formLinkageGroups(cross30.3, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross30.4 <- subset(cross30.4,chr=names(which.max(nmar(cross30.4))))
 cross30.4 <- tspOrder(cross = cross30.4, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

 pos <- as.numeric(gsub(".*:","",markernames(cross30.4)))
 map <- as.numeric(pull.map(cross30.4)[[1]])

 if(cor(pos,map, use="complete.obs") < 0) cross30.4 <- flip.order(cross30.4, 1)

 cross30.4 <- shiftmap(cross30.4, offset=0)

 mapfile <- paste0(pop,'_unmapped_all_mark_imputed_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross30.4,chr=i,filestem=filename,format="csv")
 plotit(cross30.4)
 print(paste('done with chr',i))
}

foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% mapit(i)

plotit(cross2)
plotit(cross30.2)
plotit(cross30.4)
