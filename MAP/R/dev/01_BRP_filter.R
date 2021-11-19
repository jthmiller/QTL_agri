#!/bin/R

pop <- 'BRP'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

## Map a smaller subset to id qtls
## go back and remap chrs with a qtl with a larger sample set
#################################################################################
### read in the QTL cross
cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
#################################################################################

################################################################################
### Pull sample names from plinkfile
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

## DROP PARENTS ################################################################
pars <- c('BRP_BRP1M','BRP_BRP8F','BRP_BRP1F','BRP_BRP8M')
parc <- subset(cross, ind=pars)
cross <- subset(cross, ind=!cross$pheno$ID %in% pars)
################################################################################

################################################################################
### FILTER
################################################################################
## Drop samples that have high missing data or discordant genotype
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.25))
cross <- subset(cross, ind=!cross$pheno$ID %in% toss.missing)
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X, perc) { nind(cross) * perc }
## Drop markers with greater than 12.5% missing data
mis <- misg(cross,0.125)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
print(paste('dropping',length(drop),'markers'))
cross <- drop.markers(cross,drop)
################################################################################

sum(nmar(cross))

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers on all LGs)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 5.0e-2),])
cross <- drop.markers(cross,toss)
i <- 5 ; plotit(cross,'post-dist')
################################################################################

sum(nmar(cross))

################################################################################
mapfile <- paste0(pop,'_unmapped_filtered')
fl <- file.path(mpath,mapfile)
write.cross(cross,filestem=fl,format="csv")
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

sum(nmar(cross))

################################################################################
linked_marks <- function(cross, X, LOD = 10, RF = 1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=c(1,2))
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- pull.markers(cross, unlist(linked))
i <- 95 ; plotit(cross,'cross_linked')
################################################################################

#X <- 8
#crossX <- est.rf(subset(cross, chr=X))
#LOD = 10; RF = 0.05
#crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
#i <- 8 ; plotit(crossX,'cross_linked')

################################################################################
mapfile <- paste0(pop,'_unmapped_filtered_linked')
fl <- file.path(mpath,mapfile)
#write.cross(cross,filestem=fl,format="csv")
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

sum(nmar(cross))

################################################################################
linked_marks <- function(cross, X, LOD = 8, RF = 0.1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- switchAlleles(cross, unlist(linked))
i <- 95 ; plotit(cross,'cross_linked_1')

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- switchAlleles(cross, unlist(linked))
i <- 96 ; plotit(cross,'cross_linked_2')
################################################################################


cross.sub <- subset(cross,chr=2)
cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
rf <- pull.rf(cross.sub, what='lod')
freq <- rowSums(rf, na.rm=T)/dim(rf)[1]
plot_test('chr2_brp')
hist(as.numeric(freq), breaks = 50)
dev.off()

################################################################################
switch.phase <- function(i, cross){
 cross.sub <- subset(cross,chr=i)
 cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
 rf <- pull.rf(cross.sub)
 freq <- rowSums(rf, na.rm=T)/dim(rf)[1]
 return(names(freq)[which(freq > 0.40)])
}
################################################################################

################################################################################
switch.many <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% switch.phase(i, cross = cross)
cross <- switchAlleles(cross, as.character(unlist(switch.many)))
i <- 95 ; plotit(cross,'cross_linke_sw')
################################################################################

################################################################################
linked_marks <- function(cross, X, LOD = 8, RF = 1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=c(1))
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- pull.markers(cross, unlist(linked))
i <- 95 ; plotit(cross,'cross_linked')
################################################################################


 cross.sub <- subset(cross,chr=2)
 cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
 rf <- pull.rf(cross.sub)
 freq <- rowSums(rf, na.rm=T)/dim(rf)[1]
plot_test('chr2_brp')
hist(as.numeric(freq), breaks = 50)
dev.off()















pop <- 'BRP'

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
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



##################################################################################

mfl <- file.path(mpath,'NBH_markernames.tsv')
nbh_marks <- read.table(mfl)
cross <- pull.markers(cross,nbh_marks$x)

##################################################################################
##### Switch phase and keep only parent conf markers #############################
##### ENRICH FOR AAxBB ##########################################################

#### DROP DANGEROUS ABxAB in grandparents ##################################################
## DROP1M <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1M',]
## DROP1F <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1F',]
## DROP8M <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8M',]
## DROP8F <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8F',]
##
## m <- markernames(cross)[which(DROP1M == 2 | DROP8M == 2)]
## f <- markernames(cross)[which(DROP1F == 2 | DROP8F == 2)]
## mf <- intersect(m,f)
##
## cross <- drop.markers(cross, mf)

##m <- names(DROP1M)[which(DROP1M == DROP8M)]
##f <- names(DROP1F)[which(DROP1F == DROP8F)]
##mf <- intersect(m,f)
##cross <- pull.markers(cross, mf)
##
##### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
##bfix <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8M',]
##bfix_swit1 <- names(bfix)[which(as.numeric(bfix)==1)]
##bfix <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1F',]
##bfix_swit2 <- names(bfix)[which(as.numeric(bfix)==3)]
##bfix_swit12 <- intersect(bfix_swit1 ,bfix_swit2)
##
##cross <- switchAlleles(cross, markers = bfix_swit12)
##################################################################################
##
##
##################################################################################
##### Get highly likely AB x AB markers ##########################################
##bfix1 <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8M',]
##bfix1 <- names(bfix1)[which(as.numeric(bfix1)==3)]
##bfix2 <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1F',]
##bfix2 <- names(bfix2)[which(as.numeric(bfix2)==1)]
##parABxAB <- intersect(bfix1,bfix2)
##
##gt_nopar <- geno.table(subset(cross,ind=!cross$pheno$ID %in% pars))
##parABxAB <- intersect(rownames(gt_nopar[which(gt_nopar$P.value > 0.01),]) ,parABxAB)
##cross.1 <- pull.markers(cross,parABxAB)
##################################################################################

##### TEST SAMPLE GT SIMILARITY ##################################################
##cross.1 <- subset(cross.1,ind=!cross.1$pheno$ID%in% pars)
##cpgt <- comparegeno(cross.1)
##colnames(cpgt) <- cross.1$pheno$ID
##rownames(cpgt) <- cross.1$pheno$ID
##cpgt[cpgt==NaN] <- NA
##diag(cpgt) <- NA
##cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
##################################################################################
## png(paste0('~/public_html/BRP_relat.png'))
##  hist(cpgt)
## dev.off()
################################################################################
toss.missing <- c("BRP_2535","BRP_2410","BRP_2687","BRP_2710")
################################################################################

################################################################################
#### Pvalue and Missing ##############################################
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,pars)))
bfixA <- rownames(gt[which(gt$P.value > 0.00001),])
##bfixA <- rownames(gt[which(gt$P.value > 0.00001 & gt$missing < 5),])
##bfixA <- rownames(gt[which(gt$P.value > 0.000001 & gt$missing < 5),])
################################################################################

###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,pars))
################################################################################

sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

sm <- scanone(cross, pheno.col=4, model="binary",method="mr")

plot_test('brp_mar_regression', width = 1500, height = 750)
par(mfrow=c(2,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,12), cex = 0.25)
 plot(1:length(gt[bfixA,1]), -log10(gt[bfixA,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,15), cex = 0.25)
dev.off()

################################################################################

png(paste0('~/public_html/BRP_pvals.png'))
 hist(log10(gt$missing))
dev.off()

fl <- file.path(mpath,'BRP_unmapped_filtered')
write.cross(cross,filestem=fl,format="csv")
