#!/bin/R
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

################################################################################
## read in the QTL cross
mpath <- '/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops'
fl <- 'NBH.um.unmapped.f2.csvr'
cross <- read.cross(file = file.path(mpath, fl),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
################################################################################

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

#### SWITCH ALLELES THAT ARE PROB AA x BB ######################################
m <- which(cross$pheno$ID=='NBH_NBH1M')
f <- which(cross$pheno$ID=='NBH_NBH1F')
bfixm <- pull.geno(cross)[m,]
bfix_swit1 <- names(bfixm)[which(as.numeric(bfixm)==1)]
bfixf <- pull.geno(cross)[f,]
bfix_swit2 <- names(bfixf)[which(as.numeric(bfixf)==3)]
bfix_swit12 <- unique(c(bfix_swit1 ,bfix_swit2))
cross <- switchAlleles(cross, markers = bfix_swit12)
################################################################################

################################################################################
## Parent markers
## AA, AB, BB are displayed in the colors red, blue, and green,
################################################################################
parc1 <- subset(cross,ind=c('NBH_NBH1M','NBH_NBH1F'))
plot_test('par_nbh_unfilt', width=4000, height = 5000)
par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc1, chr=i)} ; dev.off()
################################################################################

crossbk <- cross

################################################################################
## drop invariant and ABxAB cross in grand parents
m <- which(cross$pheno$ID=='NBH_NBH1M')
f <- which(cross$pheno$ID=='NBH_NBH1F')
pars <- which(cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))

bfixbk <- pull.geno(cross)
drop <- names(which(bfixbk[m,] == bfixbk[f,]))
drop.na <- names(which(is.na(bfixbk[m,]) & is.na(bfixbk[f,])))
table(bfixbk[pars,drop])
cross <- drop.markers(cross,c(drop.na,drop))
################################################################################

################################################################################
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.50))
### is "NBH_5646" another grandparent sample??
toss.missing <- c(toss.missing,"NBH_5646")
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers on all LGs)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 3.16e-4),])
cross <- drop.markers(cross,toss)
################################################################################

################################################################################
### ALL BUT 2, 13, and 18 can be filtered down to 9e-3. Truncates these LGS (see plots)
gt.sub <- geno.table(cross,chr=c(1,3:12,14:24))
toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 9e-3),])
cross <- drop.markers(cross,toss.sub)
################################################################################

crossbk2 <- cross

cross <- thin_by_distortion(cross,10)

################################################################################
### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_filtered_pvalue_NW')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
################################################################################

################################################################################
parc2 <- pull.markers(parc1,markernames(cross))
plot_test('par_nbh_unfilt', width=4000, height = 5000)
par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc2, chr=i)} ; dev.off()
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.25)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

################################################################################
### READ THE UNMAPPED MARKER ASSIGNMENT TABLE
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
move <- move[which(move$nw_marks_assign %in% markernames(cross)),]

### ASSIGN UNMAPPED MARKERS
for (i in 1:length(move[,1])){
 cross <<- movemarker(cross, marker = move[i,'nw_marks_assign'], newchr = move[i,'nw_ch'], newpos = as.numeric(move[i,'nw_pos']))
 print(i)
}
cross <- subset(cross,chr=1:24)
################################################################################

### WRITE THE ABOVE CROSS OBJECT ###############################################
mapfile <- paste0(pop,'_filtered_unphased_NW_moved')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
################################################################################







################################################################################
### Switch phase and keep only parent conf markers #############################
### ENRICH FOR AAxBB ##########################################################

###################################################################################
###### Get highly likely (Parent) AB x AB markers ##########################################
##bfix <- pull.geno(cross)
##m <- which(cross$pheno$ID=='NBH_NBH1M')
##f <- which(cross$pheno$ID=='NBH_NBH1F')
##bfixm <- names(which(bfix[m,]==3))
##bfixf <- names(which(bfix[f,]==1))
##
#### Higher confidence markers (GP genotyped as AA x BB)
##high_parABxAB <- intersect(bfixm,bfixf)
#### Lower confidence markers (at least one GP genotyped)
##low_parABxAB <- unique(bfixf,bfixm)
##
##crossh <- pull.markers(cross, high_parABxAB)
##crossl <- pull.markers(cross, low_parABxAB)
##plot_test('nbh_remove_AB_h', width=1000, height = 500); geno.image(crossh, chr=2);dev.off()
##plot_test('nbh_remove_AB_l', width=1000, height = 500); geno.image(crossl, chr=2);dev.off()
##
##high_confid <- intersect(bfixm,bfixf)
##low_confid <- low_parABxAB[!low_parABxAB %in% high_parABxAB]
##
#### PULL MARKERS THAT HAVE AT LEAST ONE PARENT GENOTYPED
##cross_highC <- pull.markers(cross,high_confid)
##
#### USE HIGH AND LOW TO GET SUFFICIENT MARKERS
##cross <- pull.markers(cross,unique(c(high_confid,low_confid)))
##
##################################################################################
##parc3 <- subset(cross_highC,ind=c('NBH_NBH1M','NBH_NBH1F'))
##plot_test('par_nbh_filt2', width=4000, height = 5000)
##par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc3, chr=i)} ; dev.off()
##################################################################################

##################################################################################
#### TOSS PARENTS AND HIGH MISSING DATA SAMPLES TO MAP
##toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.25))
####toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137","NBH_5646")
##### is "NBH_5646" another grandparent sample??
##toss.missing <- c(toss.missing,"NBH_5646")
##
##cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
##cross_highC <- subset(cross_highC, ind=!cross_highC$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
##################################################################################
##
##################################################################################
##### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
##misg <- function(X,perc) { nind(cross) * perc }
##mis <- misg(cross,0.50)
##drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
##cross <- drop.markers(cross,drop)
##
##drop <- names(which(colSums(is.na(pull.geno(cross_highC))) > mis))
##cross_highC <- drop.markers(cross_highC,drop)
##################################################################################

### PLOTS ######################################################################
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross)
plot_test('premap_nbh_mar_regression_low_conf', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 abline(h=4)
 plot(c(1,length(X)),c(0,max(Y)),type="n", ylab='physical position')
  points(X,Y)
dev.off()
################################################################################

### PLOTS ######################################################################
sm <- scanone(cross_thin10, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross_thin10))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross_thin10)
plot_test('premap_nbh_mar_regression_no_conf', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 abline(h=4)
 plot(c(1,length(X)),c(0,max(Y)),type="n", ylab='physical position')
  points(X,Y)
dev.off()
################################################################################


##cross_res <- subset(cross, ind = cross$pheno$Pheno == 0)
##cross_sen <- subset(cross, ind = cross$pheno$Pheno == 5)
##
##geno.table(cross_res,13)
##
##### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers)
##
##gt <- geno.table(cross)
##toss <- rownames(gt[which(gt[,'P.value'] < 3.16e-4),])
##cross <- drop.markers(cross,toss)
##
##### ALL BUT 2 and 18 can be filtered down to 1.0e-3
##gt.sub <- geno.table(cross,chr=c(1,3:12,14:24))
##toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 3.16e-3),])
##cross <- drop.markers(cross,toss.sub)
##
##################################################################################
#### TAKE THE LEAST DISTORTED MARKER FROM EACH READ (best marker per 500bp)

#cross <- thin_by_distortion(cross,5)

################################################################################

### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_filtered_unphased_NW')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

## WRITE THE MARKER TABLE
mfh <- file.path(mpath,'NBH_high_confid.tsv')
mfl <- file.path(mpath,'NBH_low_confid.tsv')
write.table(high_confid,mfh)
write.table(low_confid,mfl)
##########################

## TEST WHICH UNMAPPED SCAFFOLDS ARE LINKED
nw_chr <- grep('NW_',chrnames(cross), value = T)
chr <- 1:24

cross <- est.rf(cross)

##########################
ldm <- function(nw) {
 sapply(chr,function(z){ mean(pull.rf(cross, what='lod')[markernames(cross,nw),markernames(cross,z)]) })
}
##########################

library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
ld <- foreach(nw = nw_chr, .inorder = F, .packages = libs2load) %dopar% ldm(nw)
ld <- do.call(rbind,ld)
rownames(ld) <- nw_chr

nms <- which(apply(ld,1,max,na.rm=T) > 5)
reassign <- apply(ld,1,which.max)
reassign <- reassign[nms]

nw_marks_assign <- sapply(names(reassign),markernames,cross = cross)
nw_length <- sapply(nw_marks_assign,length)
nw_marks_assign <- as.character(unlist(nw_marks_assign))
nw_ch <- rep(as.numeric(reassign), times = as.numeric(nw_length))
nw_pos <- unlist(sapply(nw_length,seq,from = 1, by = 1))
nw_old <- gsub(":.*","",nw_marks_assign)

### REASSING THE UNMAPPED MARKERS
move <- data.frame(cbind(nw_old,nw_marks_assign,nw_ch,nw_pos), stringsAsFactors=F)
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
write.table(move,movefl)
################################################################################

################################################################################
### READ IN THE CROSS
fl <- paste0(pop,'_filtered_unphased_NW.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

### READ THE UNMAPPED MARKER ASSIGNMENT TABLE
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
move <- move[which(move$nw_marks_assign %in% markernames(cross)),]

a <- c('NW_012224817.1','NW_012225741.1','NW_012224621.1')
try[which(move$nw_old %in% a),]




### ASSIGN UNMAPPED MARKERS
for (i in 1:length(move[,1])){
 cross <<- movemarker(cross, marker = move[i,'nw_marks_assign'], newchr = move[i,'nw_ch'], newpos = as.numeric(move[i,'nw_pos']))
 print(i)
}
cross <- subset(cross,chr=1:24)
################################################################################

### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_filtered_unphased_NW_moved')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

 ### MAKE LOW/NO RECOMB GROUPS FOR CLEANUP #####################################
 ## high LOD initial check phase ###############################################
badmarks <- function(X){
 crossX <- est.rf(subset(cross,chr=X))
 RF <- 10/nind(crossX)
 LOD <- 12
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 linked1 <- markernames(crossX,chr=1)
 linked2 <- markernames(crossX,chr=2)
 #### REMOVE NON AB AB ###################################################

 dist <- sapply(chrnames(crossX), function(X) { mean(-log10(geno.table(crossX, chr=X)$P.value)) })

 if(any(dist > 4)){
  drops <- names(which(dist > 4))
  drops <- markernames(crossX,chr=drops)
  list(linked1,linked2,drops)
 } else { list(linked1,linked2,drops=NA) }
}

cl <- makeCluster(20)
registerDoParallel(cl)
drops <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% badmarks(X)


high <- unlist(lapply(drops,"[[",1))
low <- unlist(lapply(drops,"[[",2))

cross.sub <- pull.markers(cross,c(unlist(high),unlist(low)))
nw_chr <- grep('NW_',chrnames(cross.sub), value = T)

cross.sub <- drop.markers(cross,nw_chr)
cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
cross.sub <- tspOrder(cross = cross.sub, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross.sub <- removeDoubleXO(cross.sub)
mis <- misg(cross.sub,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross.sub))) > mis))
cross.sub <- drop.markers(cross.sub,drop)
keep <- markernames(cross.sub)

cross1 <- pull.markers(cross,c(keep,nw_chr))
cross1 <- tspOrder(cross = cross1, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross1 <- removeDoubleXO(cross1)
mis <- misg(cross1,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross1))) > mis))
cross1 <- drop.markers(cross1,drop)
cross1 <- tspOrder(cross = cross1, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

################################################################################
drop <- c('15:21481705','20:16313414','22:8254649','24:3192380')
cross1 <- drop.markers(cross1,drop)
plotit(cross1)
################################################################################

### WRITE MAP
mapfile <- paste0(pop,'_',sum(nmar(cross1)),'_noimpute_NW_remapped_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross1,filestem=filename,format="csv")
################################################################################

################################################################################
### FIX GENOTYPING ERRORS WITH GENOPROB
cross2 <- fill.geno(cross1, method="maxmarginal", error.prob = 0.05, min.prob=0.95)
mis <- misg(cross2,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross2))) > mis))
cross3 <- drop.markers(cross2,drop)

cross4 <- est.rf(cross3)
cross4 <- tspOrder(cross = cross4, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross4 <- fill.geno(cross4, method="no_dbl_XO", error.prob = 0.01, min.prob=0.99)

crs <- formLinkageGroups(cross4, max.rf = 0.15, min.lod = 15, reorgMarkers = TRUE)
drop <- rownames(geno.table(crs,chr=25:nchr(crs)))
cross5 <- drop.markers(cross4,drop)
cross5 <- tspOrder(cross = cross5, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

### WRITE IMPUTED AND CORRECTED MAP
mapfile <- paste0(pop,'_',sum(nmar(cross5)),'_imputed_NW_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross5,filestem=filename,format="csv")
################################################################################
################################################################################



################################################################################
