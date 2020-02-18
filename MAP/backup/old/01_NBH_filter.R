#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
LOD <- as.numeric(commandArgs(TRUE)[3])
RF <- as.numeric(commandArgs(TRUE)[4])
mis <- as.numeric(commandArgs(TRUE)[5])
pval <- as.numeric(commandArgs(TRUE)[6])

source("/home/jmiller1/QTL_agri/MAP/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

################################################################################
## read in the QTL cross
cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

#arg <- paste0(pop,'_all_mark_imputed_?[0-9]?[0-9]_tsp.csv')
#file_list <- list.files(mpath, arg)
#done <- gsub('_tsp.csv','',gsub('NBH_all_mark_imputed_','',file_list))
#todo <- which(!1:24 %in% done)

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
DROP1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
DROP1 <- names(DROP1)[which(as.numeric(DROP1)==2)]
DROP2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
DROP2 <- names(DROP2)[which(as.numeric(DROP2)==2)]
DROP <- intersect(DROP1,DROP2)
cross <- drop.markers(cross,DROP)
### WHAT PERCENT? ####
## table(gsub(":.*","",DROP))/table(gsub(":.*","",markernames(cross2)))
################################################################################

#### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
bfixm <- pull.geno(cross)[m,]
bfix_swit1 <- names(bfixm)[which(as.numeric(bfixm)==1)]
bfixf <- pull.geno(cross)[f,]
bfix_swit2 <- names(bfixf)[which(as.numeric(bfixf)==3)]
bfix_swit12 <- unique(c(bfix_swit1 ,bfix_swit2))
cross <- switchAlleles(cross, markers = bfix_swit12)

parc <- subset(cross,ind=c('NBH_NBH1M','NBH_NBH1F'))
##2 27500454 27504907      aip
# 2:27374265   2       0  1  0  1      0      0 0.36787944
# 2:27601321   2       0  1  0  1      0      0 0.36787944

a <- which(markernames(parc, chr=2) == '2:27374265')
b <- which(markernames(parc, chr=2) == '2:27601321')

plot_test('par_nbh', width=2000, height = 5000)
par(mfrow = c(24,1))
for(i in 1:24){
 geno.image(parc, chr=i)
 #if(i==2) abline(v=a,lwd = 10)
 #if(i==2) abline(v=b,lwd = 10)
}
dev.off()

#################################################################################
#### Get highly likely (Parent) AB x AB markers ##########################################
bfix1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
bfix1 <- names(bfix1)[which(as.numeric(bfix1)==3)]
bfix2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
bfix2 <- names(bfix2)[which(as.numeric(bfix2)==1)]
parABxAB <- intersect(bfix1,bfix2)
cross <- pull.markers(cross, parABxAB)

plot_test('par2_nbh_AB', width=1000, height = 500); geno.image(cross, chr=1);dev.off()

#gt_nopar <- geno.table(subset(cross,ind=!cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F')))
#parABxAB <- intersect(rownames(gt_nopar[which(gt_nopar$P.value > 0.01),]) ,parABxAB)
#cross.1 <- pull.markers(cross,parABxAB)
#################################################################################

#### TEST SAMPLE GT SIMILARITY ##################################################
#cross.1 <- subset(cross.1,ind=!cross.1$pheno$ID%in%c('NBH_NBH1M','NBH_NBH1F'))
#cpgt <- comparegeno(cross.1)
#colnames(cpgt) <- cross.1$pheno$ID
#rownames(cpgt) <- cross.1$pheno$ID
#cpgt[cpgt==NaN] <- NA
#diag(cpgt) <- NA
#cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
#################################################################################
#png(paste0('~/public_html/NBH_relat.png'))
# hist(cpgt)
#dev.off()
#################################################################################
#toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137","NBH_5646")
## is "NBH_5646" another grandparent sample??
################################################################################
##cross <- subset(cross, chr=c(1,2,8,18,24))
################################################################################

##### Pvalue and Missing ##############################################
#gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F')))
#bfixA <- rownames(gt[which(gt$P.value > pval & gt$missing < mis),])
#################################################################################
#cross <- subset(cross,chr=Z)

toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137","NBH_5646")
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)

# Determine what percent of markers are kept after filter
#table(gsub(":.*","",bfixA))/table(gsub(":.*","",markernames(cross)))

#################################################################################
### Get a cross object of parent genotypes
#cross.par <- subset(cross, ind=cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
#DROP1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
#DROP1 <- names(DROP1)[which(as.numeric(DROP1)==2)]
#DROP2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
#
####### FILTER #######################################################
#cross <- pull.markers(cross,bfixA)
#cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
#################################################################################
#
#### Faster filter ##############################################################
#mfl <- paste0(pop,'prefiltered_markernames.tsv')
#mfl <- file.path(mpath,mfl)
#write.table(markernames(cross), mfl)
#
#inds <- paste0(pop,'prefiltered_indnames.tsv')
#inds <- file.path(mpath,inds)
#write.table(cross$pheno, inds)
#
#swit <- paste0(pop,'prefiltered_switch.tsv')
#swit <- file.path(mpath,swit)
#write.table(bfix_swit12, swit)
#################################################################################

#################################################################################
#
#mfl <- file.path(mpath,paste0(pop,'prefiltered_markernames.tsv'))
#marks <- read.table(mfl, stringsAsFactors=F)
#
#inds <- file.path(mpath,paste0(pop,'prefiltered_indnames.tsv'))
#inds <- read.table(inds, stringsAsFactors=F)
#
#switch <- file.path(mpath,paste0(pop,'prefiltered_switch.tsv'))
#switch <- read.table(switch, stringsAsFactors=F)
#
#cross <- subset(cross, ind=cross$pheno$ID %in% inds$ID)
#cross <- switchAlleles(cross, markers = switch)
#cross <- pull.markers(cross,marks$x)
#
##################################################################################
#sapply(1:24,function(i){
#ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))
#cross <<- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
# maxit = 1, tol = 0.1, sex.sp = F)
#})
#################################################################################

### PLOTS ######################################################################
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross)
plot_test('nbh_mar_regression', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
  points(X,Y)
dev.off()
################################################################################

mapit <- function(i){
 erprob <- 0.05
 Z <- i
 cross2 <- subset(cross,chr=Z)
 toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137")
 cross2 <- subset(cross2, ind=!cross2$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
 cross2 <- est.rf(cross2)

 ### MAKE HIGHLY LINKED GROUPS FOR CLEANUP
 ## high LOD initial check phase
 RF <- 4/nind(cross2)
 LOD <- 6
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 dist <- sapply(chrnames(cross2), function(X) { mean(-log10(geno.table(cross2, chr=X)$P.value)) })
 keep <- names(which(dist < 5))
 plot_test(paste0(pop,i,'dist')); hist(dist[which(nmar(cross2) > 10)], breaks=30) ; dev.off()

 ## drop distorted lg
 cross2 <- subset(cross2,chr=keep)
 cross2 <- fill.geno(cross2, method="no_dbl_XO", error.prob = 0.08, min.prob=0.995)

 ### PVAL filt #################
 gt <- geno.table(cross2)
 pval <- 1.0e-5
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
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 swits <- markernames(cross2, chr=chrnames(cross2)[1])
 cross2 <- switchAlleles(cross2, markers = swits)
 #############################

 cross30 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ca <- checkAlleles(cross30, threshold=6)
 cross30 <- switchAlleles(cross30, markers = ca$marker)
 cross30 <- formLinkageGroups(cross30, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ### REMOVES SINGLETONS
 cross30 <- subset(cross30, chr=names(which(nmar(cross30) > 1)))
 cross30 <- tspOrder(cross = cross30, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

 ## HIGH LOD REMOVE CROSSOVERS
 cross30 <- removeDoubleXO(cross30)

 ## REMOVE MARKERS WITH HIGH MISSING DATA
 mis <- misg(cross30,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross30))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross30 <- drop.markers(cross30, bfixA)

 ## CALCULATE AVERAGE PVAL FOR EACH GROUP
 dist <- sapply(chrnames(cross30), function(X) { mean(-log10(geno.table(cross30, chr=X)$P.value)) })
 plot_test(paste0(pop,i,'dist')); hist(dist[which(nmar(cross30) > 10)], breaks=10); dev.off()

 ### Drop lod distortion greater than log10 5
 keep <- names(which(dist < 5))
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

 RF <- 0.15
 LOD <- 15
 cross30.3 <- formLinkageGroups(cross30.2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross30.3 <- subset(cross30.3, chr=names(which(nmar(cross30.3) > 1)))
 cross30.3 <- tspOrder(cross = cross30.3, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 cross30.3 <- fill.geno(cross30.3, method="no_dbl_XO", error.prob = 0.05, min.prob=0.995)

 ## REMOVE MARKERS WITH MISSING DATA
 mis <- misg(cross30.3,0.10)
 bfixA <- names(which(colSums(is.na(pull.geno(cross30.3))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross30.3 <- drop.markers(cross30.3, bfixA)

 RF <- 0.25
 LOD <- 10
 cross30.4 <- formLinkageGroups(cross30.3, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross30.4 <-subset(cross30.4,chr=1)
 cross30.4 <- tspOrder(cross = cross30.4, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

 pos <- as.numeric(gsub(".*:","",markernames(cross30.3)))
 map <- as.numeric(pull.map(cross30.3)[[1]])

 if(cor(pos,map, use="complete.obs") < 0){
  cross30.3 <<- flip.order(cross30.3, i)
 }
 cross30.3 <- shiftmap(cross30.3, offset=0)

 mapfile <- paste0(pop,'unmapped_all_mark_imputed_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross30.3,chr=i,filestem=filename,format="csv")
}


cross30.3_map <-  est.map(cross30.3, error.prob=erprob, map.function="kosambi",maxit=10000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross30.3 <- qtl:::replace.map(cross30.3,cross30.3_map)



foreach(i = todo, .inorder = F, .packages = c("qtl")) %dopar% mapit(i)

library(qtl)
library(doParallel)
cl <- makeCluster(22)
registerDoParallel(cl)


#########################################################

cross30 <- thin_by_radtag(cross30, dist = 0.5)
plotit(cross30.3)









###### FUNCTION #################################################################
thin_by_radtag <- function(cross_in = cross30, dist = 1){
 chr <- chrnames(cross_in)
 map <- pull.map(cross_in)
 newpos <- lapply(map,function(X) { setNames(as.numeric(gsub(".*:","",names(X)))/100000,names(X))  } )
 newpos <- lapply(newpos, function(X){  class(X) <- 'A'; X } )
 attr(newpos,'class') <- 'map'
 ##attr(newpos[[chr]], "loglik") <- attr(map[[chr]], "loglik")
 names(newpos) <- chr
 cross_in <- replace.map(cross_in, newpos)
 print(summary(pull.map(cross_in)))

 ### GET ONLY 1 MARKER PER RAD TAG
 mrks <- as.numeric(gsub(".*:","",markernames(cross_in)))/100
 names(mrks) <- markernames(cross_in)
 n.missing <- nmissing(cross_in, what="mar")
 wts <- -log( (n.missing+1) / (nind(cross_in)+1) )
 a <- pickMarkerSubset(mrks, dist, wts)
 cross_in <- pull.markers(cross_in,a)
 print(nmar(cross_in))
 return(cross_in)
}
####################################################################################

plotit <- function(crs){
 Y <- c(0, as.numeric(gsub(".*:","",markernames(crs))))/1000000
 X <- 1:length(Y)
 gt <- geno.table(crs)
 sm <- scanone(crs, pheno.col=4, model="binary",method="mr")

 png(paste0('~/public_html/',pop,'_gts_test',i,'.png'),height=1500,width=2500)
 par(mfrow=c(4,1))
  plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,20), cex =1)
  plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex =1)
  crs$pheno$gtps <- (as.numeric(rowSums(pull.geno(crs) == 1 | pull.geno(crs) == 1, na.rm = T))*10) + (as.numeric(rowSums(pull.geno(crs) == 3, na.rm = T))*5)
  #crs$pheno$gtps <- rowSums(pull.geno(cross))
  geno.image(crs, reorder=6, cex=2)
  plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
  abline(h=27.504907, col='red')
  points(X,Y)
 dev.off()

 plot_test(paste0(pop,'_rf_test',i))
 plotRF(crs)
 dev.off()
}

library(qtl)
library(doParallel)
cl <- makeCluster(22)
registerDoParallel(cl)


dist <- sapply(chrnames(crs), function(X) { mean(-log10(geno.table(crs, chr=X)$P.value)) })

plot_test('dist')
hist(dist[which(nmar(crs) > 10)], breaks = 50)
dev.off()

keep <- names(which(dist < 10))
crs <- subset(crs,chr=keep)





cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

cross <- removeDoubleXO(cross, chr=1:10)

dist <- sapply(chrnames(cross), function(X) { mean(-log10(geno.table(cross, chr=X)$P.value)) })

plot_test('dist')
hist(dist[which(nmar(cross) > 10)])
dev.off()

keep <- names(which(dist < 5))
cross <- subset(cross,chr=keep)

cross1 <- cross

cross <- formLinkageGroups(cross, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
keep <- names(which(nmar(cross) > 2))
cross <- subset(cross,chr=keep)
cross <- fill.geno(cross, method="maxmarginal", error.prob = 0.08, min.prob=0.995)

swits <- markernames(cross, chr=chrnames(cross)[1])
cross <- switchAlleles(cross, markers = swits)
cross <- formLinkageGroups(cross, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
swits <- markernames(cross, chr=chrnames(cross)[1])
cross <- switchAlleles(cross, markers = swits)
cross <- formLinkageGroups(cross, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')






## TOSS INVARIANT
###cross_hom <- subset(cross,chr=c(1,2))
#cross <- subset(cross,chr=chrnames(cross)[3]:nchr(cross))
#kp <- names(which(nmar(cross) > 2))
#cross <- subset(cross,chr=kp)
#plm(cross)

rf <- pull.rf(reorg.2)
lod <- pull.rf(reorg.2, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


cross_all_um <- cross


library(qtl)
library(doParallel)
cl <- makeCluster(22)
registerDoParallel(cl)

test_link <- function(Z){

  all <- subset(cross,chr=Z)

RF <- 0.01
LOD <- 30

  reorg.2 <- formLinkageGroups(all, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

  RF <- 0.02
  LOD <- 35
  reorg.2 <- formLinkageGroups(reorg.2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
  reorg.2b <- reorg.2
swits <- markernames(reorg.2b, chr=chrnames(reorg.2b)[3])
reorg.2b <- switchAlleles(reorg.2b, markers = swits)
reorg.2b <- formLinkageGroups(reorg.2b, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)




  reorg.2a <- tspOrder(cross = subset(reorg.2,chr=1:10), hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
  reorg.2b <- fill.geno(reorg.2a, method="maxmarginal", error.prob = 0.08, min.prob=0.995)

plot_test('gi', width=1000)
geno.image(reorg.2d,chr=1:2)
dev.off()

swits <- markernames(reorg.2b, chr=chrnames(reorg.2b)[1])
reorg.2b <- switchAlleles(reorg.2b, markers = swits)
reorg.2b <- formLinkageGroups(reorg.2b, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
swits <- markernames(reorg.2b, chr=chrnames(reorg.2b)[1])
reorg.2b <- switchAlleles(reorg.2b, markers = swits)
reorg.2b <- formLinkageGroups(reorg.2b, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

  reorg.2c <- tspOrder(cross = subset(reorg.2b,chr=1:4), hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

swits <- markernames(reorg.2c, chr=chrnames(reorg.2c)[3])
reorg.2c <- switchAlleles(reorg.2c, markers = swits)
reorg.2c <- formLinkageGroups(reorg.2c, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

likely_drop <- markernames(reorg.2c, chr=chrnames(reorg.2c)[2])

reorg.2d <- tspOrder(cross = subset(reorg.2c,chr=1:2), hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

reorg.2d <- tspOrder(cross = subset(reorg.2d,chr=1:2), hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')



  reorg.2a <- reorg.2

  ##plm(reorg.2a)
  ##ro <- subset(reorg.2a,chr=c(1,2))
  ##reorg.2a <- subset(reorg.2a,chr=c(2:nchr(reorg.2a))

  reorg.2a <- removeDoubleXO(reorg.2a, chr = 1:6)
  drop <- names(which(colSums(is.na(pull.geno(reorg.2a))) > (nind(reorg.2a)*0.15)))
  length(drop)
  reorg.2 <- drop.markers(reorg.2,drop)
  print(nmar(reorg.2))






   ## switch it
  swits <- markernames(reorg.2a, chr=1)
  reorg.2a <- switchAlleles(reorg.2a, markers = swits)
  reorg.2a <- formLinkageGroups(reorg.2a, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

  print(paste('initial switch of', Z))

   # switch it back
  swits <- markernames(reorg.2a, chr=1)
  reorg.2a <- switchAlleles(reorg.2a, markers = swits)
  reorg.2a <- formLinkageGroups(reorg.2a, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

  ## added to chr 1 by switches
  orig1 <- markernames(reorg.2, chr=1)
  final <- markernames(reorg.2a, chr=1)
  added <- final[!final %in% orig1]

  new_gts <- as.matrix(reorg.2a$geno[['1']]$data[,added])
  orig_gts <- as.matrix(all$geno[[as.character(Z)]]$data[,added])

  new_gts[new_gts == 2] <- NA
  orig_gts[orig_gts == 2] <- NA

  switched <- colnames(new_gts)[which(colSums(new_gts == orig_gts, na.rm =T) == 0)]

  drops <- markernames(all)[!markernames(all) %in% final]

  print(paste('done with group', Z))

  list(switched=switched,drops=drops)

}

link <- foreach(Z = seq(along=1:24), .inorder = T, .packages = c("qtl")) %dopar% test_link (Z)

switched <- unlist(lapply(link,'[[',1))
drops <- unlist(lapply(link,'[[',2))

cross <- drop.markers(cross, drops)
cross <- switchAlleles(cross, switched)

print('done with linkage groups')
################################################################################
################################################################################


################################################################################
################################################################################
#fl <- file.path(paste0(pop,'_unmapped_filtered.csv'))
#cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

inds <- which(crossbk$pheno$ID %in% cross$pheno$ID)
old.gt <- as.matrix(pull.geno(crossbk))
new.gt <- as.matrix(pull.geno(cross))
rownames(old.gt) <- crossbk$pheno$ID
rownames(new.gt) <- cross$pheno$ID

old.gt[old.gt == 2] <- NA
new.gt[new.gt == 2] <- NA

marks <- intersect(colnames(old.gt), colnames(new.gt))
inds <- intersect(rownames(old.gt), rownames(new.gt))

old.gt <- old.gt[inds, marks]
new.gt <- new.gt[inds, marks]

switched <- colnames(old.gt)[which(colSums(old.gt == new.gt, na.rm = T) == 0)]

crossbk <- switchAlleles(crossbk, switched)

crossbk <- pull.markers(crossbk, marks)

NBH_NBH1M <- subset(crossbk, ind = crossbk$pheno$ID == 'NBH_NBH1M')
NBH_NBH1F <- subset(crossbk, ind = crossbk$pheno$ID == 'NBH_NBH1F')
marks <- intersect(markernames(NBH_NBH1M), markernames(NBH_NBH1F))
NBH_NBH1M <- as.matrix(pull.geno(NBH_NBH1M))
NBH_NBH1F <- as.matrix(pull.geno(NBH_NBH1F))


m <- lapply(1:24, function(i){
 mr <- markernames(crossbk,i)[which(markernames(crossbk,i) %in% marks)]
 names(which(NBH_NBH1M[,mr] == 1))
})

f <- lapply(1:24, function(i){
 mr <- markernames(crossbk,i)[which(markernames(crossbk,i) %in% marks)]
 names(which(NBH_NBH1F[,mr] == 3))
})

m <- unlist(m)
f <- unlist(f)

cross <- drop.markers(cross,unique(c(m,f),'18:2903774'))

################################################################################

fl <- file.path(mpath,paste0(pop,'_unmapped_filtered'))
write.cross(cross,filestem=fl,format="csv")

fl.par <- file.path(paste0(pop,'_parents_filtered'))
fl.par <- file.path(mpath,fl.par)
write.cross(cross.par,filestem=fl.par,format="csv")

##system('sbatch -J "NBH" 02_map.sh "NBH"')

png(paste0('~/public_html/',pop,'_RF_physpo.png'), width=2000, height=2000)
par(mfrow=c(4,6))
for(i in 1:24){
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross,i))))
 X <- 1:length(Y)
 plot(X,Y, xlab=paste('chr',i), ylab='physical position')
}
dev.off()

################################################################################
