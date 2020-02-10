#!/bin/R

#pop <- 'ELR'
#
#source("/home/jmiller1/QTL_agri/MAP/control_file.R")
#mpath <- '/home/jmiller1/QTL_agri/data'
#
#################################################################################
### read in the QTL cross
#cross <- read.cross.jm(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
#format = "csvr", geno = c(1:3), estimate.map = FALSE)
#################################################################################
#
#### Pull names from plinkfile ##################################################
#path <- file.path(mpath, paste(pop, ".ped", sep = ""))
#popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
#indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
#cross$pheno$ID <- paste(popname, indname, sep = "_")
#################################################################################
#
##### PHENO #####################################################################
#cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
#cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
#################################################################################

##### Write raw table ###########################################################
#gts <- file.path(mpath, paste0(pop,'_gts.tsv'))
## gt <- geno.table(cross)
###write.table(gt,gts)
#gt <- read.table(gts)
#################################################################################

################################################################################
### Switch phase and keep only parent conf markers##############################
### is 10869 the real ELR mother?
##cross.bk <- cross
################################################################################

### DROP DANGEROUS ABxAB cross ##################################################
##DROP <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
##DROP <- names(DROP)[which(as.numeric(DROP)==2)]
#flt <- paste0(pop,'_ABxAB_markernames.tsv')
#flt <- file.path(mpath,flt)
##write.table(flt, DROP)
#DROP <- read.table(flt)
#cross <- drop.markers(cross,marks$x)
#################################################################################

#################################################################################
#### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
#bfix <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
#bfix_swit <- names(bfix)[which(as.numeric(bfix)==1)]
#bfix_swit <- intersect(rownames(gt[which(gt$P.value > 0.05),]) ,bfix_swit)
#swt <- file.path(mpath,paste0(pop,'_switch_markernames.tsv'))
##write.table(bfix_swit, swt)
#bfix_swit <- read.table(flt)
#cross <- switchAlleles(cross, markers = bfix_swit$x)
#################################################################################

####################################################################################
####### Get highly likely AB x AB markers ##########################################
##gt <- geno.table(cross)
##bfix <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
##bfix <- names(bfix)[which(as.numeric(bfix)==3)]
##parABxAB <- intersect(rownames(gt[which(gt$P.value > 0.05),]) ,bfix)
##cross.1 <- pull.markers(cross,parABxAB)
####################################################################################
####
####### TEST SAMPLE GT SIMILARITY ##################################################
##cross.1 <- subset(cross.1,ind=!cross$pheno$ID%in%c('BLI_BI1124M'))
##cpgt <- comparegeno(cross.1)
##colnames(cpgt) <- cross.1$pheno$ID
##rownames(cpgt) <- cross.1$pheno$ID
##cpgt[cpgt==NaN] <- NA
##diag(cpgt) <- NA
##cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
####################################################################################
####
####################################################################################
########## Remove the samples related by more than 80% of genotypes #####
##wh <- which(cpgt > 0.7, arr=TRUE)
##wh <- wh[wh[,1] < wh[,2],]
##mats <- cbind(rownames(wh),colnames(cpgt)[as.numeric(wh[,2])])
##toss.relat <- unique(apply(mats,1,function(X){
## X[which.max(c(nmissing(cross.1)[X[1]],nmissing(cross.1)[X[2]]))]
##}))
##
####USE MIS_IDd samples for map, but not QTL
##
######### toss related and then #################################################
####### FILTER BY Pvalue and Missing ##############################################
## toss.related <- c(toss.relat,"ELR_10987")
## cross.1 <- subset(cross, ind=!cross$pheno$ID %in% c(toss.related,'BLI_BI1124M','ELR_ER1124F'))
##
## gt.pmiss <- geno.table(cross.1)
## gtpm <- file.path(mpath,paste0(pop,'_gtpmiss.tsv'))
## write.table(gt.pmiss, gtpm)
## gt.pmiss <- read.table(gtpm)
#### bfixA <- rownames(gt.pmiss[which(gt.pmiss$P.value > 0.001 & gt.pmiss$missing < 5),])
##
## ind <- which(gt.pmiss$missing < 10)
## png('~/public_html/elr_mis.png')
## hist(gt.pmiss$missing[ind], pch=21,breaks=50)
## dev.off()
################################################################################
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
LOD <- as.numeric(commandArgs(TRUE)[3])
RF <- as.numeric(commandArgs(TRUE)[4])

#LOD <- 14
#RF <- 0.15

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

cross <- read.cross.jm(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)

path <- file.path(mpath, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))

par.genos <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]

sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

################################################################################

flt <- file.path(mpath,paste0(pop,'_ABxAB_markernames.tsv'))
DROP <- read.table(flt)$x

## Switch AAxBB
swt <- file.path(mpath,paste0(pop,'_switch_markernames.tsv'))
bfix_swit <- read.table(swt)$x

## Pval and missing
gtpm <- file.path(mpath,paste0(pop,'_gtpmiss.tsv'))
gt.pmiss <- read.table(gtpm)
bfixA <- rownames(gt.pmiss[which(gt.pmiss$P.value > 0.001 & gt.pmiss$missing < 5),])

## Bad data individuals
toss.related <- c("ELR_10978","ELR_10977","ELR_10982","ELR_10974","ELR_10980","ELR_10973","ELR_10971","ELR_10979","ELR_10987")
##toss.badata <- c("ELR_10869","ELR_10967","ELR_11592","ELR_11115","ELR_11103","ELR_10981","ELR_11593")
toss.badata <- c("ELR_10869","ELR_10987")

cross <- drop.markers(cross,DROP)
cross <- switchAlleles(cross, markers = bfix_swit)
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.related,toss.badata,'BLI_BI1124M','ELR_ER1124F'))
################################################################################

sm <- scanone(cross, pheno.col=4, model="binary",method="mr")

plot_test('elr_mar_regression', width = 1500, height = 750)
par(mfrow=c(2,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,8), cex = 0.25)
 plot(1:length(gt[bfixA,1]), -log10(gt[bfixA,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,5), cex = 0.25)
dev.off()

################################################################################

for(Z in 1:24){

 all <- subset(cross,chr=Z)
 reorg.2 <- formLinkageGroups(all, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 reorg.2a <- reorg.2

  ## switch it
 swits <- markernames(reorg.2a, chr=1)
 reorg.2a <- switchAlleles(reorg.2a, markers = swits)
 reorg.2a <- formLinkageGroups(reorg.2a, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

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

 cross <<- drop.markers(cross, drops)
 cross <<- switchAlleles(cross, switched)

}

################################################################################
## write #######################################################################
################################################################################
fl <- file.path(mpath,paste0(pop,'_unmapped_filtered'))
write.cross(cross,filestem=fl,format="csv")
################################################################################

system('sbatch 02_map.sh "ELR"')

png(paste0('~/public_html/',pop,'_RF_physpo.png'),width=2000, height=2000)
par(mfrow=c(4,6))
for(i in 1:24){
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross,i))))/1000000
 X <- 1:length(Y)
 plot(X,Y, xlab=paste('chr',i), ylab='physical position')
}
dev.off()

#cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

#################################################################################
#inds <- which(crossbk$pheno$ID %in% cross$pheno$ID)
#old.gt <- as.matrix(pull.geno(crossbk))
#new.gt <- as.matrix(pull.geno(cross))
#old.gt[old.gt == 2] <- NA
#new.gt[new.gt == 2] <- NA
#marks <- intersect(colnames(old.gt), colnames(new.gt))
#old.gt <- old.gt[inds,marks]
#new.gt <- new.gt[,marks]
#
#switched <- colnames(old.gt)[which(colSums(old.gt == new.gt, na.rm =T) == 0)]
#
#crossbk <- switchAlleles(crossbk, switched)
#BLI_BI1124M <- pull.geno(crossbk)[which(crossbk$pheno$ID == 'BLI_BI1124M'),marks]
#
#new_gts <- as.matrix(reorg.2a$geno[['1']]$data[,added])
#orig_gts <- as.matrix(all$geno[[as.character(Z)]]$data[,added])
#
#new_gts[new_gts == 2] <- NA
#orig_gts[orig_gts == 2] <- NA
#
#switched <- colnames(new_gts)[which(colSums(new_gts == orig_gts, na.rm =T) == 0)]
