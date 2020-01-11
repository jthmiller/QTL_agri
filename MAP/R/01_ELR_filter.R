#!/bin/R

pop <- 'ELR'

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

### Switch phase and keep only parent conf markers##############################
### is 10869 the real ELR mother?
### ENRICH FOR AAxBB
##cross.bk <- cross
## DROP DANGEROUS ABxAB cross
#DROP <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
#DROP <- names(DROP)[which(as.numeric(DROP)==2)]
#cross <- drop.markers(cross,DROP)
################################################################################

### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
bfix <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
bfix_swit <- names(bfix)[which(as.numeric(bfix)==1)]
gt <- geno.table(cross)
bfix_swit <- intersect(rownames(gt[which(gt$P.value > 0.05),]) ,bfix_swit)
cross <- switchAlleles(cross, markers = bfix_swit)
################################################################################

### Get highly likely AB x AB markers ##########################################
bfix <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
bfix <- names(bfix)[which(as.numeric(bfix)==3)]
parABxAB <- intersect(rownames(gt[which(gt$P.value > 0.05),]) ,bfix)
cross.1 <- pull.markers(cross,parABxAB)
################################################################################

### TEST SAMPLE GT SIMILARITY ##################################################
cross.1 <- subset(cross.1,ind=!cross$pheno$ID%in%c('BLI_BI1124M'))
cpgt <- comparegeno(cross.1)
colnames(cpgt) <- cross.1$pheno$ID
rownames(cpgt) <- cross.1$pheno$ID
cpgt[cpgt==NaN] <- NA
diag(cpgt) <- NA
cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
################################################################################

################################################################################
###### Remove the samples related by more than 80% of genotypes #####
wh <- which(cpgt > 0.7, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
mats <- cbind(rownames(wh),colnames(cpgt)[as.numeric(wh[,2])])
toss.missing <- apply(mats,1,function(X){
 X[which.max(c(nmissing(cross.1)[X[1]],nmissing(cross.1)[X[2]]))]
})

## USE MIS_ID'd samples for map, but not QTL
toss.missing <- c(mats[,1],"ELR_10869")
################################################################################

################################################################################
#### FILTER BY Pvalue and Missing ##############################################
cross.1 <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'BLI_BI1124M','ELR_ER1124F'))
gt <- geno.table(cross.1)
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 5),])
################################################################################
## cross.bk <- cross
###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'BLI_BI1124M','ELR_ER1124F'))
################################################################################
LOD <- 17
RF <- 0.15

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

fl <- file.path(mpath,paste0(pop,'_unmapped_filtered'))
write.cross(cross,filestem=fl,format="csv")
