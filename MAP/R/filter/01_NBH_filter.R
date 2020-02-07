#!/bin/R

pop <- 'NBH'

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

### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
bfix <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
bfix_swit1 <- names(bfix)[which(as.numeric(bfix)==1)]
bfix <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
bfix_swit2 <- names(bfix)[which(as.numeric(bfix)==3)]
bfix_swit12 <- intersect(bfix_swit1 ,bfix_swit2)

cross <- switchAlleles(cross, markers = bfix_swit12)
################################################################################

################################################################################
### Get highly likely AB x AB markers ##########################################
bfix1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
bfix1 <- names(bfix1)[which(as.numeric(bfix1)==3)]
bfix2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
bfix2 <- names(bfix2)[which(as.numeric(bfix2)==1)]
parABxAB <- intersect(bfix1,bfix2)

gt_nopar <- geno.table(subset(cross,ind=!cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F')))
parABxAB <- intersect(rownames(gt_nopar[which(gt_nopar$P.value > 0.01),]) ,parABxAB)
cross.1 <- pull.markers(cross,parABxAB)
################################################################################

### TEST SAMPLE GT SIMILARITY ##################################################
cross.1 <- subset(cross.1,ind=!cross.1$pheno$ID%in%c('NBH_NBH1M','NBH_NBH1F'))
cpgt <- comparegeno(cross.1)
colnames(cpgt) <- cross.1$pheno$ID
rownames(cpgt) <- cross.1$pheno$ID
cpgt[cpgt==NaN] <- NA
diag(cpgt) <- NA
cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
################################################################################
png(paste0('~/public_html/NBH_relat.png'))
 hist(cpgt)
dev.off()
################################################################################
toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137","NBH_5646")
################################################################################

################################################################################
#### Pvalue and Missing ##############################################
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F')))
bfixA <- rownames(gt[which(gt$P.value > 0.00001 & gt$missing < 4),])
##bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 5),])

################################################################################
## Determine what percent of markers are kept after filter
table(gsub(":.*","",bfixA))/table(gsub(":.*","",markernames(cross)))

################################################################################
## Get a cross object of parent genotypes
cross.par <- subset(cross, ind=cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
DROP1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
DROP1 <- names(DROP1)[which(as.numeric(DROP1)==2)]
DROP2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]

###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################

sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

sm <- scanone(cross, pheno.col=4, model="binary",method="mr")

plot_test('nbh_mar_regression', width = 1500, height = 750)
par(mfrow=c(2,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[bfixA,1]), gt[bfixA,'P.value'], pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
dev.off()

###### Retain markers that are linked ########
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

fl.par <- file.path(paste0(pop,'_parents_filtered'))
fl.par <- file.path(mpath,fl.par)
write.cross(cross.par,filestem=fl.par,format="csv")

system('sbatch -J "NBH" 02_map.sh "NBH"')

png(paste0('~/public_html/',pop,'_RF_physpo.png'))
par(mfrow=c(4,6))
for(i in 1:24){
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross,i))))
 X <- 1:length(Y)
 plot(X,Y, xlab=paste('chr',i), ylab='physical position')
}
dev.off()
