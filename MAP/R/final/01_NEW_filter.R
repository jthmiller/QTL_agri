#!/bin/R

pop <- 'NEW'

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
mfl <- file.path(mpath,'NBH_markernames.tsv')
nbh_marks <- read.table(mfl)
cross <- pull.markers(cross,nbh_marks$x)
################################################################################

################################################################################
### Switch phase and keep only parent conf markers #############################
### ENRICH FOR AAxBB ##########################################################

#### DROP DANGEROUS ABxAB cross ##################################################
##DROP1 <- pull.geno(cross)[cross$pheno$ID=="NEW_NEW911M",]
##DROP1 <- names(DROP1)[which(as.numeric(DROP1)==2)]
##DROP2 <- pull.geno(cross)[cross$pheno$ID=="NEW_NEW911F",]
##DROP2 <- names(DROP2)[which(as.numeric(DROP2)==2)]
##DROP <- intersect(DROP1,DROP2)
##cross <- drop.markers(cross,DROP)
##################################################################################
##
##### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
##bfix <- pull.geno(cross)[cross$pheno$ID=="NEW_NEW911M",]
##bfix_swit1 <- names(bfix)[which(as.numeric(bfix)==1)]
##bfix <- pull.geno(cross)[cross$pheno$ID=="NEW_NEW911F",]
##bfix_swit2 <- names(bfix)[which(as.numeric(bfix)==3)]
##bfix_swit12 <- intersect(bfix_swit1 ,bfix_swit2)
##
##cross <- switchAlleles(cross, markers = bfix_swit12)
##################################################################################
##
####ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))
####cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
#### maxit = 1, tol = 0.1, sex.sp = F)
####cross <- calc.errorlod(cross, err=0.05)
####ca <- checkAlleles(cross,threshold=2)
####cross <- switchAlleles(cross,ca$marker)
##
##################################################################################
##### Get highly likely AB x AB markers ##########################################
##bfix1 <- pull.geno(cross)[cross$pheno$ID=="NEW_NEW911M",]
##bfix1 <- names(bfix1)[which(as.numeric(bfix1)==3)]
##bfix2 <- pull.geno(cross)[cross$pheno$ID=="NEW_NEW911F",]
##bfix2 <- names(bfix2)[which(as.numeric(bfix2)==1)]
##parABxAB <- intersect(bfix1,bfix2)
##
##gt_nopar <- geno.table(subset(cross,ind=!cross$pheno$ID %in% c("NEW_NEW911M","NEW_NEW911F")))
##parABxAB <- intersect(rownames(gt_nopar[which(gt_nopar$P.value > 0.01),]) ,parABxAB)
##cross.1 <- pull.markers(cross,parABxAB)
##################################################################################

##### TEST SAMPLE GT SIMILARITY ##################################################
##cross.1 <- subset(cross.1,ind=!cross.1$pheno$ID%in%c("NEW_NEW911M","NEW_NEW911F"))
##cpgt <- comparegeno(cross.1)
##colnames(cpgt) <- cross.1$pheno$ID
##rownames(cpgt) <- cross.1$pheno$ID
##cpgt[cpgt==NaN] <- NA
##diag(cpgt) <- NA
##cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
##################################################################################
##png(paste0('~/public_html/NBH_relat.png'))
## hist(cpgt)
##dev.off()
################################################################################
toss_missing <- c('NEW_11365','NEW_11398','NEW_11403','NEW_11388','NEW_11228',' NEW_11298','NEW_11407','NEW_11370')
################################################################################

################################################################################
#### Pvalue and Missing ##############################################
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss_missing,"NEW_NEW911M","NEW_NEW911F")))
bfixA <- rownames(gt[which(gt$P.value > 0.00001),])

################################################################################
png(paste0('~/public_html/NEW_pval.png'))
 hist(-log10(gt$P.value),breaks=30)
dev.off()
################################################################################

###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss_missing,"NEW_NEW911M","NEW_NEW911F"))
################################################################################
gt <- geno.table(cross)

sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

sm <- scanone(cross, pheno.col=4, model="binary",method="mr")

plot_test('new_mar_regression', width = 1500, height = 750)
par(mfrow=c(2,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,10), cex = 0.25)
 plot(1:length(gt[bfixA,1]), -log10(gt[bfixA,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
dev.off()


##for(Z in 1:24){
## reorg.1 <- formLinkageGroups(subset(cross,chr=Z), max.rf = 0.15, min.lod = 12, reorgMarkers = TRUE)
## swits <- markernames(reorg.1, chr=2)
## reorg.1 <- switchAlleles(reorg.1, markers = markernames(reorg.1,chr=2))
## reorg.2 <- formLinkageGroups(reorg.1, max.rf = 0.15, min.lod = 12, reorgMarkers = TRUE)
## r2gt <- geno.table(reorg.2,chr=2)
##
## if(mean(r2gt$P.value) > 0.0001 & length(r2gt$P.value) > 10 ){
##   swits2 <- markernames(reorg.2, chr=2)
##   reorg.2 <- switchAlleles(reorg.2, markers = markernames(reorg.2,chr=2))
##   reorg.2 <- formLinkageGroups(reorg.2, max.rf = 0.15, min.lod = 12, reorgMarkers = TRUE)
##   swits <<- c(swits,swits2[swits2 %in% markernames(reorg.2, chr=1)])
##  }
##
## subs <- markernames(reorg.2, chr=1)
## drops <- markernames(reorg.1)[!markernames(reorg.1) %in% subs]
## cross <<- switchAlleles(cross, swits)
## cross.par <<- switchAlleles(cross.par, swits)
## cross <<- drop.markers(cross, drops)
##}

fl <- file.path(mpath,'NEW_unmapped_filtered')
write.cross(cross,filestem=fl,format="csv")
