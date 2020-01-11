#!/bin/R

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

pars <- c('BRP_BRP1M','BRP_BRP8F','BRP_BRP1F','BRP_BRP8M')

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
bfixA <- rownames(gt[which(gt$P.value > 0.0001),])
##bfixA <- rownames(gt[which(gt$P.value > 0.00001 & gt$missing < 5),])
##bfixA <- rownames(gt[which(gt$P.value > 0.000001 & gt$missing < 5),])
################################################################################

###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,pars))
################################################################################

png(paste0('~/public_html/BRP_pvals.png'))
 hist(log10(gt$missing))
dev.off()

fl <- file.path(mpath,'BRP_unmapped_filtered')
write.cross(cross,filestem=fl,format="csv")
