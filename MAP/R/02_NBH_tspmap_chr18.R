#!/bin/R

i <- 18
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
toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137")
################################################################################

################################################################################
#### Pvalue and Missing ##############################################
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F')))
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 5),])
################################################################################

cross.par <- subset(cross, ind=cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
DROP1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
DROP1 <- names(DROP1)[which(as.numeric(DROP1)==2)]
DROP2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]

###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################









check_crosslinks <- function(X,Y){
 reorg <- subset(cross,chr=c(Y,X))
 formLinkageGroups(reorg, max.rf = 0.1, min.lod = 10, reorgMarkers = TRUE)
}


check_18 <- lapply(c(1:17,19:24),check_crosslinks,Y=18)

################################################################################
reorg <- formLinkageGroups(cross, max.rf = 0.1, min.lod = 10, reorgMarkers = TRUE)
fl <- file.path(mpath,'NBH_unmapped_reassigned_markers')
write.cross(reorg,filestem=fl,format="csv")
################################################################################


sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

################################################################################

ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))
cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)

#### CHR8 #######



cross <- subset(cross,ind=!cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))


png(paste0('~/public_html/NBH_gts_preclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross, chr=i, cex=2)
dev.off()

################################################################################
cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
cross <- removeDoubleXO(cross, chr=i)
cross <- calc.errorlod(cross, err=0.05)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
################################################################################

################################################################################
gt <- geno.table(cross)
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 6),])
cross <- pull.markers(cross,bfixA)
################################################################################

################################################################################
cross <-tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

if(cor(pos,map) < 0){
 cross <<- flip.order(cross, i)
}

cross <- shiftmap(cross, offset=0)

cross_map <-  est.map(cross, error.prob=0.04,map.function="kosambi",maxit=100000,tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

write.cross(cross,chr=i,filestem=filename,format="csv")

png(paste0('~/public_html/NBH_RF_concord',i,'_tsp.png'))
  plotRF(cross)
dev.off()
