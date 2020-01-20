#!/bin/R

#pop <- 'ELR'
#
#source("/home/jmiller1/QTL_agri/MAP/control_file.R")
#mpath <- '/home/jmiller1/QTL_agri/data'
#

pop <- 'ELR'
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

## NEW GENOs
fl <- file.path(paste0(pop,'_unmapped_filtered.csv'))
cross_new <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
new.inds <- cross_new$pheno$ID
new.marks <- markernames(cross_new)


## Old genos
cross <- read.cross.jm(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
path <- file.path(mpath, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")

cross <- pull.markers(cross,new.marks)
cross <- subset(cross,ind=cross$pheno$ID %in% c(as.character(cross_new$pheno$ID),'BLI_BI1124M'))

## old parent genos
par.genos <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
################################################################################

inds <- intersect(cross$pheno$ID,cross_new$pheno$ID)
old.gt <- as.matrix(pull.geno(cross))
new.gt <- as.matrix(pull.geno(cross_new))
rownames(old.gt) <- cross$pheno$ID
rownames(new.gt) <- cross_new$pheno$ID
old.gt[old.gt == 2] <- NA
new.gt[new.gt == 2] <- NA
marks <- intersect(colnames(old.gt), colnames(new.gt))
old.gt <- old.gt[inds,marks]
new.gt <- new.gt[inds,marks]
switched <- colnames(old.gt)[which(colSums(old.gt == new.gt, na.rm =T) == 0)]
cross <- switchAlleles(cross, switched)
BLI_BI1124M <- pull.geno(cross)[which(cross$pheno$ID == 'BLI_BI1124M'),marks]


################################################################################

chrs <- gsub(":.*","",names(BLI_BI1124M))
pos <- as.numeric(gsub(".*:","",names(BLI_BI1124M)))
BLI_BI1124M[is.na(BLI_BI1124M)] <- 0
colrs <- c('grey','blue','purple','red')
names(colrs) <- c('0','1','2','3')

png("/home/jmiller1/public_html/ELR_par_BLI_BI1124M_genos.png", width=1500, height=750)
par(mfrow=c(4,6))

for (i in 1:24){

 ind <- which(chrs == i)
 X <- pos[ind]
 Y <- as.numeric(BLI_BI1124M)[ind]

 plot(X,Y,col='red', ylim=c(-0.5,3.5),xlim=c(0,max(X)), xlab=paste('chr',i), ylab='genotype', cex.axis=2,pch=19,cex.main=2)


 }
dev.off()

################################################################################
################################################################################
################################################################################

pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

## NEW GENOs
fl <- file.path(paste0(pop,'_unmapped_filtered.csv'))
cross_new <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
new.inds <- cross_new$pheno$ID
new.marks <- markernames(cross_new)

## Old genos
cross <- read.cross.jm(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
path <- file.path(mpath, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")

cross <- pull.markers(cross,new.marks)
cross <- subset(cross,ind=cross$pheno$ID %in% c(as.character(cross_new$pheno$ID),'NBH_NBH1F','NBH_NBH1M'))
################################################################################

inds <- intersect(cross$pheno$ID,cross_new$pheno$ID)
old.gt <- as.matrix(pull.geno(cross))
new.gt <- as.matrix(pull.geno(cross_new))
rownames(old.gt) <- cross$pheno$ID
rownames(new.gt) <- cross_new$pheno$ID
old.gt[old.gt == 2] <- NA
new.gt[new.gt == 2] <- NA
marks <- intersect(colnames(old.gt), colnames(new.gt))
old.gt <- old.gt[inds,marks]
new.gt <- new.gt[inds,marks]
switched <- colnames(old.gt)[which(colSums(old.gt == new.gt, na.rm =T) == 0)]
cross <- switchAlleles(cross, switched)
NBH_NBH1F <- pull.geno(cross)[which(cross$pheno$ID == 'NBH_NBH1F'),marks]
NBH_NBH1M <- pull.geno(cross)[which(cross$pheno$ID == 'NBH_NBH1M'),marks]
NBH_NBH1F[is.na(NBH_NBH1F)] <- 0
NBH_NBH1M[is.na(NBH_NBH1M)] <- 0


which(NBH_NBH1M == 2 & NBH_NBH1F == 2)
################################################################################
#colrs <- c('grey','blue','purple','red')
#names(colrs) <- c('0','1','2','3')


chrs <- gsub(":.*","",names(NBH_NBH1F))
pos <- as.numeric(gsub(".*:","",names(NBH_NBH1F)))


png("/home/jmiller1/public_html/NBH_par_genos.png", width=1500, height=750)
par(mfrow=c(4,6))

for (i in 1:24){

 ind <- which(chrs == i)
 X <- pos[ind]
 Yf <- as.numeric(NBH_NBH1F)[ind]
 Ym <- as.numeric(NBH_NBH1M)[ind]

 plot(X,Yf-0.1, col='blue', ylim=c(-0.5,3.5),xlim=c(0,max(X)), xlab=paste('chr',i), ylab='genotype', cex.axis=2,pch=19,cex.main=2)
 points(X,Ym+0.1,col='red',pch=19)

 }
dev.off()

################################################################################

chrs <- gsub(":.*","",names(BLI_BI1124M))
pos <- as.numeric(gsub(".*:","",names(BLI_BI1124M)))
BLI_BI1124M[is.na(BLI_BI1124M)] <- 0
colrs <- c('grey','blue','purple','red')
names(colrs) <- c('0','1','2','3')

png("/home/jmiller1/public_html/NBH_par_NBH1M_genos.png", width=1500, height=750)
par(mfrow=c(4,6))

for (i in 1:24){

 ind <- which(chrs == i)
 X <- pos[ind]
 Y <- as.numeric(BLI_BI1124M)[ind]

 plot(X,Y,col=colrs[as.character(BLI_BI1124M)[ind]], ylim=c(-0.5,3.5),xlim=c(0,max(X)), xlab=paste('chr',i), ylab='genotype', cex.axis=2,pch=19,cex.main=2)


 }
dev.off()
################################################################################
