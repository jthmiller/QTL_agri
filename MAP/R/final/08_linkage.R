#!/bin/R
pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- as.numeric(commandArgs(TRUE)[2])

################################################################################
##load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################
fl <- paste0(pop,'_imp.mapped.tsp.csv')
fl <- file.path(mpath,fl)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

## DOWNSAMPLED
fl <- file.path(paste0(pop,'_downsampled.csv'))
cross <- read.cross(file=fl , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross$pheno <- as.data.frame(cross$pheno)
################################################################################

ahr <- get_AHR(cross)

##rf <- subset(cross, chr = c(1:4,6:24))
rf <- est.rf(cross, maxit=1000, tol=1e-6)

#############################################
### test locus interaction seg distortion

probs <- c(0.0625,0.125,0.25)
gts <- c('AA','AB','BB')

homs <- c('AA','BB')
hets <- 'AB'

tr.table <- matrix(NA, ncol=3, nrow=3)
rownames(tr.table) <- colnames(tr.table) <- gts

tr.table[homs,homs] <- 0.0625
tr.table[hets,homs] <- 0.125
tr.table[homs,hets] <- 0.125
tr.table[hets,hets] <- 0.25

gtf <- c('AA','AB','BB')
gt_gt <- cbind(rep(gtf,3),c(rep('AA',3),rep('AB',3),rep('BB',3)))
gt_names <- paste0(gt_gt[,1],gt_gt[,2])
gt_probs <- setNames(tr.table[gt_gt], gt_names)

rf.gts <- pull.geno(rf)

csq <- function(mara, marb) {
 test <- factor(paste0(factor(mara, labels = gtf), factor(marb, labels = gtf)), levels = gt_names)
 chisq.test(table(test), p = gt_probs)$p.value
}

csq.each <- function(X){ apply(rf.gts, 2, csq, marb = X) }

library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
csq.pval  <- foreach(marb = iter(rf.gts, by='column'), .inorder = F, .packages = libs2load) %dopar% csq.each(marb)
colnames(csq.pval) <- rownames(csq.pval) <- colnames(rf.gts)

#### without parallel ################################

### applying over columns
#csq.pval <- apply(rf.gts, 2, function(X){
# apply(rf.gts, 2, csq, marb = X)
#})
#colnames(csq.pval) <- rownames(csq.pval) <- colnames(rf.gts)
#csq.bk <- csq.pval



########################################

TeaTasting <-
matrix(c(3, 1, 1, 3),
       nrow = 2,
       dimnames = list(Guess = c("Milk", "Tea"),
                       Truth = c("Milk", "Tea")))

fisher.test(TeaTasting, alternative = "greater")
