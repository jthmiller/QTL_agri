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
## HIGH CONFID IMPUTED
mapfile <- paste0(pop,'_',sum(nmar(cross10)),'_imputed_high_confidence_tsp')
filename <- file.path(mpath,mapfile)
cross <- read.cross(
 file = filename,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)


## DOWNSAMPLED
#fl <- file.path(paste0(pop,'_downsampled.csv'))
#cross <- read.cross(file=fl , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#cross$pheno <- as.data.frame(cross$pheno)
################################################################################

cross <- cross11

ahr_genes <- get_AHR(cross)
gt <- geno.table(cross)
ahr_genes$segdist <- -log10(gt[ahr_genes$close_marker,'P.value'])
ahr_genes_sub <- ahr_genes[!is.na(ahr_genes$PATH),]




#############################################
### test 2 locus interaction seg distortion
##rf <- subset(cross, chr = c(1:4,6:24))
rf <- est.rf(cross, maxit=1000, tol=1e-6)

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

### WITH PARALELLE #########################################
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
csq.pval  <- foreach(marb = iter(rf.gts, by='column'), .inorder = F, .packages = libs2load) %dopar% csq.each(marb)
csq.pval <- do.call(rbind,csq.pval)
colnames(csq.pval) <- rownames(csq.pval) <- colnames(rf.gts)
csq.pval.bk <- csq.pval
#############################################################

rf.plots <- rf

## Set within chromosomes to zero #####
for (i in chrnames(cross)){
 mars <- markernames(rf, i)
 csq.pval[mars,mars] <- NA
 rf.plots$rf[mars,mars] <- NA
}

########################################
### TABLE OF THE HIGHEST LOD SCORES OF LINKAGE FOR EACH CHR
maxdist <- lapply(chrnames(cross), function(i) {
  mars <- markernames(rf, i)
  a <- which(csq.pval[mars,] == min(csq.pval[mars,], na.rm = T), arr.ind=T)
  b <- markernames(rf)[a[,'col']][1]
  a <- rownames(a)[1]
  cbind(a, b, -log10(csq.pval[cbind(a,b)]))
})
maxdist <- do.call(rbind,maxdist)
maxdist <- maxdist[order(as.numeric(maxdist[,3])),]
maxdist <- data.frame(maxdist, stringsAsFactors = F)

##rownames(maxdist)  <- as.character(unique(bin.em.2$map$chr))

plot_test('CHRxCHR_LOD_scores',height=1000,width=1000)
plotRF(rf.plots,zmax=8,col.scheme="redblue")
dev.off()


###############

csq.pval[ahr_genes_sub$close_marker,ahr_genes_sub$close_marker]


crs <- pull.markers(rf.plots, ahr_genes$close_marker)

plot_test('qtlxqtl_SegLOD_scores')
plotRF(crs,zmax=7,col.scheme="redblue")
dev.off()
