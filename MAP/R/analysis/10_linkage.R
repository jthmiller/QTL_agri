#!/bin/R
mpath <- '/home/jmiller1/QTL_agri/data'
library('qtl')
library('snow')

cores <- as.numeric(commandArgs(TRUE)[3])
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

print(commandArgs(TRUE))

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

################################################################################
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
### HEATMAP WITH INTERACTION LOD AND TWO_LOCUS SEG DISTORTION
##load(file.path(mpath,paste0(pop,'_scan2_bin_mr.rsave')))
################################################################################

names(cross$geno) <- ifelse(names(cross$geno) == "X","5",names(cross$geno))
attr(cross$geno[["5"]], 'class') <- 'A'

#############################################
ahr_genes <- get_AHR(cross)
gt <- geno.table(cross)
ahr_genes$segdist <- -log10(gt[ahr_genes$close_marker,'P.value'])
ahr_genes_sub <- ahr_genes[!is.na(ahr_genes$PATH),]
#############################################

cross <- est.rf(cross, maxit=100, tol=1e-6)

#############################################
### test 2 locus interaction seg distortion
##rf <- subset(cross, chr = c(1:4,6:24))
rf <- cross

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

gtf_mod <- c('AA','AB','BB')
gt_names_mod <- c('AAAA','AABB','BBAA','BBBB')
gt_prob_mod <- setNames(c(0.25,0.25,0.25,0.25),gt_names_mod)

rf.gts <- pull.geno(rf)

################################################################################
################################################################################

csq <- function(mara, marb) {
 test <- factor(paste0(factor(mara, labels = gtf), factor(marb, labels = gtf)), levels = gt_names)
 chisq.test(table(test), p = gt_probs)$p.value
}

csq.each <- function(X){ apply(rf.gts, 2, csq, marb = X) }

################################################################################
################################################################################

csq_mod <- function(mara, marb) {
 test <- factor(paste0(factor(mara, labels = gtf_mod), factor(marb, labels = gtf_mod)), levels = gt_names_mod)
 chisq.test(table(test), p = gt_prob_mod)$p.value
}

csq_mod.each <- function(X){ apply(rf.gts, 2, csq_mod, marb = X) }
################################################################################
################################################################################

### WITH PARALLEL #########################################
#cores <- as.numeric(commandArgs(TRUE)[2])
cores <- 20
library(doParallel)
cl <- makeCluster(cores)
registerDoParallel(cl)
csq.pval  <- foreach(marb = iter(rf.gts, by='column'), .inorder = F, .packages = libs2load) %dopar% csq.each(marb)
csq.pval <- do.call(rbind,csq.pval)
colnames(csq.pval) <- rownames(csq.pval) <- colnames(rf.gts)
csq.pval.bk <- csq.pval
#############################################################

### HOMOZYGOTE DISTORTION
cores <- 20
library(doParallel)
cl <- makeCluster(cores)
registerDoParallel(cl)
csq_mod.pval  <- foreach(marb = iter(rf.gts, by='column'), .inorder = F, .packages = libs2load) %dopar% csq_mod.each(marb)
csq_mod.pval <- do.call(rbind,csq_mod.pval)
colnames(csq_mod.pval) <- rownames(csq_mod.pval) <- colnames(rf.gts)
csq_mod.pval.bk <- csq_mod.pval

################################################################################
save.image(file.path(mpath,paste0(pop,'_csq_scan.rsave')))
################################################################################

rf.plots <- rf

## Set within chromosomes to zero #####
for (i in chrnames(cross)){
 mars <- markernames(rf, i)
 csq.pval[mars,mars] <- NA
 rf.plots$rf[mars,mars] <- NA
}

## Set within chromosomes to zero #####
for (i in chrnames(cross)){
 mars <- markernames(rf, i)
 csq_mod.pval[mars,mars] <- NA
}

mf1 <- file.path(mpath,paste0(pop,'_csq.pval.tsv'))
write.table(csq.pval,mf1)

mf2 <- file.path(mpath,paste0(pop,'_csq_mod.pval.tsv'))
write.table(csq_mod.pval,mf2)

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
######################################################################

### TABLE OF THE HIGHEST HOMOZYGOTE LOD SCORES OF LINKAGE FOR EACH CHR
maxdist_mod <- lapply(chrnames(cross)[-17], function(i) {
  mars <- markernames(rf, i)
  a <- which(csq_mod.pval[mars,] == min(csq_mod.pval[mars,], na.rm = T), arr.ind=T)
  b <- markernames(rf)[a[,'col']][1]
  a <- rownames(a)[1]
  cbind(a, b, -log10(csq_mod.pval[cbind(a,b)]))
})
maxdist_mod <- do.call(rbind,maxdist_mod)
maxdist_mod <- maxdist_mod[order(as.numeric(maxdist_mod[,3])),]
maxdist_mod <- data.frame(maxdist_mod, stringsAsFactors = F)
######################################################################
