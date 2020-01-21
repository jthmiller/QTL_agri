#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtlDesign')
library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################

print(paste(cores,'cores'))
erp <-

################################################################################

################################################################################
load(file.path(mpath,paste0(pop,'_downsampled.rsave')))

map10 <- pull.map(cross)
n.sim <- 10000
res0 <- rep(NA, n.sim)
for(i in 1:n.sim) {
x <- sim.cross(map10[c(1:4,6:24)], n.ind=91, type="f2", error.prob=0.0025)
x <- calc.genoprob(x, step=1)
out <- scanone(x, method="hk")
res0[i] <- max(out[,3])
print(i)
}

print(thr <- quantile(res0, 0.95))
print(G <- sum(summary(map10)[c(1:4,6:24),"length"]))

# d is marker spacing
thresh(G, "f2", d=1, p=0.05)
## 4.053364

plot_test('no_qtl_lod_dist_NBH')
hist(res0, breaks=100, xlab="Genome-wide maximum LOD score")
rug(res0)
dev.off()

## simulations with QTL explaining 25% of the pheno at 54 cms on chr 1
alpha <- sqrt(2*0.25/(1-0.08))
n.sim <- 10000
loda <- est <- lo <- hi <- rep(NA, n.sim)
for(i in 1:n.sim) {
 x <- sim.cross(map10[1], n.ind=91, type="f2", model=c(1, 54, alpha, 0))
 x <- calc.genoprob(x, step=1)
 out <- scanone(x, method="hk")
 loda[i] <- max(out[,3])
 temp <- out[out[,3]==loda[i],2]
if(length(temp) > 1) temp <- sample(temp, 1)
 est[i] <- temp
 li <- lodint(out)
 lo[i] <- li[1,2]
 hi[i] <- li[nrow(li),2]
 print(i)
}

mean(loda >= thr)

#######
# Distribution of the chromosome-wide maximum LOD scores in the presence of a
# single QTL responsible for 25% of the phenotypic variance, for the case of an
# intercross with 91 individuals and with equally spaced markers at a 10 cM spacing.
plot_test('single_qtl_lod_dist_NBH')
hist(loda, breaks=100, xlab="Maximum LOD score")
dev.off()
#######
plot_test('single_qtl_lod_location_NBH')
hist(est, breaks=100, xlab="Estimated QTL location (cM)")
rug(map10[[1]])
dev.off()

save.image(file.path(mpath,paste0(pop,'_powercalc_NBH.rsave')))

##Note that we use theta=0.09, the approximate recombination fraction between two markers (from the Haldane map function, for a genetic distance of 10 cM).
## change to 0.009
powercalc("f2", 91, sigma2=1, effect=c(alpha,0), thresh=thr,theta=0.009)
