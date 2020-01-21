#!/bin/R


map10 <- pull.map(cross)
n.sim <- 1000
res0 <- rep(NA, n.sim)
for(i in 1:n.sim) {
x <- sim.cross(map10[c(1:4,6:24)], n.ind=91, type="f2")
x <- calc.genoprob(x, step=1)
out <- scanone(x, method="hk")
res0[i] <- max(out[,3])
print(i)
}

print(thr <- quantile(res0, 0.95))
print(G <- sum(summary(map10)[c(1:4,6:24),"length"]))

# d is marker spacing
thresh(G, "f2", d=1, p=0.05)

plot_test('no_qtl_lod_dist_NBH')
hist(res0, breaks=100, xlab="Genome-wide maximum LOD score")
rug(res0)
dev.off()

## simulations with QTL explaining 25% of the pheno at 54 cms on chr 1
alpha <- sqrt(2*0.25/(1-0.08))
n.sim <- 1000
loda <- est <- lo <- hi <- rep(NA, n.sim)
for(i in 1:n.sim) {
 x <- sim.cross(map10[1], n.ind=91, type="f2", model=c(1, 54, alpha, 0))
 x <- calc.genoprob(x, step=1)
 out <- scanone(x, method="hk")
 loda[i] <- max(out[,3])
 temp <- out[out[,3]==loda[i],2]if(length(temp) > 1)
 temp <- sample(temp, 1)
 est[i] <- temp
 li <- lodint(out)
 lo[i] <- li[1,2]
 hi[i] <- li[nrow(li),2]
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

save.image(file.path(mpath,paste0(pop,'_powercalc_NBH.rsave')))

powercalc("f2", 91, sigma2=1, effect=c(alpha,0), thresh=thr,theta=0.09)


hist(est, breaks=100, xlab="Estimated QTL location (cM)")
rug(map10[[1]])
