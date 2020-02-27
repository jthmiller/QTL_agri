#!/bin/R

################################################################################
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

################################################################################
fl <- "NBH_4822_imputed_NW_tsp.csv"
cross <- read.cross(file=fl , format = "csv", dir=mpath,
 genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

################################################################################

#################################################################################
png(paste0('~/public_html/',pop,'high_conf_imputed_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
  plotRF(cross5,chr=B)
}
dev.off()
################################################################################
################################################################################
plot_test('nbh_high_confid_rf')
 plotRF(cross5,zmax=10,col.scheme="redblue")
dev.off()
################################################################################
################################################################################

### ESTIMATE MAP ###############################################################

loglik <- err <- c(0.0001, 0.0005, 0.001, 0.01, 0.05)

cross.sub <- subset(cross,chr=c(1,3,4,7,9,12,15,17,20,23))

update.lik <- function(z){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross.sub,maxit=1000, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

library(qtl)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z)

loglik <- unlist(loglik)

lod <- (loglik - max(loglik))/log(10)

erprob <- err[which.max(lod)]

print(paste('error lod =',erprob))

cross_map <-  est.map(cross, error.prob=erprob, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

direc <- sapply(1:24,function(i) {
 pos <- as.numeric(gsub(".*:","",markernames(cross,i)))
 map <- as.numeric(pull.map(cross)[[i]])
 cor(pos,map, use="complete.obs")
})

if(any(direc < 0)) cross <- flip.order(cross,which(direc < 0))

mapfile <- paste0(pop,'_',sum(nmar(cross)),'_imputed_high_confidence_tsp_mapped')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

print(paste(filename, 'cross written'))
