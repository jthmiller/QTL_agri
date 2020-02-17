#!/bin/R

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

##source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
mapfile <- paste0(pop,'_unmapped_all_mark_imputed_',i,'_tsp.csv')
filename <- file.path(mpath,mapfile)

#libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
#suppressMessages(sapply(libs2load, require, character.only = TRUE))

library(qtl)
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

################################################################################

cross <- read.cross(file=mapfile , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

loglik <- err <- c(0.0001, 0.001, 0.01, 0.08)

update.lik <- function(z){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross,maxit=10000, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z)

loglik <- unlist(loglik)

lod <- (loglik - max(loglik))/log(10)

erprob <- err[which.max(lod)]

print(paste('error lod =',erprob))

cross_map <-  est.map(cross, error.prob=erprob, map.function="kosambi",maxit=10000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

write.cross(cross,chr=i,filestem=filename,format="csv")

print(paste(pop, 'cross written'))
################################################################################

Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)

png(paste0('~/public_html/',pop,'_RF_physpo_concord',i,'_tsp.png'),width=1000,height=500)
par(mfrow=c(1,3))
 plotRF(cross,main=NULL)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
 points(X,Y)
 plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))
dev.off()
