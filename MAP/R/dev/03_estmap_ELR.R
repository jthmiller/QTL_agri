#!/bin/R

pop <- 'ELR'
##source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
library(qtl)
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

################################################################################
fl <- paste0(pop,'_imp.mapped.tsp.csv')
fl <- file.path(mpath,fl)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
################################################################################

###############################################################################
## REMOVE MARKERS WITH HIGH MISSING DATA
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.125)
bfixA <- names(which(colSums(is.na(pull.geno(cross))) > mis))
print(paste('dropped',length(bfixA),'markers due to missing data'))
cross <- drop.markers(cross, bfixA)
print(sum(nmar(cross)))
################################################################################

################################################################################
loglik <- err <- c(0.0001, 0.001, 0.01, 0.05)

update.lik <- function(z){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross,maxit=100, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z)

loglik <- unlist(loglik)

lod <- (loglik - max(loglik))/log(10)

erprob <- err[which.max(lod)]

print(paste('error lod =',erprob))

##erprob <- 1e-04
par.est.map <- function(X) { est.map(cross, chr = X, error.prob=erprob, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE) }

### MAP ########################################################################
cl <- makeCluster(20)
registerDoParallel(cl)
maps <- foreach(X = 1:24, .inorder = T, .packages = libs2load) %dopar% par.est.map(X)
################################################################################

maps <- lapply(maps,"[[",1)
names(maps) <- c(1:24)
attr(maps, 'class') <- 'map'

cross <- qtl:::replace.map(cross,maps)

direc <- sapply(1:24,function(i) {
 pos <- as.numeric(gsub(".*:","",markernames(cross,i)))
 map <- as.numeric(pull.map(cross)[[i]])
 cor(pos,map, use="complete.obs")
})

if(any(direc < 0)) cross <- flip.order(cross,which(direc < 0))

mapfile <- paste0(pop,'_',sum(nmar(cross)),'_imputed_tsp_mapped')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

print(paste(filename, 'cross written'))

cross <- est.rf(cross)

################################################################################
save.image(file.path(mpath,paste0(pop,'_estrf_imputed.rsave')))
################################################################################

################################################################################
png(paste0('~/public_html/',pop,'high_conf_imputed_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
  plotRF(cross,chr=B)
}
dev.off()
################################################################################

################################################################################
plot_test(paste0(pop,'_high_confid_rf'))
 plotRF(cross,zmax=10,col.scheme="redblue")
dev.off()
################################################################################
################################################################################
