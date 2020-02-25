#!/bin/R

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

##source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'
#mapfile <- paste0(pop,'_unmapped_all_mark_imputed_',i,'_tsp.csv')
#filename <- file.path(mpath,mapfile)
#mapfile <- paste0(pop,'_order_impute_NW_',i,'_tsp.csv')
#mapfile <- paste0(pop,'_order_impute_',i,'_tsp.csv')
#mapfile <- paste0(pop,'_order_impute_',i,'_tsp.csv')
#filename <- file.path(mpath,mapfile)

mapfile <- paste0(pop,'_',sum(nmar(cross10)),'_imputed_high_confidence_tsp')
filename <- file.path(mpath,mapfile)

#mapfile <- paste0(pop,'_imputed_',i,'_tsp.csv')
#filename <- file.path(mpath,mapfile)

#libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
#suppressMessages(sapply(libs2load, require, character.only = TRUE))

library(qtl)
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

################################################################################
fl <- paste0(pop,'_imp.mapped.tsp.csv')
fl <- file.path(mpath,mapfile)

cross2 <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

################################################################################

loglik <- err <- c(0.0001, 0.001, 0.01, 0.05)

update.lik <- function(z){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross,maxit=100, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

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

if(any(direc < 0) cross <- flip.order(cross,which(direc < 0))

mapfile <- paste0(pop,'_',sum(nmar(cross)),'_imputed_high_confidence_tsp_mapped')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

print(paste(filename, 'cross written'))
################################################################################
