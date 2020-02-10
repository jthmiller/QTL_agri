#!/bin/R

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

source("/home/jmiller1/QTL_agri/MAP/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'
fl <- file.path(paste0(pop,'_unmapped_filtered.csv'))
mapfile <- paste0(pop,'_all_mark_',i,'_tsp')
filename <- file.path(mpath,mapfile)
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

print(pop)
print(i)
################################################################################

cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross <- subset(cross,chr=i)

sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

################################################################################

ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))

cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)

png(paste0('~/public_html/ELR_gts_preclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross, chr=i, cex=2)
dev.off()

## 10932
################################################################################

#cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.1)
cross <- removeDoubleXO(cross, chr=i)



cross <- fill.geno(cross,"no_dbl_XO")

ill.geno(cross, method=c("imp","argmax", "no_dbl_XO", "maxmarginal"), error.prob=0.0001,
Arguments
cross method
error.prob map.function min.prob
map.function=c("haldane","kosambi","c-f","morgan"), min.prob=0.95)


#cross <- calc.errorlod(cross, err=0.05)
#cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
################################################################################

################################################################################
gt <- geno.table(cross)
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 6),])
cross <- pull.markers(cross,bfixA)
################################################################################

################################################################################
png(paste0('~/public_html/ELR_gts_postclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross, chr=i, cex=2)
dev.off()
################################################################################
##cross <- drop.markers(cross,drop)

cross <- tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

pos <- as.numeric(gsub(".*:","",markernames(cross)))
map <- as.numeric(pull.map(cross)[[1]])

if(cor(pos,map, use="complete.obs") < 0){
 cross <<- flip.order(cross, i)
}

cross <- shiftmap(cross, offset=0)

write.cross(cross,chr=i,filestem=filename,format="csv")

################################################################################
### Determine error rate
loglik <- err <- c(0.001, 0.0025, 0.005, 0.01, 0.015, 0.02)
for(z in seq(along=err)) {
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

lod <- (loglik - max(loglik))/log(10)

erpob <- err[which.max(lod)]

cross_map <-  est.map(cross, error.prob=erpob,map.function="kosambi",maxit=100000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

write.cross(cross,chr=i,filestem=filename,format="csv")

print(paste(pop, 'cross written'))
################################################################################

Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))
X <- 1:length(Y)

png(paste0('~/public_html/',pop,'_RF_physpo_concord',i,'_tsp.png'),width=1000,height=500)
par(mfrow=c(1,3))
 plotRF(cross,main=NULL)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
 points(X,Y)
 plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))
dev.off()
