#!/bin/R

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

mapfile <- paste0(pop,'_all_mark_',i,'_tsp')
filename <- file.path(mpath,mapfile)

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

################################################################################

fl <- file.path(paste0(pop,'_unmapped_filtered.csv'))
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

cross  <- subset(cross,chr=i)
sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

################################################################################

ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))
cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)

cross <- subset(cross,ind=!cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))

png(paste0('~/public_html/NBH_gts_preclean',i,'.png'),height=2500,width=4500)
 geno.image(cross, chr=i, cex=2)
dev.off()

################################################################################

cross <- removeDoubleXO(cross, chr=i)
cross <- fill.geno(cross, method="no_dbl_XO")
cross <- calc.errorlod(cross, err=0.05)

xos <- locateXO(cross, full.info=T)
xos <- xos[which(unlist(lapply(xos, is.matrix)))]
indx <- sapply(xos,function(X){
 if(any( X[,'nTypedBetween'] < 4 | is.na(X[,'nTypedBetween']))){
  a <- which(X[,'nTypedBetween'] < 4 | is.na(X[,'nTypedBetween']))
  l <- as.list(X[a,'ileft'])
  r <- as.list(X[a,'iright'])

  a <- mapply(function(Z,Y) { seq(Z,Y) }, Z = l, Y = r )
  as.numeric(unlist(a))
 }
})
ind <- as.numeric(sapply(names(indx), function(x) { which(cross$pheno$ID == x) } ))
a <- rep(ind, times = unlist(lapply(indx,length)))
b <- as.numeric(unlist(indx))
ab <- cbind(a,b)

ch <- as.character(i)
mat <- cross$geno[[ch]]$data
mat[cbind(a,b)] <- NA

cross$geno[[ch]]$data <- mat

cross <- removeDoubleXO(cross, chr=i)
cross <- fill.geno(cross, method="no_dbl_XO")
cross <- calc.errorlod(cross, err=0.05)

#cross <- calc.errorlod(cross, err=0.05)
#cross <- calc.genoprob(cross, step=0, off.end=0, error.prob=0.05, map.function=c("kosambi"),stepwidth=c("fixed"))
#cross <- cleanGeno_jm_2(cross, chr=i, maxdist=25, maxmark=4, verbose=TRUE)
#cross <- fill.geno(cross, min.prob = 0.95 ,method="maxmarginal")
#cross <- fill.geno(cross, error.prob=0.001, method="argmax")
#cross <- removeDoubleXO(cross, chr=i)
#cross <- calc.genoprob(cross, step=0, off.end=0, error.prob=0.01, map.function=c("kosambi"),stepwidth=c("fixed"))
#cross <- fill.geno(cross, min.prob = 0.99 ,method="maxmarginal")
#crossz <- cleanGeno_jm_2(cross, chr=i, maxdist=25, maxmark=4, verbose=TRUE)
#cross <- calc.errorlod(cross, err=0.05)
################################################################################

################################################################################
##gt <- geno.table(cross)
##bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 6),])
##cross <- pull.markers(cross,bfixA)
################################################################################

################################################################################
png(paste0('~/public_html/NBH_gts_postclean',i,'.png'),height=2500,width=4500)
 geno.image(cross, chr=i, cex=2)
dev.off()
################################################################################

cross <- tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

pos <- as.numeric(gsub(".*:","",markernames(cross)))
map <- as.numeric(pull.map(cross)[[1]])

if(cor(pos,map, use="complete.obs") < 0){
 cross <<- flip.order(cross, i)
}

cross <- shiftmap(cross, offset=0)

write.cross(cross,chr=i,filestem=filename,format="csv")

################################################################################

###loglik <- err <- c(0.0001,0.001, 0.0025)
###for(z in seq(along=err)) {
###  cat(z, "of", length(err), "\n")
###  tempmap <- est.map(cross, error.prob=err[z])
###  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
###}
###
###lod <- (loglik - max(loglik))/log(10)
###
###erpob <- err[which.max(lod)]
###
###print(paste('error lod =',erprob))
###
###cross_map <-  est.map(cross, error.prob=erpob,map.function="kosambi",maxit=10000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
###
###cross <- qtl:::replace.map(cross,cross_map)
###
###write.cross(cross,chr=i,filestem=filename,format="csv")
###
###print(paste(pop, 'cross written'))
###################################################################################
###
###Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))
###X <- 1:length(Y)
###
###png(paste0('~/public_html/',pop,'_RF_physpo_concord',i,'_tsp.png'),width=1000,height=500)
###par(mfrow=c(1,3))
### plotRF(cross,main=NULL)
### plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
### points(X,Y)
### plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))
###dev.off()
###
