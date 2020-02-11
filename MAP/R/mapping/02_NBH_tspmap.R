#!/bin/R

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

mapfile <- paste0(pop,'_all_mark_',i,'_tsp')
filename <- file.path(mpath,mapfile)

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

library(scales)
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

png(paste0('~/public_html/',pop,'_gts_preclean',i,'.png'),height=2500,width=4500)
 cross$pheno$gtps <- as.numeric(rowSums(pull.geno(cross) == 2, na.rm = T))
 geno.image(cross, chr=i, reorder=6, cex=2)
dev.off()

## SET MAP TO RESONABLE DIST TO CLEAN
chr <- as.character(i)
map <- pull.map(cross)

newpos <- lapply(map,function(X) { setNames(rescale(as.numeric(X),to = c(1,150)),markernames(cross))  } )
attr(newpos,'class') <- 'map'
class(newpos[[chr]]) <- 'A'
attr(newpos[[chr]], "loglik") <- attr(map[[chr]], "loglik")
names(newpos) <- chr
cross <- replace.map(cross,newpos)
print(summary(pull.map(cross)))

### REMOVE SINGLE CROSSOVERS
cross <- removeDoubleXO(cross, chr=chr)
cross <- fill.geno(cross, method="no_dbl_XO", error.prob = 0.08)

### REMOVE PROBLEM MARKERS LEADING TO A CROSSOVER IN > 50% of individuals
xos <- locateXO(cross, full.info=T)
xos <- xos[which(unlist(lapply(xos, is.matrix)))]
indx <- sapply(xos,function(X){
 if(any( X[,'nTypedBetween'] < 3 | is.na(X[,'nTypedBetween']))){
  a <- which(X[,'nTypedBetween'] < 3 | is.na(X[,'nTypedBetween']))
  l <- as.list(X[a,'ileft'])
  r <- as.list(X[a,'iright'])

  a <- mapply(function(Z,Y) { seq(Z,Y) }, Z = l, Y = r )
  as.numeric(unlist(a))
 }
})
ind <- as.numeric(sapply(names(indx), function(x) { which(cross$pheno$ID == x) } ))
a <- rep(ind, times = unlist(lapply(indx,length)))
b <- as.numeric(unlist(indx))

prob.marks <- as.numeric(names(which(table(b) > (nind(cross)/2))))
prob.marks <- markernames(cross)[prob.marks]
cross <- drop.markers(cross,prob.marks)
#################################################################################

### SET VERY SHORT (< 3 markers) CROSSOVERS TO NA
xos <- locateXO(cross, full.info=T)
xos <- xos[which(unlist(lapply(xos, is.matrix)))]
indx <- sapply(xos,function(X){
 if(any( X[,'nTypedBetween'] < 2 | is.na(X[,'nTypedBetween']))){
  a <- which(X[,'nTypedBetween'] < 2 | is.na(X[,'nTypedBetween']))
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

## SET SHORT XO TO NA
ch <- as.character(i)
mat <- cross$geno[[ch]]$data
mat[cbind(a,b)] <- NA
cross$geno[[ch]]$data <- mat

## CLEANUP
cross <- removeDoubleXO(cross, chr=chr)
cross <- fill.geno(cross, method="no_dbl_XO", error.prob = 0.08)
#################################################################################

## REMOVE INDIVIDUAL GTs THAT LEAD TO HIGH ERRORLOD
cross <- calc.errorlod(cross, error.prob = 0.08, version="new", map.function="kosambi")

print('done with errorlod calculation')

if ( length(top.errorlod(cross, cutoff=5)[1,]) > 0 ) {
 toperr <- top.errorlod(cross, cutoff=5)
 for(z in 1:nrow(toperr)) {
  chr <- toperr$chr[z]
  id <- toperr$id[z]
  mar <- toperr$marker[z]
  cross$geno[[chr]]$data[cross$pheno$id==id, mar] <- NA
 }
 cross <- removeDoubleXO(cross, chr=chr)
 cross <- fill.geno(cross, method="no_dbl_XO", , error.prob = 0.08)
}
#################################################################################

################################################################################

png(paste0('~/public_html/',pop,'_gts_postclean',i,'.png'),height=2500,width=4500)
 cross$pheno$gtps <- order(colSums(pull.geno(cross) == 2, na.rm = T))
 geno.image(cross, chr=i, reorder=6, cex=2)
dev.off()

png(paste0('~/public_html/',pop,'_gts_phenosort',i,'.png'),height=2500,width=4500)
 geno.image(cross, chr=i, reorder=1, cex=2)
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
## about 8% error rate
##cross$pheno$gtps <- as.numeric(rowSums(pull.geno(cross) == 2, na.rm = T))/nmar(cross)
##ind <- cross$pheno$ID[which(cross$pheno$gtps > 0.90)]
##
##(0.907113463 + 0.918871252 + 0.941211052)/3
