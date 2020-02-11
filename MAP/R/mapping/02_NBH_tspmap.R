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

png(paste0('~/public_00html/',pop,'_gts_preclean',i,'.png'),height=2500,width=4500)
 cross$pheno$gtps <- as.numeric(rowSums(pull.geno(cross) == 2, na.rm = T))
 geno.image(cross, chr=i, reorder=6, cex=2)
dev.off()

################################################################################
map <- pull.map(cross)
newpos <- lapply(map,function(X) { rescale(as.numeric(X),to = c(1,150)) } )
chr <- as.character(i)
cross$geno[[chr]]$map <- as.numeric(unlist(newpos))


### REMOVE SINGLE CROSSOVERS
cross <- removeDoubleXO(cross, chr=i)
cross <- fill.geno(cross, method="no_dbl_XO")

## REMOVE GTs that lead to high errorlof
cross <- calc.errorlod(cross, err=0.05,version="new",map.function="kosambi")
toperr <- top.errorlod(cross, cutoff=5)

print(top.errorlod(cross, cutoff=5))

if ( length(top.errorlod(cross, cutoff=5)[1,]) > 0 ) {
 for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  cross$geno[[chr]]$data[cross$pheno$id==id, mar] <- NA
 }
 cross <- removeDoubleXO(cross, chr=i)
 cross <- fill.geno(cross, method="no_dbl_XO")
}

### REMOVE SHORT DOUBLE CROSSOVERS
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
ab <- cbind(a,b)

ch <- as.character(i)
mat <- cross$geno[[ch]]$data
mat[cbind(a,b)] <- NA

cross$geno[[ch]]$data <- mat

## CLEANUP
cross <- removeDoubleXO(cross, chr=i)
cross <- fill.geno(cross, method="no_dbl_XO")

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
 cross$pheno$gtps <- order(colSums(pull.geno(cross) == 2, na.rm = T))
 geno.image(cross, chr=i, reorder=6, cex=2)
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
