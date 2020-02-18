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

cross <- removeDoubleXO(cross, chr=chr)

cross <- fill.geno(cross, method="maxmarginal", error.prob = 0.08, min.prob=0.9975)

drop <- names(which(colSums(is.na(pull.geno(cross))) > 5))

cross <- drop.markers(cross,drop)

png(paste0('~/public_html/',pop,'_gts_postclean_mapped',i,'.png'),height=2500,width=4500)
 geno.image(cross, reorder=1, cex=2)
dev.off()

#### MAP #######################################################################

cross <- tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

pos <- as.numeric(gsub(".*:","",markernames(cross)))
map <- as.numeric(pull.map(cross)[[1]])

if(cor(pos,map, use="complete.obs") < 0){
 cross <<- flip.order(cross, i)
}

cross <- shiftmap(cross, offset=0)

write.cross(cross,chr=i,filestem=filename,format="csv")

png(paste0('~/public_html/',pop,'_gts_phenosort_mapped',i,'.png'),height=2500,width=4500)
 geno.image(cross, chr=i, reorder=1, cex=2)
dev.off()

################################################################################
