#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'temp.imp.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
## put chromosomes together
###############################################################################
arg <- paste0(pop,'_imputed_estmap_?[0-9]?[0-9]_tsp.csv')
#arg <- paste0(pop,'_all_mark_imputed_?[0-9]?[0-9]_tsp.csv')
file_list <- list.files(mpath, arg)
#file_list <- c(file_list, 'NBH_imputed_2unmapped_tsp.csv')

cross <- lapply(file_list,function(X){ read.cross(file=X,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)})

gnos <- lapply(cross,function(X){
  data.frame(X[[1]][[1]][['data']],stringsAsFactors=F)
})

ph <- c('Pheno','sex','ID','bin','pheno_norm')

gnos <- do.call(cbind,gnos)
gnos <- cbind(cross[[1]]$pheno[,ph],gnos)
gnos$ID <- as.character(gnos$ID)

m_names <- unlist(sapply(cross,function(X){ markernames(X) }))

colnames(gnos) <- c('Pheno','sex','ID','bin','pheno_norm',m_names)
rownames(gnos) <- cross[[1]]$pheno$ID

ph <- c('Pheno','sex','ID','bin','pheno_norm')
map <- c(colnames(cross[[1]]$pheno[,ph]),unname(unlist(sapply(cross,pull.map))))
chr <- c(colnames(cross[[1]]$pheno[,ph]),gsub(":.*","",m_names))
info <- c(colnames(cross[[1]]$pheno[,ph]),m_names)
headers <- rbind(info,chr,map)
colnames(headers) <- headers[1,]
headers[2:3,1:5] <- ''

headers.u <- unname(data.frame(headers,row.names=NULL,stringsAsFactors=FALSE))
gnos.u <- unname(data.frame(lapply(gnos, as.character),row.names=NULL,stringsAsFactors=FALSE))
colnames(headers.u) <- colnames(gnos.u) <- headers.u[1,]
to_write <- rbind(headers.u,gnos.u)

write.table(to_write, fl, sep=',',row.names=F,quote=F,col.names = F)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

### PLOTS ######################################################################
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross)
plot_test('nbh_mar_regression_hi_confid', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 abline(h=6)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
  points(X,Y)
dev.off()
##############################################################################
if(pop == 'NBH') {
 mfl <- paste0(pop,'_markernames.tsv')
 mfl <- file.path(mpath,mfl)
 write.table(markernames(cross), mfl)
}
################################################################################
