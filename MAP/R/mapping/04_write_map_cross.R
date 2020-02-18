#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'temp.imp.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
## put chromosomes together
###############################################################################
arg <- paste0(pop,'_imputed_?[0-9]?[0-9]_tsp.csv')
#arg <- paste0(pop,'_all_mark_imputed_?[0-9]?[0-9]_tsp.csv')
file_list <- list.files(mpath, arg)
file_list <- c(file_list, 'NBH_imputed_2unmapped_tsp.csv')

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

##############################################################################
if(pop == 'NBH') {
 mfl <- paste0(pop,'_markernames.tsv')
 mfl <- file.path(mpath,mfl)
 write.table(markernames(cross), mfl)
}
################################################################################
