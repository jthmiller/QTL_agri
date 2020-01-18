#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
## put chromosomes together
###############################################################################

arg <- paste0(pop,'_all_mark_?[0-9]?[0-9]_tsp.csv')
file_list <- list.files(mpath, arg)

arg2 <- paste0(pop,'_all_mark_')
chr <- gsub(arg2,'',file_list)
chr <- as.numeric(gsub("_tsp.csv",'',chr))

cross <- lapply(file_list,function(X){ read.cross(file=X,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)})

gnos <- lapply(cross,function(X){
  data.frame(X[[1]][[1]][['data']],stringsAsFactors=F)
})
gnos <- do.call(cbind,gnos)
gnos <- cbind(cross[[1]]$pheno,gnos)
gnos$ID <- as.character(gnos$ID)

m_names <- unlist(sapply(cross,function(X){
  markernames(X)
}))

colnames(gnos) <- c('Pheno','sex','ID','bin','pheno_norm',m_names)
rownames(gnos) <- cross[[1]]$pheno$ID

map <- c(colnames(cross[[1]]$pheno),unname(unlist(sapply(cross,pull.map))))
chr <- c(colnames(cross[[1]]$pheno),gsub(":.*","",m_names))
info <- c(colnames(cross[[1]]$pheno),m_names)
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
