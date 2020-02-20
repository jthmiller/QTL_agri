
###### ADD UNMAPPED TO TEST FOR LINKAGE
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- 'NBH_imputed_2unmapped_tsp.csv'
fl <- file.path(mpath,fl)
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

 cross <- formLinkageGroups(cross, max.rf = 0.5, min.lod = 0.01, reorgMarkers = TRUE)
 names(cross$geno) <- 'NW'

 mpath <- '/home/jmiller1/QTL_agri/data'
 mapfile <- paste0(pop,'_imputed_2unmapped_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross,filestem=filename,format="csv")
####################################
## COMBINE NBH_imputed_2unmapped_tsp.csv and CHR10
## unmapped NW_012234311.1 are on chr10 and

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137","NBH_6125")
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))

 RF <- 6/nind(cross)
 LOD <- 16
 cross <- formLinkageGroups(cross, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

i <- 'NW'

 cross2 <- cross
 ### REMOVE NON AB AB ###################################################
 dist <- sapply(chrnames(cross2), function(X) { mean(-log10(geno.table(cross2, chr=X)$P.value)) })
 keep <- names(which(dist < 5))
 cross2 <- subset(cross2,chr=keep)
 ####################################################################

## REMOVE MARKERS WITH HIGH RATE OF MISSING DATA
mis <- misg(cross2,0.10)
bfixA <- names(which(colSums(is.na(pull.geno(cross2))) > mis))
print(paste('dropped',length(bfixA),'markers due to missing data'))
cross2 <- drop.markers(cross2, bfixA)
###############################################################################

 ### PVAL filt #################################################################
 gt <- geno.table(cross2)
 pval <- 1.0e-5
 mis <- misg(cross2,0.15)
 bfixA <- rownames(gt[which(gt$P.value > pval & gt$missing < mis),])
 cross2 <- pull.markers(cross2, bfixA)

 swits <- markernames(cross2, chr=chrnames(cross2)[1])
 cross2 <- switchAlleles(cross2, markers = swits)
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)





####################
mpath <- '/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops'

fl <- 'NBH.um.unmapped.f2.csvr'

################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(mpath, fl),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
### Pull names from plinkfile
path <- file.path(mpath, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
################################################################################

 toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137")
 cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))

 cross2 <- cross

 ### PVAL filt #################################################################
 gt <- geno.table(cross2)
 pval <- 1.0e-5
 mis <- misg(cross2,0.15)
 bfixA <- rownames(gt[which(gt$P.value > pval & gt$missing < mis),])
 cross2 <- pull.markers(cross2, bfixA)


 sm <- scanone(cross2, pheno.col=4, model="binary",method="mr")
 chr <- summary(sm)$chr[which(summary(sm)$lod > 3)]
 chr <- unique(c(1:24,as.character(chr)))
 cross2 <- subset(cross2, chr=chr)


### PLOTS ######################################################################
sm <- scanone(cross2, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross2))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross2)
plot_test('nbh_mar_regression_hi_confid', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 abline(h=6)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
  points(X,Y)
dev.off()

######################################################################

unmp <- grep('NW',chrnames(cross2), value=T)
cross.write <- subset(cross2,chr = unmp)


cross.write <- formLinkageGroups(cross.write, max.rf = 0.5, min.lod = 0.1, reorgMarkers = TRUE)
names(cross.write$geno) <- 'NW'


 mpath <- '/home/jmiller1/QTL_agri/data'
 mapfile <- paste0(pop,'_order_impute_NW_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross.write,filestem=filename,format="csv")


cross.write2 <- formLinkageGroups(cross.write, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)



while read -r file
do
cp $file NW_${file}

done <<< "NBH_order_impute_*[0-9]_tsp.csv"






#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

################################################################################
## put chromosomes together
###############################################################################
arg <- paste0(pop,'_imputed_estmap_?[0-9]?[0-9]_tsp.csv')
#arg <- paste0(pop,'_all_mark_imputed_?[0-9]?[0-9]_tsp.csv')
file_list <- list.files(mpath, arg)
file_list <- c(file_list[2], 'NBH_imputed_2unmapped_tsp.csv')

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

fl <- paste0(pop,'_10_unmapped.tsp.csv')
fl <- file.path(mpath,fl)

write.table(to_write, fl, sep=',',row.names=F,quote=F,col.names = F)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

 RF <- 6/nind(cross)
 LOD <- 16
 cross2 <- formLinkageGroups(cross, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)



plotit(cross2)
