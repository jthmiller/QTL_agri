###############
AHR.bed <- read.table(file.path(mpath,"lift_AHR_genes.bed"), stringsAsFactors = F, header = F)
colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
AHR.bed$str <- as.numeric(AHR.bed$str)
AHR.bed$stp <- as.numeric(AHR.bed$stp)
AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]
AHR.bed$gene <- gsub(":158640", "", AHR.bed$gene)
AHR.bed <- AHR.bed[!AHR.bed$chr == 5,]
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
ahr_genes <- cnv.ahrs(cross, AHRdf = AHR.bed, EXP = F)
###############
##2 27500454 27504907      aip

## skips aip
crossT <- crossbk
aip <- c("2:27373969","2:27374040","2:27374166","2:27374218","2:27374233","2:27374243","2:27374265","2:27374287","2:27600701","2:27600733","2:27600769","2:27600770","2:27600796","2:27600825","2:27600841","2:27600881","2:27601072","2:27601121","2:27601212","2:27601321")
aip <-  which(markernames(crossT) %in% aip)
aipL <- which(markernames(crossT) == "2:27374287")
aipR <- which(markernames(crossT) == "2:27600701")
round(sm$lod[aip])


##8 16483144 16483822     ARNT
## 18 20388317 20468133    AHR2b
## ALL GENES
genes.bed <- read.table(file.path(mpath,"lifted_genes.bed"), stringsAsFactors = F, header = T)
genes.bed$chr <- gsub('chr','',genes.bed$chr)
genes.bed <- genes.bed[genes.bed$chr %in% c(1:24),]
genes.bed$mid <- round(apply(genes.bed[c('start','end')],1,mean))

get_genes <- function(chr,pos){
 whole_chrom <- genes.bed[which(genes.bed$chr == chr),]
 min.ind <- which.min(abs(whole_chrom$mid - pos))
 near <- c(min.ind -1,min.ind,min.ind+1)
 dis <- whole_chrom[near,'mid'] - pos
 cbind(whole_chrom[near,],dis)
}

aip_peak2 <- get_genes(37728488,2)
## bad assembly
bad_marker <- get_genes(14,12556225)
bad_marker <- get_genes(19,2435634)

## 18 20388317 20468133    AHR2b
max18 <- '18:17538823'
get_genes(18,17538823)

## close to above
18:20367708   1 1060 5.2823253
18:20367780   1 1070 5.2823253
