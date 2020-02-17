pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

source("/home/jmiller1/QTL_agri/MAP/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

mpath <- '/home/jmiller1/QTL_agri/data'
#fl <- paste0(pop,'imp.mapped.tsp.csv')
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

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
ahr_genes$mid_phy <- apply(ahr_genes[,c('str','stp')],1,mean,na.rm=T)
###############

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


get_marks <- function(chr,pos,cross_in = cross){
 phy_vec <- as.numeric(gsub(".*:","",markernames(cross_in,chr)))
 markernames(cross_in,chr)[which.min(abs(phy_vec - pos))]
}

ahr_genes$close_marker <- mapply(get_marks,chr=ahr_genes$chr ,pos=ahr_genes$mid_phy)
ahr_genes$dist <- abs(as.numeric(gsub(".*:","",ahr_genes$close_marker)) - ahr_genes$mid_phy)


########## NOTES ########################
## Segregation distortion at chr 13



##2 27500454 27504907      aip

## skips aip in NBH
crossT <- crossbk
aip <- c("2:27373969","2:27374040","2:27374166","2:27374218","2:27374233","2:27374243","2:27374265","2:27374287","2:27600701","2:27600733","2:27600769","2:27600770","2:27600796","2:27600825","2:27600841","2:27600881","2:27601072","2:27601121","2:27601212","2:27601321")
aip <-  which(markernames(crossT) %in% aip)
aipL <- which(markernames(crossT) == "2:27374287")
aipR <- which(markernames(crossT) == "2:27600701")
round(sm$lod[aip])


##8 16483144 16483822     ARNT
## 18 20388317 20468133    AHR2b
## ALL GENES

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



a <- '/home/jmiller1/QTL_agri/data/out.ldepth'
a <- read.table(a,header=T)

plot_test('depth_2', width=1000, height=10000)
par(mfrow = c(24,1))
for(i in 1:24){
chr <- paste0('chr',i)
dp <- a[which(a$CHROM == chr),]
plot(1:length(dp[,1]),dp[,3], pch=19, col='red', cex=0.5 )
points(1:length(dp[,1]),dp[,4], pch=19, col='blue' , cex=0.5)
}
dev.off()


b <- '/home/jmiller1/QTL_agri/data/filt.out.frq.count'
b <- read.table(b)
b$V5 <- gsub(".*:","",b$V5)
b$V6 <- gsub(".*:","",b$V6)

ords <- order(as.numeric(gsub("chr","",b$V1)))
b <- b[ords,]

plot_test('depth', width=10000, height=1000)
par(mfrow = c(2,1))
plot(1:length(b[,1]),b[,5], pch=19, col=as.factor(b$V1), cex=0.5)
plot(1:length(b[,1]),b[,6], pch=19, col=as.factor(b$V1), cex=0.5)
dev.off()


plot_test('depth', width=1000, height=10000)
par(mfrow = c(24,1))
for(i in 1:24){
chr <- paste0('chr',i)
dp <- b[which(b$V1 == chr),]
plot(1:length(dp[,1]),dp[,5], pch=19, col='red', cex=0.5 )
points(1:length(dp[,1]),dp[,6], pch=19, col='blue' , cex=0.5)
}
dev.off()


b <- '/home/jmiller1/QTL_agri/data/out.ldepth.mean'
b <- read.table(b,header=T)

plot_test('mean_depth', width=10000, height=1000)
plot(1:length(b[,1]),b[,4], pch=19, col=as.factor(b$CHROM), cex=0.5 )

#points(1:length(b[,1]),b[,4], pch=19, col='red' , cex=0.5)
dev.off()
