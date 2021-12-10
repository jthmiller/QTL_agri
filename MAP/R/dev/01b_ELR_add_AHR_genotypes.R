#!/bin/R
### Map QTLs 1 of 3
pop <- 'ELR'
debug.cross <- T
#source("/home/jmiller1/QTL_agri/MAP/control_file.R")
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

library('qtl')
mpath <- '/home/jmiller1/QTL_agri/data'

################################################################################
## ADD AHR GENOTYPES ##
################################################################################

fl <- file.path(mpath,'ELR_unmapped_filtered.csv')
cross.df <- read.csv(fl,header=FALSE,stringsAsFactors=F)

marks_nms <- cross.df[1,6:length(cross.df[1,])]
gts <- cross.df[4:length(cross.df[,1]),6:length(cross.df[1,])]
## 88 X 19856
rownames(gts) <- cross.df[c(4:length(cross.df[,1])),'V3']
colnames(gts) <- as.character(cross.df[1,6:length(cross.df[1,])])
## gts 88x19856
nmars <- length(gts[1,])
phenotpyes <- cross.df[,1:5]
rownames(phenotpyes ) <- c('info','chr','map',c(cross.df[4:length(cross.df[,1]),'V3']))
### phen 89x5

################################################################################
fla <-file.path(mpath, 'ER_ahr_aip_whoi_gt.csv')
cross.df.ahr <- read.csv(fla,header=FALSE,stringsAsFactors=F)
cross.df.ahr[which(cross.df.ahr$V3 == 'ELR_10869'),1:3] <- cross.df[which(cross.df$V3 == "ELR_11115"),1:3]
cross.df.ahr[which(cross.df.ahr$V3 == "ELR_11115"),4:6] <- cbind("-","-","-")

add_ahr <- cbind(cross.df[which(cross.df$V3 == "ELR_11115"),1:3],"-","-","-")
rownames(add_ahr) <- "ELR_11115"
rbind(cross.df.ahr,add_ahr)
ahr_mark_nms <- cross.df.ahr[1,4:length(cross.df.ahr[1,])]

## AHR genotypes
ahr_gts <- cross.df.ahr[4:length(cross.df.ahr[,1]),4:length(cross.df.ahr[1,])]
rownames(ahr_gts) <- cross.df.ahr[c(4:length(cross.df.ahr[,1])),'V3']
colnames(ahr_gts) <- cross.df.ahr[1,c(4:6)]

##ahr_gts <- ahr_gts[which(!rownames(ahr_gts) == "ELR_10869" ),]
### ahr.gts 87x3

##pheno data with blank rows included
phen.ah <- cross.df.ahr[,1:3]
rownames(phen.ah) <- c('info','chr','map',phen.ah[4:length(phen.ah[,1]),'V3'])
## phen.ah 90x3
################################################################################

### add ind to ahr data to make df even number of ind
all_ind <- rownames(gts)[which(!rownames(gts) %in% rownames(ahr_gts))]
b <- matrix('-',ncol=ncol(ahr_gts),nrow=length(all_ind))
rownames(b) <- all_ind
colnames(b) <- colnames(ahr_gts)
ahr_gts <- rbind(ahr_gts, b)
ahr_ind <- rownames(ahr_gts)[which(!rownames(ahr_gts) %in% rownames(gts))]
##new <- rownames(ahr_gts)[!which(rownames(ahr_gts) %in% rownames(gts.2))]

a <- matrix('-',ncol=nmars,nrow=length(ahr_ind))
rownames(a) <- ahr_ind
colnames(a) <- colnames(gts)
gts.2 <- rbind(gts, a)
## add ahr genotypes
gts.2 <- cbind(gts.2, ahr_gts[rownames(gts.2),])
gts.2 <- cbind(phen.ah[rownames(gts.2),],NA,NA,gts.2)
colnames(gts.2)[1:5] <- c('Pheno','sex','ID','bin','pheno_norm')

final.gts <- gts.2

row1 <- colnames(gts.2)
row2 <- c(cross.df[2,],c(1,2,2))
names(row2) <- colnames(gts.2)
row3 <- c(cross.df[3,],c(0,0,0))
names(row3) <- colnames(gts.2)

final.gts <- rbind(row1,row2,row3,final.gts)

fl <- file.path(mpath,'ELR_unmapped_added_markers.csv')
write.table(final.gts, fl,col.names=F,row.names=F,quote=F,sep=',')

################################################################################
################################################################################

fl <- file.path(mpath,'ELR_unmapped_added_markers.csv')

cross2 <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

all_marks <- c(markernames(cross),"AHR2a_del","AIP_252","AIP_261")

cross2_unmapped <- pull.markers(cross2, all_marks)

cross2_unmapped <- subset(cross2_unmapped, ind = cross2_unmapped$pheno$ID %in% cross$pheno$ID)


## map with deletion
noperms  <- est.rf(cross2_unmapped)
noperms  <- tspOrder(cross = noperms , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
newmap <- est.map(noperms, error.prob = erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
noperms <- replace.map(noperms, newmap)

i <- 100 ; plotit(noperms, 'cross_noimp_exact_nopar')


## is AHRdel linked to CHR1? What about chr18?

## perform marker regression --> pvalue 

erprob <- 0.001 #determined from mapping step
noperms <- sim.geno(noperms, n.draws = 150, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
noperms <- calc.genoprob(noperms, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
################################################################################
noperms$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)

sone_bin <- scanone(noperms, pheno.col=4, method="hk", model="bin")
sone_reg <- scanone(noperms, pheno.col=4, method="mr", model="bin")
########
