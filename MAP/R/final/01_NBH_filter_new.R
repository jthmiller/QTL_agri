#!/bin/R
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

################################################################################
## read in the QTL cross
umpath <- '/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops'
fl <- 'NBH.um.unmapped.f2.csvr'
cross <- read.cross(file = file.path(umpath, fl),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

#################################################################################
### read in the QTL cross
#cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
#format = "csvr", geno = c(1:3), estimate.map = FALSE)
#################################################################################

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

#### SEX #######################################################################
sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec
################################################################################

#### SWITCH ALLELES THAT ARE PROB AA x BB ######################################
m <- which(cross$pheno$ID=='NBH_NBH1M')
f <- which(cross$pheno$ID=='NBH_NBH1F')
bfixm <- pull.geno(cross)[m,]
bfix_swit1 <- names(bfixm)[which(as.numeric(bfixm)==1)]
bfixf <- pull.geno(cross)[f,]
bfix_swit2 <- names(bfixf)[which(as.numeric(bfixf)==3)]
bfix_swit12 <- unique(c(bfix_swit1 ,bfix_swit2))
cross <- switchAlleles(cross, markers = bfix_swit12)
################################################################################

################################################################################
## Parent markers
## AA, AB, BB are displayed in the colors red, blue, and green,
################################################################################
parc1 <- subset(cross,ind=c('NBH_NBH1M','NBH_NBH1F'))
plot_test('par_nbh_unfilt', width=4000, height = 5000)
par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc1, chr=i)} ; dev.off()
################################################################################

################################################################################
## drop invariant and ABxAB cross in grand parents
m <- which(cross$pheno$ID=='NBH_NBH1M')
f <- which(cross$pheno$ID=='NBH_NBH1F')
pars <- which(cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))

bfixbk <- pull.geno(cross)
drop <- names(which(bfixbk[m,] == bfixbk[f,]))
drop.na <- names(which(is.na(bfixbk[m,]) & is.na(bfixbk[f,])))
table(bfixbk[pars,drop])
cross <- drop.markers(cross,c(drop.na,drop))
################################################################################

################################################################################
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.50))
### is "NBH_5646" another grandparent sample??
### after filtering, NBH_6137 appears to have high allelic dropout
toss.missing <- c(toss.missing,"NBH_5646","NBH_6137")
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.125)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
print(paste('dropping',length(drop),'markers'))
cross <- drop.markers(cross,drop)
################################################################################

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers on all LGs)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-3),])
cross <- drop.markers(cross,toss)
################################################################################

################################################################################
### Sliding window filter to target non-distorted regions ######################
mean.dist <- function(X, phys, gtb){
  ind <- which(phys > X[1] & phys < X[2])
  mean(gtb[ind,'P.value'])
  }

seg.window <- lapply(1:24,function(Z) {
 gta <- geno.table(cross, chr = Z)
 phys <- as.numeric(gsub(".*:",'',rownames(gta)))
 X1 <- seq(0, max(phys), by = 100000)
 Y1 <- c(seq(0, max(phys), by = 100000)[-1], max(phys))

 pos <- apply(cbind(X1,Y1), 1, mean)
 #print(head(gta))
 #pval <- apply(mean.dist, X = X1, Y = Y1, phys = phys, gtb = gta)
 pval <- apply(cbind(X1,Y1),1, mean.dist,phys = phys, gtb = gta)
cbind(pos,pval)
})

### CHR2 < 13150000
seg.window[[2]]
chr2 <- which(as.numeric(gsub(".*:",'', markernames(cross, 2))) > 13150000)
chr2 <- markernames(cross, 2)[chr2]

### CHR13 < 12350000
seg.window[[13]]
chr13 <- which(as.numeric(gsub(".*:",'', markernames(cross, 13))) > 12350000)
chr13 <- markernames(cross, 13)[chr13]

################################################################################
### ALL BUT 2, 13 can be filtered down to 9e-3. Truncates these LGS (see plots)
gt.sub <- geno.table(cross)
gt.sub <- gt.sub[!rownames(gt.sub) %in% c(chr2,chr13),]
toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
cross <- drop.markers(cross,toss.sub)
################################################################################

#################################################################################
### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_switched_filtered')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
#################################################################################

crossbk <- cross

### ONE MARKER PER RAD SITE ####################################################
cross <- thin_by_distortion(cross,5)
################################################################################

################################################################################
### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_switched_filtered_thinned_NW')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
################################################################################

################################################################################
parc2 <- pull.markers(parc1,markernames(cross))
plot_test('par_nbh_unfilt', width=4000, height = 5000)
par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc2, chr=i)} ; dev.off()
################################################################################

#################################################################################
#################################################################################
### ASSIGNING UNMAPPED SCAFFOLDS ################################################
### TEST WHICH UNMAPPED SCAFFOLDS ARE LINKED ####################################
#nw_chr <- grep('NW_',chrnames(cross), value = T)
#cross <- removeDoubleXO(cross,nw_chr)
#
#chr <- 1:24
#cross <- est.rf(cross)
#
###########################
#ldm <- function(nw) {
# sapply(chr,function(z){ mean(pull.rf(cross, what='lod')[markernames(cross,nw),markernames(cross,z)]) })
#}
###########################
#
#ld <- foreach(nw = nw_chr, .inorder = F, .packages = libs2load) %dopar% ldm(nw)
#ld <- do.call(rbind,ld)
#rownames(ld) <- nw_chr
#
#nms <- which(apply(ld,1,max,na.rm=T) > 5)
#reassign <- apply(ld,1,which.max)
#reassign <- reassign[nms]
#
#nw_marks_assign <- sapply(names(reassign),markernames,cross = cross)
#nw_length <- sapply(nw_marks_assign,length)
#nw_marks_assign <- as.character(unlist(nw_marks_assign))
#nw_ch <- rep(as.numeric(reassign), times = as.numeric(nw_length))
#nw_pos <- unlist(sapply(nw_length,seq,from = 1, by = 1))
#nw_old <- gsub(":.*","",nw_marks_assign)
#
#### WRTIE THE TABLE TO REASSIGN THE UNMAPPED MARKERS ###########################
#move <- data.frame(cbind(nw_old,nw_marks_assign,nw_ch,nw_pos), stringsAsFactors=F)
#movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
#write.table(move,movefl)
#################################################################################
#################################################################################

################################################################################
#### READ IN THE CROSS
#fl <- paste0(pop,'_filtered_pvalue_thinned_NW')
###fl <- paste0(pop,'_filtered_unphased_NW.csv')
#cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#################################################################################

################################################################################
### READ THE UNMAPPED MARKER ASSIGNMENT TABLE
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
move <- move[which(move$nw_marks_assign %in% markernames(cross)),]
################################################################################

### ASSIGN UNMAPPED MARKERS ####################################################
for (i in 1:length(move[,1])){
 cross <<- movemarker(cross, marker = move[i,'nw_marks_assign'], newchr = move[i,'nw_ch'], newpos = as.numeric(move[i,'nw_pos']))
 print(i)
}
cross <- subset(cross,chr=1:24)
################################################################################

### LATER SHOWN TO BE BAD MARKERS
###################################################################################
##drop <- c('15:21481705','20:16313414','22:8254649','24:3192380','NW_012234461.1:562796','NW_012234461.1:1059423','NW_012234494.1:872524','8:12257474','NW_012234558.1:739273','NW_012234311.1:628953','NW_012234311.1:166994','NW_012234311.1:628953','NW_012234311.1:166994','NW_012234520.1:169476','NW_012225526.1:51452','NW_012234326.1:397833','NW_012234326.1:380172')
##cross <- drop.markers(cross,drop)
##
##drop2 <- c('10:9960308','NW_012234311.1:3468744', 'NW_012234311.1:3445534','NW_012234311.1:3383007','NW_012234326.1:2162519','NW_012224981.1:180658','NW_012224824.1:3425')
##cross <- drop.markers(cross,drop2)
###################################################################################

#################################################################################
### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_filtered_pvalue_thinned_NW_moved')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
#################################################################################

################################################################################
### READ IN THE CROSS
fl <- paste0(pop,'_filtered_pvalue_thinned_NW_moved.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
linked_marks <- function(X, LOD = 12, RF = 1){
 crossX <- est.rf(subset(cross,chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}
linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(X)
cross <- pull.markers(cross,unlist(linked))
################################################################################

################################################################################
switched_marks <- function(X){
 checkAlleles(subset(cross,chr=X), threshold = 2)
}
switched <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% switched_marks(X)
switched <- switched[!sapply(switched,is.null)]
switched <- do.call(rbind,switched)
cross <- switchAlleles(cross, as.character(switched$marker))
################################################################################

################################################################################
drop.unkinked <- function(i, cross){
 cross.sub <- subset(cross,chr=i)
 nw_chr <- grep('NW_',chrnames(cross.sub), value = T)
 cross.sub <- drop.markers(cross.sub, nw_chr)
 cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
 cross.sub <- tspOrder(cross = cross.sub, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 cross.sub <- removeDoubleXO(cross.sub)
 ### HIGH MISSING DATA DUE TO REMOVING XOs ######################################
 mis <- misg(cross.sub,0.125)
 drop <- names(which(colSums(is.na(pull.geno(cross.sub))) > mis))
 cross.sub <- drop.markers(cross.sub,drop)
 return(c(markernames(cross.sub),nw_chr))
}
################################################################################

################################################################################
keep <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% drop.unkinked(i, cross = cross)
cross <- pull.markers(cross,unlist(keep))
################################################################################

################################################################################
cross <- est.rf(cross)
cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross <- removeDoubleXO(cross)
cross <- est.rf(cross)
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_estrf_no_impute.rsave')))
################################################################################

################################################################################
### WRITE IMPUTED AND CORRECTED MAP
mapfile <- paste0(pop,'_',sum(nmar(cross5)),'_estrf_no_impute')
filename <- file.path(mpath,mapfile)
write.cross(cross5,filestem=filename,format="csv")
################################################################################

################################################################################
### FIX GENOTYPING ERRORS WITH GENOPROB
cross <- fill.geno(cross, method="maxmarginal", error.prob = 0.05, min.prob=0.95)
mis <- misg(cross,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

################################################################################
cross <- est.rf(cross)
cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross <- fill.geno(cross, method="no_dbl_XO", error.prob = 0.01, min.prob=0.99)
################################################################################

################################################################################
### WRITE IMPUTED AND CORRECTED MAP
mapfile <- paste0(pop,'_',sum(nmar(cross)),'_estrf_imputed')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_estrf_no_impute.rsave')))
################################################################################
