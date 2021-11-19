#!/bin/R

### Split dataset: A: mapping cross B: Genotype validation for exploring causative variation

pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

## Map a smaller subset to id qtls
## go back and remap chrs with a qtl with a larger sample set

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
### Pull sample names from plinkfile
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

### RAW DATA SET ###############################################################
mapfile <- paste0(pop,'_raw')
filename <- file.path(mpath,mapfile)
write.cross(cross, filestem=filename, format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################
i <- 1 ; plotit(cross,'_raw')

################################################################################
cg <- comparegeno(cross)
plot_test('raw_rela_nbh')
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
dev.off()
################################################################################

################################################################################
### FILTER
################################################################################
## Drop samples that have high missing data or discordant genotype
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.25))
### is "NBH_5646" another grandparent sample??
### after filtering, NBH_6137 appears to have high allelic dropout
toss.missing <- c(toss.missing,"NBH_5646")
cross <- subset(cross, ind=!cross$pheno$ID %in% toss.missing)
################################################################################

################################################################################
## Parent markers
## AA, AB, BB are displayed in the colors red, blue, and green,
################################################################################
parc1 <- subset(cross,ind=c('NBH_NBH1M','NBH_NBH1F'))
plot_test('par_nbh_prefilt', width=4000, height = 5000)
par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc1, chr=i)} ; dev.off()
################################################################################

################################################################################
## drop invariant and ABxAB cross in grand parents
m <- which(cross$pheno$ID=='NBH_NBH1M')
f <- which(cross$pheno$ID=='NBH_NBH1F')
pars <- which(cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
bfixbk <- pull.geno(cross)
drop.many <- names(which(bfixbk[m,] == bfixbk[f,]))
drop.fewer <- names(which(bfixbk[m,] == bfixbk[f,] & bfixbk[m,] %in% c(1,3) ))
drop.na <- names(which(is.na(bfixbk[m,]) & is.na(bfixbk[f,])))
table(bfixbk[pars,drop.many])
################################################################################

## KEEP ABxAB markers ##########################################################
cross.many <- drop.markers(cross,drop.fewer)
mapfile <- paste0(pop,'_keep_het')
filename <- file.path(mpath,mapfile)
write.cross(cross.many, filestem=filename, format="csv")
i <- 1 ; plotit(cross.many,'_lib')
################################################################################


## CONSERVATIVE ################################################################
cross.cons <- drop.markers(cross,c(drop.na,drop.many))
mapfile <- paste0(pop,'_cons_filt')
filename <- file.path(mpath,mapfile)
write.cross(cross.cons, filestem=filename, format="csv")
i <- 1 ; plotit(cross.cons,'_cons')
################################################################################

#mapfile <- paste0(pop,'_keep_het')
#cross.many <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

#mapfile <- paste0(pop,'_cons_filt')
#cross.cons <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

################################################################################
## Parent markers
## AA, AB, BB are displayed in the colors red, blue, and green,
################################################################################
parc2 <- subset(cross.many, ind=c('NBH_NBH1M','NBH_NBH1F'))
plot_test('par_nbh_keep_many', width=4000, height = 5000)
par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc2, chr=i)} ; dev.off()
################################################################################

## DROP PARENTS ################################################################
cross.many <- subset(cross.many, ind=!cross.many$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
cross.cons <- subset(cross.cons, ind=!cross.cons$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
################################################################################

sum(nmar(cross.many))
sum(nmar(cross.cons))

#### REFINE GENOTYPES ##########################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X, perc) { nind(cross.cons) * perc }
## Drop markers with greater than 12.5% missing data
mis <- misg(cross.cons,0.125)
drop <- names(which(colSums(is.na(pull.geno(cross.cons))) > mis))
print(paste('dropping',length(drop),'markers'))
cross.cons <- drop.markers(cross.cons,drop)
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X, perc) { nind(cross.many) * perc }
## Drop markers with greater than 12.5% missing data
mis <- misg(cross.many,0.125)
drop <- names(which(colSums(is.na(pull.geno(cross.many))) > mis))
print(paste('dropping',length(drop),'markers'))
cross.many <- drop.markers(cross.many,drop)
################################################################################

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers on all LGs)
gt <- geno.table(cross.cons)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-3),])
cross.cons <- drop.markers(cross.cons,toss)
i <- 5 ; plotit(cross.cons,'cons')
################################################################################

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers on all LGs)
gt <- geno.table(cross.many)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-3),])
cross.many <- drop.markers(cross.many,toss)
i <- 6 ; plotit(cross.many,'many')
################################################################################

################################################################################
### ALL BUT 2, 7, 13, 24 can be filtered down to 9e-3. Truncates these LGS (see plots)
gt.sub <- geno.table(cross.cons, c(1,3:6,8:12,14:23))
##gt.sub <- gt.sub[!gt.sub$chr %in% c(2,7,13) | is.na(as.numeric(as.character(gt.sub$chr))),]
toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
cross.cons <- drop.markers(cross.cons,toss.sub)
i <- 7 ; plotit(cross.cons,'cross.cons')
################################################################################

################################################################################
### ALL BUT 2, 7, 13 can be filtered down to 9e-3. Truncates these LGS (see plots)
gt.sub <- geno.table(cross.many, c(1,3:6,8:12,14:23))
##gt.sub <- gt.sub[!gt.sub$chr %in% c(2,7,13) | is.na(as.numeric(gt.sub$chr)),]
toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
cross.many <- drop.markers(cross.many,toss.sub)
i <- 8 ; plotit(cross.many,'cross.many')
################################################################################

################################################################################
mapfile <- paste0(pop,'_many_filt_refined')
filename <- file.path(mpath,mapfile)
write.cross(cross.many, filestem=filename, format="csv")

mapfile <- paste0(pop,'_cons_filt_refined')
filename <- file.path(mpath,mapfile)
write.cross(cross.cons, filestem=filename, format="csv")
################################################################################

sum(nmar(cross.many))
sum(nmar(cross.cons))

#mapfile <- paste0(pop,'_many_filt_refined')
#cross.many <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

#mapfile <- paste0(pop,'_cons_filt_refined')
#cross.cons <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
nw_marks <- grep('^NW*', unique(c(markernames(cross.many), markernames(cross.cons))), value = T)

linked_marks <- function(cross, X, LOD = 8, RF = 1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross.many, X)
cross.many <- pull.markers(cross.many, c(unlist(linked),nw_marks))
i <- 95 ; plotit(cross.many,'cross_linked')

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross.cons, X)
cross.cons <- pull.markers(cross.cons, c(unlist(linked),nw_marks))
i <- 95 ; plotit(cross.cons,'cross_linked')
################################################################################

################################################################################
switch.phase <- function(i, cross){
 cross.sub <- subset(cross,chr=i)
 cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
 rf <- pull.rf(cross.sub)
 freq <- rowSums(rf, na.rm=T)/dim(rf)[1]
 return(names(freq)[which(freq > 0.65)])
}
################################################################################

################################################################################
switch.many <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% switch.phase(i, cross = cross.many)
cross.many <- switchAlleles(cross.many, as.character(unlist(switch.many)))
i <- 95 ; plotit(cross.many,'cross_linke_sw')

switch.cons <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% switch.phase(i, cross = cross.cons)
cross.cons <- switchAlleles(cross.cons, as.character(unlist(switch.cons)))
i <- 95 ; plotit(cross.cons,'cross_linked_sw')
################################################################################

################################################################################
mapfile <- paste0(pop,'_many_filt_refined_switch')
filename <- file.path(mpath,mapfile)
write.cross(cross.many, filestem=filename, format="csv")

mapfile <- paste0(pop,'_cons_filt_refined_switch')
filename <- file.path(mpath,mapfile)
write.cross(cross.cons, filestem=filename, format="csv")
################################################################################

################################################################################

## READ IN
mapfile <- paste0(pop,'_many_filt_refined_switch')
cross.many <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
sum(nmar(cross.many))


mapfile <- paste0(pop,'_cons_filt_refined_switch')
cross.cons <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
sum(nmar(cross.cons))
################################################################################

##cross.many_rm <- findDupMarkers(cross.many, exact.only=FALSE, adjacent.only=TRUE) # finds 6 pairs
##cross.many <- drop.markers(cross.many, unlist(cross.many_rm))
##cross.cons_rm <- findDupMarkers(cross.cons, exact.only=FALSE, adjacent.only=TRUE) # finds 6 pairs
##cross.cons <- drop.markers(cross.cons, unlist(cross.cons_rm))
##i <- 9 ; plotit(cross.many,'cross.many')
##i <- 9 ; plotit(cross.cons,'cross.cons')
##cross_cons_noDXO <- removeDoubleXO(cross.cons)
##cross_many_noDXO <- removeDoubleXO(cross.many)

#cross.many_imp <- removeDoubleXO(cross.many_imp)
#cross.cons_imp  <- removeDoubleXO(cross.cons_imp)
#i <- 9 ; plotit(cross.many_imp,'cross.cons')
#i <- 9 ; plotit(cross.many_imp,'cross.many')

################################################################################
## IMPUTATION
################################################################################
################################################################################
################################################################################

cross.many_imp <- fill.geno(cross.many, error.prob = 0.0001, map.function= 'kosambi')
i <- 9 ; plotit(cross.many_imp,'cross.many')

cross.cons_imp <- fill.geno(cross.cons, error.prob = 0.0001, map.function= 'kosambi')
i <- 9 ; plotit(cross.cons_imp,'cross.cons')

cross.many_rm <- findDupMarkers(cross.many_imp, exact.only=FALSE, adjacent.only=TRUE) # finds 6 pairs
cross.many_imp <- drop.markers(cross.many_imp, unlist(cross.many_rm))
i <- 10 ; plotit(cross.many_imp,'cross.many')

cross.cons_rm <- findDupMarkers(cross.cons_imp, exact.only=FALSE, adjacent.only=TRUE) # finds 6 pairs
cross.cons_imp <- drop.markers(cross.cons_imp, unlist(cross.cons_rm))
i <- 10 ; plotit(cross.cons_imp,'cross.cons')

################################################################################
### READ THE UNMAPPED MARKER ASSIGNMENT TABLE
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
move.cons <- move[which(move$nw_marks_assign %in% markernames(cross.cons_imp)),]
move.many <- move[which(move$nw_marks_assign %in% markernames(cross.many_imp)),]
################################################################################

### ASSIGN UNMAPPED MARKERS ####################################################
for (i in 1:length(move.many[,1])){
 cross.many_imp <<- movemarker(cross.many_imp, marker = move.many[i,'nw_marks_assign'], newchr = move.many[i,'nw_ch'], newpos = as.numeric(move.many[i,'nw_pos']))
 print(i)
}
for (i in 1:length(move.cons[,1])){
 cross.cons_imp <<- movemarker(cross.cons_imp, marker = move.cons[i,'nw_marks_assign'], newchr = move.cons[i,'nw_ch'], newpos = as.numeric(move.cons[i,'nw_pos']))
 print(i)
}
cross.many_imp <- subset(cross.many_imp,chr=1:24)
cross.cons_imp <- subset(cross.cons_imp,chr=1:24)
################################################################################

################################################################################
linked_marks <- function(cross, X, LOD = 10, RF = 1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross.many_imp, X)
cross.many_imp <- pull.markers(cross.many_imp, unlist(linked))
i <- 95 ; plotit(cross.many_imp,'many_cross_linked')

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross.cons_imp, X)
cross.cons_imp <- pull.markers(cross.cons_imp, unlist(linked))
i <- 95 ; plotit(cross.cons_imp,'cons_cross_linked')
################################################################################

################################################################################
switch.phase <- function(i, cross){
 cross.sub <- subset(cross,chr=i)
 cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
 rf <- pull.rf(cross.sub)
 freq <- rowSums(rf, na.rm=T)/dim(rf)[1]
 return(names(freq)[which(freq > 0.65)])
}
################################################################################

################################################################################
switch.cons <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% switch.phase(i, cross = cross.cons_imp)
cross.cons_imp <- switchAlleles(cross.cons_imp, as.character(unlist(switch.cons)))
i <- 95 ; plotit(cross.cons_imp,'cross_cons_sw')

switch.many <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% switch.phase(i, cross = cross.many_imp)
cross.many_imp <- switchAlleles(cross.many_imp, as.character(unlist(switch.many)))
i <- 95 ; plotit(cross.many_imp,'cross_many_sw')
################################################################################

################################################################################
## MAP
################################################################################
cross.cons_imp  <- est.rf(cross.cons_imp)
cross.cons_imp  <- tspOrder(cross = cross.cons_imp , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(cross.cons_imp,'cross.cons_imp')

cross.many_imp  <- est.rf(cross.many_imp)
cross.many_imp  <- tspOrder(cross = cross.many_imp , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(cross.many_imp,'cross.many_imp')

save.image(file.path(mpath,paste0(pop,'_dups_dropped_tsp.rsave')))
################################################################################
################################################################################
################################################################################
################################################################################












################################################################################

################################################################################
cross.cons_imp  <- est.rf(cross.cons_imp)
cross.cons_imp  <- tspOrder(cross = cross.cons_imp , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(cross.cons_imp,'cross.cons_imp')

cross.many_imp  <- est.rf(cross.many_imp)
cross.many_imp  <- tspOrder(cross = cross.many_imp , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(cross.many_imp,'cross.many_imp')

save.image(file.path(mpath,paste0(pop,'_dups_dropped_tsp.rsave')))


cross.cons_imp .map <- est.map(cross.cons_imp , error.prob=0.001, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross.cons_imp  <- replace.map(cross.cons_imp , cross.cons_imp .map)

cross_noDXO <- removeDoubleXO(cross)
i <- 1 ; plotit(cross_noDXO,'dups_dropped_tsp')


################################################################################


i <- 95 ; plotit(cross_linked_switched,'cross_linked_switched')


cross.many_noDXO <- removeDoubleXO(cross.many)
cross.cons_noDXO <- removeDoubleXO(cross.cons)
i <- 9 ; plotit(cross.cons_noDXO,'cross.cons')
i <- 9 ; plotit(cross.many_noDXO,'cross.many')

#save.image(file.path(mpath,paste0(pop,'_noOX.rsave')))
load(file.path(mpath,paste0(pop,'_noOX.rsave')))







cg <- comparegeno(cross)
plot_test('raw_rela_new')
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
dev.off()








cross.many_0.01 <- fill.geno(cross.many, error.prob = 0.01, map.function= 'kosambi')
i <- 9 ; plotit(cross.many_0.01,'cross.many')

cross.cons_0.01 <- fill.geno(cross.cons, error.prob = 0.01, map.function= 'kosambi')
i <- 9 ; plotit(cross.cons_0.01,'cross.cons')








################################################################################
cross.cons_0.05 <- fill.geno(cross.cons, error.prob = 0.05, map.function= 'kosambi')
i <- 9 ; plotit(cross.cons_0.05,'filled_cons_0.05a')
gt.sub <- geno.table(cross.cons_0.05)
gt.sub <- gt.sub[!gt.sub$chr %in% c(2,7,13,24),]
toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
cross.cons_0.05 <- drop.markers(cross.cons_0.05,toss.sub)
i <- 10 ; plotit(cross.cons_0.05,'filled_cons_0.05b')

cross.many_0.05 <- fill.geno(cross.many, error.prob = 0.05, map.function= 'kosambi')
i <- 11 ; plotit(cross.many_0.05,'filled_many_0.05a')
gt.sub <- geno.table(cross.many_0.05)
gt.sub <- gt.sub[!gt.sub$chr %in% c(2,6,13,24),]
toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
cross.many_0.05 <- drop.markers(cross.many_0.05,toss.sub)
i <- 12 ; plotit(cross.many_0.05,'filled_many_0.05b')
################################################################################


cross_noDXO <- removeDoubleXO(cross)





## DO NOT USE for final map --> tosses much info
### ONE MARKER PER RAD SITE ####################################################
cross_cons_RADSITES <- thin_by_distortion(cross.cons_0.05,5)
i <- 13 ; plotit(cross_cons_RADSITES,'cons_RADSITES')
cross_many_RADSITES <- thin_by_distortion(cross.many_0.05,5)
i <- 14 ; plotit(cross_many_RADSITES,'many_RADSITES')
################################################################################


## cross.cons_0.01 retains markers along the entire chrm
## cross.cons is the pre-impute object (plt cross.cons 7)
cross.cons_0.01 <- fill.geno(cross.cons, error.prob = 0.01, map.function= 'kosambi')
i <- 9 ; plotit(cross.cons_0.01,'filled_cons_0.01a')
#gt.sub <- geno.table(cross.cons_0.01)
#gt.sub <- gt.sub[!gt.sub$chr %in% c(2,7,13,24),]
#toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
#cross.cons_0.01 <- drop.markers(cross.cons_0.01,toss.sub)
#i <- 10 ; plotit(cross.cons_0.01,'filled_cons_0.01b')


################################################################################
### WRITE THE ABOVE CROSS OBJECT
cross_imp <- cross.cons_0.01
mapfile <- paste0(pop,'_cons_filter_imputed')
filename <- file.path(mpath,mapfile)
write.cross(cross_imp,filestem=filename,format="csv")

cross <- cross.cons
mapfile <- paste0(pop,'_conservative_filter')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
################################################################################

mapfile <- paste0(pop,'_cons_filter_imputed')
cross_imp <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

mapfile <- paste0(pop,'_conservative_filter')
cross_cons <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)




cross_noDXO <- removeDoubleXO(cross)


##################################################################################
##### Sliding window filter to target non-distorted regions ######################
##mean.dist <- function(X, phys, gtb){
##  ind <- which(phys > X[1] & phys < X[2])
##  mean(gtb[ind,'P.value'])
##  }
##
##seg.window <- lapply(1:24,function(Z) {
## gta <- geno.table(cross, chr = Z)
## phys <- as.numeric(gsub(".*:",'',rownames(gta)))
## X1 <- seq(0, max(phys), by = 100000)
## Y1 <- c(seq(0, max(phys), by = 100000)[-1], max(phys))
##
## pos <- apply(cbind(X1,Y1), 1, mean)
## #print(head(gta))
## #pval <- apply(mean.dist, X = X1, Y = Y1, phys = phys, gtb = gta)
## pval <- apply(cbind(X1,Y1),1, mean.dist,phys = phys, gtb = gta)
##cbind(pos,pval)
##})
##
##### CHR2 < 13150000
##seg.window[[2]]
##chr2 <- which(as.numeric(gsub(".*:",'', markernames(cross, 2))) > 13150000)
##chr2 <- markernames(cross, 2)[chr2]
##
##### CHR13 < 12350000
##seg.window[[13]]
##chr13 <- which(as.numeric(gsub(".*:",'', markernames(cross, 13))) > 12350000)
##chr13 <- markernames(cross, 13)[chr13]
##
##################################################################################
##### ALL BUT 2, 13 can be filtered down to 9e-3. Truncates these LGS (see plots)
##gt.sub <- geno.table(cross)
##gt.sub <- gt.sub[!rownames(gt.sub) %in% c(chr2,chr13),]
##toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
##cross <- drop.markers(cross,toss.sub)
##################################################################################
##
## plotted i=1

##################################################################################
#### WRITE THE ABOVE CROSS OBJECT
#mapfile <- paste0(pop,'_',sum(nmar(cross)),'_switched_filtered')
#filename <- file.path(mpath,mapfile)
#write.cross(cross,filestem=filename,format="csv")
##################################################################################

#crossbk <- cross

### ONE MARKER PER RAD SITE ####################################################
#cross <- thin_by_distortion(cross,5)
################################################################################

#i <- 2
#plotit(cross,'all_filters')
#
#
#cross_0.15 <- fill.geno(cross, method="maxmarginal", error.prob = 0.15, min.prob=0.95)
### plotit(cross_0.15,'cross_0.15')
#
#cross_0.15 <- fill.geno(cross, error.prob = 0.15, min.prob=0.95)
#
#
#file <- 'geno.image.0.15.png'
#filename <- file.path('~/public_html/',file)
#png(filename,width=2500,height=500)
#geno.image(cross_0.15)
#dev.off()

################################################################################
### WRITE THE ABOVE CROSS OBJECT
#mapfile <- paste0(pop,'_switched_filtered_thinned_NW')
#filename <- file.path(mpath,mapfile)
#write.cross(cross,filestem=filename,format="csv")
################################################################################

################################################################################
parc2 <- pull.markers(parc1,markernames(cross))
plot_test('par_nbh_filt', width=4000, height = 5000)
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

cross <- cross_imp

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

i <- 99 ; plotit(cross,'imputed_assigned')

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
#fl <- paste0(pop,'_filtered_pvalue_thinned_NW_moved.csv')
#cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
linked_marks <- function(X, LOD = 12, RF = 1){
 crossX <- est.rf(subset(cross,chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}
linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(X)
cross_linked <- pull.markers(cross, unlist(linked))
i <- 95 ; plotit(cross_linked,'cross_linked')
################################################################################

################################################################################
switched_marks <- function(X){
 checkAlleles(subset(cross,chr=X), threshold = 2)
}
switched <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% switched_marks(X)
switched <- switched[!sapply(switched,is.null)]
switched <- do.call(rbind,switched)
cross_linked_switched <- switchAlleles(cross_linked, as.character(switched$marker))
i <- 95 ; plotit(cross_linked_switched,'cross_linked_switched')
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
keep <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% drop.unkinked(i, cross = cross_linked_switched)
cross_linked_switched_drop <- pull.markers(cross_linked_switched, unlist(keep))
i <- 17 ; plotit(cross_linked_switched_drop,'cross_linked_switched_drop')
################################################################################

#################################################################################
### WRITE THE ABOVE CROSS OBJECT
mapfile <- paste0(pop,'_cross_linked_switched_drop')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
#################################################################################

################################################################################
### READ IN THE CROSS
mapfile <- paste0(pop,'_cross_linked_switched_drop.csv')
cross_linked_switched_drop <- read.cross(file = mapfile, format = "csv", dir = mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

# cross_linked_switched_drop_tsp <- tspOrder(cross = cross_linked_switched_drop, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
# i <- 19 ; plotit(cross_linked_switched_drop,'cross_linked_switched_drop_tsp')
#cross_linked_switched_drop <- sim.geno(cross_linked_switched_drop, n.draws=160, error.prob=0.001, map.function="kosambi", stepwidth="fixed")
#cross_linked_switched_drop <- calc.genoprob(cross_linked_switched_drop, error.prob=0.001, map.function="kosambi", stepwidth="fixed")
#sone.bin <- scanone(cross_linked_switched_drop, pheno.col=4, model="bin", method="em")

## MAPPING SUBSET

cross <- cross_linked_switched_drop

dupmar.nexact <- findDupMarkers(cross, exact.only=FALSE, adjacent.only=TRUE) # finds 6 pairs
# one might consider dropping the extra markers totmar(hyper) # 173 markers
cross <- drop.markers(cross, unlist(dupmar.nexact))

i <- 1 ; plotit(cross,'dups_dropped')


################################################################################
cross <- est.rf(cross)
cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
save.image(file.path(mpath,paste0(pop,'_dups_dropped_tsp.rsave')))
i <- 1 ; plotit(cross,'dups_dropped_tsp')

cross.map <- est.map(cross, error.prob=0.001, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross <- replace.map(cross,cross.map)

cross_noDXO <- removeDoubleXO(cross)
i <- 1 ; plotit(cross_noDXO,'dups_dropped_tsp')

cross <- est.rf(cross_noDXO)
cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross.map <- est.map(cross, error.prob=0.001, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross <- replace.map(cross,cross.map)
summary(pull.map(cross))


cross.el <- calc.errorlod(cross,error.prob=0.01)

drops <- top.errorlod(cross.el, cutoff=3)
cross <- drop.markers(cross, unlist(drops))

#cross <- removeDoubleXO(cross)
#cross <- est.rf(cross)
################################################################################
i <- 5 ; plotit(cross,'mapped_tsp')
################################################################################
################################################################################

cross <- drop.markers(cross,c('NW_012234461.1:1059423','NW_012234461.1:562796','NW_012234557.1:172904'))
cross <- drop.markers(cross,c('NW_012234520.1:340401','NW_012234520.1:169476','NW_012224862.1:55723'))
cross <- drop.markers(cross,c('NW_012234558.1:836779','NW_012234558.1:520885','NW_012224981.1:180658','18:387984'))


cross <- est.rf(cross)
cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross.map <- est.map(cross, error.prob=0.001, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross <- replace.map(cross,cross.map)
summary(pull.map(cross))


i <- 5 ; plotit(cross,'mapped_tsp_prunned')


################################################################################
### WRITE IMPUTED AND CORRECTED MAP
mapfile <- paste0(pop,'_',sum(nmar(cross)),'_estrf_no_impute')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
################################################################################

cross.cons_0.01 <- fill.geno(cross, error.prob = 0.01, map.function= 'kosambi')
i <- 18 ; plotit(cross.cons_0.01,'filtered_filled_cons_0.01a')

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
