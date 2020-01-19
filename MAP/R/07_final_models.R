#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

pop <- 'ELR'

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

# single scans
load(file.path(mpath,'single_scans.elr.rsave'))
sw_elr <- full.norm.add_only
pni_elr <- perms.norm.imp
pbe_elr <- perms.bin.em
sbe_elr <- scan.bin.em
sbm_elr <- scan.bin.mr
elr_cross <- cross

load(file.path(mpath,'single_scans.nbh.rsave'))
sw_nbh <- full.norm.add_only
pni_nbh <- perms.norm.imp
pbe_nbh <- perms.bin.em
sbe_nbh <- scan.bin.em
sbm_nbh <- scan.bin.mr
nbh_cross <- cross


load(file.path(mpath,'scantwo.scans.nbh.short.rsave'))
full.norm.add_only
sc2_normal_imp
sc2_normal_imp_perms

load(file.path(mpath,'single_scans.new.rsave'))
sw_new <- full.norm.add_only
pni_new <- perms.norm.imp
pbe_new <- perms.bin.em
sbe_new <- scan.bin.em
sbm_new <- scan.bin.mr
new_cross <- cross


load(file.path(mpath,'single_scans.brp.rsave'))
sw_brp <- full.norm.add_only
pni_brp <- perms.norm.imp
pbe_brp <- perms.bin.em
sbe_brp <- scan.bin.em
sbm_brp <- scan.bin.mr
brp_cross <- cross

get_phenos <- function(crs,pheno){
 index <- as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
 subset(crs,ind=index)
}

pheno_ind <- function(crs,pheno){
 as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
}



pdf("/home/jmiller1/public_html/nbh_18.pdf", width=20)
plot(sbe_nbh)
dev.off()

pdf("/home/jmiller1/public_html/elr_bin.pdf", width=20)
plot(sbe_elr)
dev.off()




############################################################
## NBH
summary(full.norm.imp)
      name chr pos n.gen

Q1  2@43.0   2  43     3
Q2  2@87.0   2  87     3
Q3  3@31.0   3  31     3
Q4  3@38.0   3  38     3
Q5 13@31.0  13  31     3
Q6 18@51.0  18  51     3
Q7 19@18.0  19  18     3

Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q3:Q5 + Q1:Q6

pLOD:  49.014

fitqtl(gg_step2, pheno.col=5, qtl=full.norm.imp, method="imp",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q3:Q5 + Q1:Q6")

############################################################

summary(full.norm.hk)
 name chr pos n.gen
Q1  2@103.0   2 103     3
Q2   3@17.0   3  17     3
Q3   3@84.0   3  84     3
Q4  11@54.0  11  54     3
Q5  13@28.0  13  28     3
Q6 15@114.0  15 114     3
Q7  18@41.0  18  41     3
Q8  19@49.0  19  49     3

  Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q7 + Q2:Q5 + Q4:Q8 +
                Q3:Q6

  pLOD:  30.266

fitqtl(gg_step2, pheno.col=5, qtl=full.norm.hk, method="hk",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q7 + Q2:Q5 + Q4:Q8 + Q3:Q6")

############################################################
