#!/bin/R
### first run combine pops for multi-pop cross objects

pop <- 'NBH'

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')

mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

load(file.path(mpath,'supplemental_plot_env.rsave'))

erp <- 0.0025

cross_NBH <- sim.geno(cross_NBH, step=1, error.prob=erp, off.end=5, map.function="kosambi", n.draws=100, stepwidth="fixed")
cross_NBH <- calc.genoprob(cross_NBH, step=1, error.prob=erp, off.end=5, map.function="kosambi", stepwidth="fixed")
cross_grid <- reduce2grid(cross_NBH)


##plot.scanone, use incl.markers=FALSE
