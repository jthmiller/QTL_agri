## Directories
plotdir <- file.path(basedir, "rQTL/plots")
indpops <- file.path(basedir, "plinkfiles/ind.pops")
popdir <- file.path(basedir, "rQTL", pop, "REMAPS")
qtldir <- file.path(basedir, "rQTL/remap_out")
errfile <- file.path(qtldir, "genotyping_error_rate.txt")
dirso <- "/home/jmiller1/QTL_agri/data"

## Funtions for processing rQTL map data
source(file.path("MAP/R/source_file.R"))
# source(file.path(basedir, 'rQTL/scripts/QTL_remap/QTL/model_source_file.R'))

## Libraries
flib <- "/share/apps/rmodules"
fpacks <- c("devtools", "httr", "RColorBrewer", "qtl")
lapply(fpacks, require, character.only = TRUE, lib.loc = flib)

mylib <- "/home/jmiller1/R/x86_64-pc-linux-gnu-library/3.5"
##mpacks <- c("qtl", "foreach", "qtl2", "qtlTools", "plyr")
mpacks <- c("qtl", "qtl2", "qtlTools", "plyr")
#mpacks <- c("qtl", "foreach", "qtl2", "qtlTools", "doParallel", "plyr")
lapply(mpacks, require, character.only = TRUE, lib.loc = mylib)

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
