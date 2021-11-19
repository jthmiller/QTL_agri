#!/bin/R

#### MAKE QTL LIST

qtls <- read.table('QTLs-Table 1.tsv',stringsAsFactors=F, sep='\t')
qtls <- qtls[!is.na(qtls$V1),]
qtls <- qtls[!duplicated(qtls$V1),]
rownames(qtls) <- qtls$V1
qtls$V14 <- gsub(" ","",qtls$V14)

s1 <- read.table('~/Desktop/PTC_uploads/S_tables/TABLE_S1_JM.csv',stringsAsFactors=F, sep=',',header=T)



pheotypes <- rowSums(!is.na(qtls[,8:13]))



chr_pos <- do.call(rbind,lapply(s1$loc, function(X){ cbind( gsub(":.*","",X), as.numeric(gsub(".*:","",X)))  } ))


chr <- lapply(s1$loc, function(X){  gsub(":.*","",X)  })
pos <- lapply(s1$loc, function(X){  as.numeric(gsub(".*:","",X))  })

mark_chr <- lapply( chr, function(X){ rownames(qtls) [which(qtls$V14 == X)] })
near_marker_dis <- mapply(function(X,Y){ min(abs(qtls[which(qtls$V14 == X),'V5'] - as.numeric(Y)),na.rm=T) }, X=chr, Y=pos)
near_marker_ind <- mapply(function(X,Y){ which.min(abs(qtls[which(qtls$V14 == X),'V5'] - as.numeric(Y))) }, X=chr, Y=pos)
near_marker <- mapply(function(X,Y){ X[Y] } , X = mark_chr, Y = near_marker_ind)
near_marker <- lapply(near_marker, function(X){ if(length(X) < 1 ) { NA } else { X } })

mark_chr_qtl <- lapply( chr, function(X){ rownames(qtls) [which(qtls$V14 == X & pheotypes == 1)] })
near_qtl_dis <- mapply(function(X,Y){ min(abs(qtls[which(qtls$V14 == X & pheotypes == 1 ),'V5'] - as.numeric(Y)),na.rm=T) }, X=chr, Y=pos)
near_qtl_ind <- mapply(function(X,Y){ which.min(abs(qtls[which(qtls$V14 == X & pheotypes == 1 ),'V5'] - as.numeric(Y))) }, X=chr, Y=pos)
near_qtl <-  mapply(function(X,Y){ X[Y] } , X = mark_chr_qtl, Y = near_qtl_ind)
near_qtl <-  lapply(near_qtl, function(X){ if(length(X) < 1 ) { NA } else { X } })

ptc_qtl <- data.frame(
  chrn=qtls[unlist(near_marker),'V3'],
  chr=unlist(chr),
  pos=unlist(pos),
  near_marker=unlist(near_marker),
  near_marker_dis=unlist(near_marker_dis),
  near_qtl=unlist(near_qtl),
  near_qtl_dis=unlist(near_qtl_dis))

rownames(ptc_qtl) <- s1$loc



write.table(ptc_qtl,"~/Des
