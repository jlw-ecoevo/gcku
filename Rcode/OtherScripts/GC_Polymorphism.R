# JLW - 2019
# Calculate GC content of ATGC polymorphisms
# ***Must first download and unpack ATGC database (Download_Unpack_ATGCDatabase.R) - but output of this script already available in /Data folder (for downstream analyses)

library(data.table)
library(seqinr)
library(ggplot2)
library(dplyr)
library(vioplot)

compareSeq <- function(x,a,b){sum(x[[a]]==x[[b]])/length(x[[a]])}

compareAllSeq <- function(x){
  seq_dist <- matrix(0,nrow=length(x),ncol=length(x))
  for(i in 1:(length(x)-1)){
    for(j in (i+1):length(x)){
      seq_dist[i,j] <- compareSeq(x,i,j)
    }
  }
  seq_dist <- seq_dist + t(seq_dist)
  diag(seq_dist) <- NA
  return(seq_dist)
}

polyGC <- function(x,a,b,outgroup){
  p1 <- which(x[[a]]!=x[[outgroup]])
  p2 <- which(x[[b]]!=x[[outgroup]])
  pshared <- intersect(p1,p2)
  p1 <- setdiff(p1,pshared)
  p2 <- setdiff(p2,pshared)
  
  nGCtoAT <- sum((x[[a]][p1] %in% c("a","t")) & (x[[b]][p1] %in% c("c","g")))
  nATtoGC <- sum((x[[a]][p1] %in% c("c","g")) & (x[[b]][p1] %in% c("a","t")))
  nAT <- sum(x[[outgroup]][-pshared] %in% c("a","t"))
  nGC <- sum(x[[outgroup]][-pshared] %in% c("c","g"))
  
  if(length(p1)==0){
    return(c(NA,nGCtoAT+nATtoGC))
  } else {
    #m <- (nATtoGC/nAT)/(nGCtoAT/nGC)
    m <- (nGCtoAT/nGC)/(nATtoGC/nAT)
    return(c(m,nGCtoAT+nATtoGC))
  }
}

allPolyGC <- function(x,outgroup,pairs){
  seq_dist <- matrix(nrow=nrow(pairs),ncol=2)
  for(i in 1:nrow(pairs)){
    seq_dist[i,] <- polyGC(x,pairs[i,1],pairs[i,2],outgroup)
  }
  return(seq_dist)
}

polyGCBiasVsDist <- function(atgc_id){
  print(atgc_id)
  foldername <- paste("~/gcku/Shareable/ATGC/ATGC_unpacked/ATGC",
                      atgc_id,"/",sep="")
  #x <- read.fasta(paste(foldername,"ATGC",atgc_id,"snp.fasta",sep=""))
  x_full <- read.fasta(paste(foldername,"ATGC",atgc_id,".fasta",sep=""))
  #
  xsim <- compareAllSeq(x_full)
  outgroup <- which.min(colMeans(xsim,na.rm=T))
  pairs <- cbind(1:nrow(xsim),apply(xsim,1,which.max))
  pairs <- pairs[-which(pairs[,1]==outgroup | pairs[,2]==outgroup),]
  poly_GC <- allPolyGC(x_full,outgroup,pairs)
  seq_GC <- unlist(lapply(x_full,GC))[pairs[,1]]
  dist_Pairs <- xsim[pairs]
  #
  return(data.frame(Seq=names(x_full)[pairs[,1]],
                    Divergence=1-dist_Pairs,
                    m=poly_GC[,1],
                    nInformative=poly_GC[,2],
                    GCseq=seq_GC,
                    stringsAsFactors = F))
}

tryPolyGCBiasVsDist <- function(atgc_id){return(try(polyGCBiasVsDist(atgc_id)))}
  
setwd("~/gcku/Shareable/Data/")
atgc_metadata <- read.csv("atgcdata.csv")
to_poly <- setdiff(gsub("ATGC","",as.character(atgc_metadata$ATGC.[atgc_metadata$numgenomes>=4])),c("071","078"))
# 
poly <- lapply(to_poly,tryPolyGCBiasVsDist)
x <- lapply(poly,inherits,"try-error")
poly_df <- as.data.frame(do.call("rbind",poly[which(!unlist(x))]),stringsAsFactors=F)

setwd("~/gcku/Shareable/Data/")
save(poly_df,file="PolyGCBiasRate.RData")
