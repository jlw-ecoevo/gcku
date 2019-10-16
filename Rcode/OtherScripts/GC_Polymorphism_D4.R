# JLW - 2019
# Calculate GC content of ATGC polymorphisms at fourfold degenerate sites
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
  xb <- x[[b]][x[[a]]!="-"]
  xo <- x[[outgroup]][x[[a]]!="-"]
  xa <- x[[a]][x[[a]]!="-"]
  xa3 <- xa[seq(3,length(xa),3)]
  xb3 <- xb[seq(3,length(xb),3)]
  xo3 <- xo[seq(3,length(xo),3)]
  degenerate <- D4(xa)
  xad <- xa3[degenerate]
  xbd <- xb3[degenerate]
  xod <- xo3[degenerate]
  
  p1 <- which(xad!=xod)
  p2 <- which(xbd!=xod)
  pshared <- intersect(p1,p2)
  p1 <- setdiff(p1,pshared)
  p2 <- setdiff(p2,pshared)
  
  nGCtoAT <- sum((xad[p1] %in% c("a","t")) & (xbd[p1] %in% c("c","g")))
  nATtoGC <- sum((xad[p1] %in% c("c","g")) & (xbd[p1] %in% c("a","t")))
  nAT <- sum(xod[-pshared] %in% c("a","t"))
  nGC <- sum(xod[-pshared] %in% c("c","g"))

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

GC4 <- function(x){
  x <- x[x!="-"]
  x1 <- x[seq(1,length(x),3)]
  x2 <- x[seq(2,length(x),3)]
  x3 <- x[seq(3,length(x),3)]
  
  D4_sites <- x3[(x1=="g" & x2=="c") |
    (x1=="g" & x2=="g") |
    (x1=="g" & x2=="t") |
    (x1=="c" & x2=="g") |
    (x1=="c" & x2=="c") |
    (x1=="c" & x2=="t") |
    (x1=="a" & x2=="c") |
    (x1=="t" & x2=="c")]
  
  return(GC(D4_sites))
}

D4 <- function(x){
  #needs gapless sequence
  x1 <- x[seq(1,length(x),3)]
  x2 <- x[seq(2,length(x),3)]
  x3 <- x[seq(3,length(x),3)]
  
  D4 <- which((x1=="g" & x2=="c") |
                   (x1=="g" & x2=="g") |
                   (x1=="g" & x2=="t") |
                   (x1=="c" & x2=="g") |
                   (x1=="c" & x2=="c") |
                   (x1=="c" & x2=="t") |
                   (x1=="a" & x2=="c") |
                   (x1=="t" & x2=="c"))
  
  return(D4)
}

polyGCBiasVsDist <- function(atgc_id){
  print(atgc_id)
  foldername <- paste("~/gcku/Shareable/ATGC/ATGC_unpacked/ATGC",
                      atgc_id,"/",sep="")
  # x <- read.fasta(paste(foldername,"ATGC",atgc_id,"snp.fasta",sep=""))
  x_full <- read.fasta(paste(foldername,"ATGC",atgc_id,".fasta",sep=""))
  #
  xsim <- compareAllSeq(x_full)
  outgroup <- which.min(colMeans(xsim,na.rm=T))
  pairs <- cbind(1:nrow(xsim),apply(xsim,1,which.max))
  pairs <- pairs[-which(pairs[,1]==outgroup | pairs[,2]==outgroup),]
  poly_GC <- allPolyGC(x_full,outgroup,pairs)
  seq_GC <- unlist(lapply(x_full,GC))[pairs[,1]]
  seq_GC4 <- unlist(lapply(x_full,GC4))[pairs[,1]]
  dist_Pairs <- xsim[pairs]
  #
  return(data.frame(Seq=names(x_full)[pairs[,1]],
                    Divergence=1-dist_Pairs,
                    m=poly_GC[,1],
                    nInformative=poly_GC[,2],
                    GCseq=seq_GC,
                    GC4seq=seq_GC4,
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
save(poly_df,file="PolyGC4BiasRate_GC4expected.RData")


