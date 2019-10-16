#JLW -2019
#Finds GC content of bases flanking fake restriction sites on genome (generated with Generate_Null_RestrictionEnzymes.R)
#Written to be run on multiple cores (change mc.cores argument to mcmapply), taking the distance away from sites as an argument
# ***Before running you must download the assemblies from RefSeq listed in enzymetype_rm$File (FileRMRec_EnzymeType_NullAltSeq.RData)


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (for flank distance)", call.=FALSE)
} else {
  flank_dist <- as.numeric(as.character(args[1]))
}

library(dplyr)
library(seqinr)
library(stringr)
library(stringi)
library(parallel)

setwd("~/gcku/Shareable/Data")
load("FileRMRec_EnzymeType_NullAltSeq.RData")

illegalchars <- paste("[",paste(setdiff(unique(unlist(strsplit(enzymetype_rm$PredTargetSeq,split=""))),letters),collapse="]|["),"]",sep="")
enzymetype_rm <- enzymetype_rm[-grep(illegalchars,enzymetype_rm$PredTargetSeq),]


frm <- enzymetype_rm  %>% subset(EnzymeType=="R") %>% group_by(File) %>% 
  summarise(RecSeqs=paste(PredTargetSeq,collapse = ","))

recToRegex <- function(x){
  return(x %>% gsub("n",".",.) %>%
           gsub("b","[cgt]",.) %>%
           gsub("d","[agt]",.) %>%
           gsub("h","[act]",.) %>%
           gsub("k","[gt]",.) %>%
           gsub("m","[ac]",.) %>%
           gsub("r","[ag]",.) %>%
           gsub("s","[cg]",.) %>%
           gsub("v","[acg]",.) %>%
           gsub("w","[at]",.) %>%
           gsub("y","[ct]",.))
}

flankGC <- function(recognition_seq,seq,flankdist){
  sites <- str_locate_all(paste(seq,collapse=""),recToRegex(recognition_seq))
  
  first_flank <- sites[[1]][,1]-flankdist
  second_flank <- sites[[1]][,2]+flankdist
  first_flank <- first_flank[first_flank>0]
  second_flank <- second_flank[second_flank<=length(seq)]
  seq <- seq[unlist(c(first_flank,second_flank))]
  flank_gc <- try(GC(seq))
  if(inherits(flank_gc,"try-error")){
    flank_gc <- NA
  }
  
  return(c(recognition_seq, flank_gc, length(sites[[1]][,2])))
}

getFlanks <- function(file,rec_seqs,flankdist){
  seq <- unlist(read.fasta(file))
  seq_gc <- GC(seq)
  recseq_list <- unlist(strsplit(rec_seqs,split=","))
  xgc <- unlist(lapply(recseq_list,flankGC,seq=seq,flankdist=flankdist))
  return(list(File=file,GC=seq_gc,FlankGC=xgc))
}


fileFlanks <- function(file,recseq,flankdist){
  print(file)
  x <- try(getFlanks(file,recseq,flankdist))
  if(!inherits(x,"try-error")){
    write(paste(x$File,
                x$GC,
                x$FlankGC[seq(1,length(x$FlankGC),3)],
                x$FlankGC[seq(2,length(x$FlankGC),3)],
                x$FlankGC[seq(3,length(x$FlankGC),3)],
                sep="\t"),file=paste("FlankGCRfast_NULLAltSeq_flankdist",flankdist,".txt",sep=""),append=T)
  }
  return(x)
}

x <- mcmapply(fileFlanks,frm$File,frm$RecSeqs,MoreArgs=list(flankdist=flank_dist),mc.cores=19)
save(x,file=paste("FlankGCRfast_NULLAltSeq_flankdist",flank_dist,".RData",sep=""))






