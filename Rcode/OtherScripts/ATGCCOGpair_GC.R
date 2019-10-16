# JLW - 2019
# finds average GC content of alignment for each cluster-gene pair

library(seqinr)
library(dplyr)


getSeqGC <- function(file){
  seqs <- read.fasta(file)
  gcseqs <- lapply(seqs,GC)
  ATGC.COG <- names(gcseqs) %>% strsplit(split="[|]") %>% lapply("[",2) %>% unlist()
  Spp.Acc <- names(gcseqs) %>% strsplit(split="[|]") %>% lapply("[",6) %>% unlist()
  ATGC <- ATGC.COG %>% gsub(pattern="[.].*",replace="")
  COG <- ATGC.COG %>% gsub(pattern=".*[.]",replace="")
  Spp <- Spp.Acc %>% gsub(pattern="[.].*",replace="")
  Acc <- Spp.Acc %>% strsplit(split="[.]") %>% lapply("[",2) %>% unlist()
  return(data.frame(ACTG=ATGC,
                    COG=COG,
                    Species=Spp,
                    Accession=Acc,
                    GC=unlist(gcseqs),
                    stringsAsFactors = F))
}


#system("cd ~/gcku/Shareable/ATGCATGC_nucleotidealignments; ls *.fa > ../ATGC_alginment_list.txt")
files <- readLines("~/gcku/Shareable/ATGCATGC_alginment_list.txt") %>% paste0("~/gcku/Shareable/ATGCATGC_nucleotidealignments/",.)


atgc_gc <- data.frame(ACTG=character(),
                       COG=character(),
                       Species=character(),
                       Accession=character(),
                       GC=numeric(),
                       stringsAsFactors = F)
for(file in files){
  print(file)
  atgc_gc <- rbind(atgc_gc,
                   getSeqGC(file))
}

setwd("~/gcku/Shareable/Data")
save(atgc_gc,file="ATGC_GC_Full.Rdata")

atgc_acc <- atgc_gc %>% subset(select=c(Accession,ACTG)) %>% unique()
atgccog_gc <- atgc_gc %>% group_by(ACTG,COG) %>% summarize(GC=mean(GC))

setwd("~/gcku/Shareable/Data")
save(atgc_acc,atgccog_gc,file="ATGC_GC.Rdata")





