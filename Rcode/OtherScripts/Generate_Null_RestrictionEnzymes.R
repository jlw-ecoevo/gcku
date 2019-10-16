#JLW -2019
# Takes known restriction target sequences and scrambles to generate "fake" enzymes for null 

library(stringi)
library(dplyr)


setwd("~/gcku/Shareable/Data")
load("REBASEGenomeEnzymes.RData")

enzymetype_rm <- enzymetype_rm %>% subset(enzymetype_rm$PredTargetSeq!="cc")

enzymetype_rm$PredTargetSeqOld <- enzymetype_rm$PredTargetSeq
enzymetype_rm$PredTargetSeq <- stri_rand_shuffle(enzymetype_rm$PredTargetSeq)

enzymetype_rm$PredTargetSeq==enzymetype_rm$PredTargetSeqOld

while(sum(enzymetype_rm$PredTargetSeq==enzymetype_rm$PredTargetSeqOld)>0){
  print(enzymetype_rm$PredTargetSeq[enzymetype_rm$PredTargetSeq==enzymetype_rm$PredTargetSeqOld])
  enzymetype_rm$PredTargetSeq[enzymetype_rm$PredTargetSeq==enzymetype_rm$PredTargetSeqOld] <- 
    stri_rand_shuffle(enzymetype_rm$PredTargetSeq[enzymetype_rm$PredTargetSeq==enzymetype_rm$PredTargetSeqOld])
}

save(enzymetype_rm,file="FileRMRec_EnzymeType_NullAltSeq.RData")


