# JLW - 2019
# Code to generate violin plot of GC content and tree w/ Ku mapped on

library(ape)
library(phylolm)
library(phytools)
library(vioplot)
library(dplyr)
library(data.table)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

#Load in Tree + Mapping to Data + GC Data for full dataset
setwd("~/gcku/Shareable/Data")
LTP <- read.tree("LTPs128_SSU_tree.newick") #Tree
load("SILVAGCKu.RData") #Match tree to Ku data
load("GCKU.RData") #GC Content + Ku incidence across RefSeq

#Violin plot (All refseq)
setwd("~/gcku/Shareable/Figs")
pdf(file="GCKU_Violin.pdf",width=5,height=5)
vioplot(gc_ku$GC[gc_ku$Ku==0 & !is.na(gc_ku$GC)],gc_ku$GC[gc_ku$Ku==1 & !is.na(gc_ku$GC)],col="gray",names=c("No Ku","Ku"))
title(ylab="GC Content")
dev.off()

#Line up tree w/ Ku data
silva_gcku_full <- silva_gcku_full %>% subset(!is.na(Ku))
species <- silva_gcku_full$NodeNames
pruned.LTP <- drop.tip(LTP,LTP$tip.label[!(LTP$tip.label %in% species)])
silva_gcku_full <- silva_gcku_full %>% subset(NodeNames %in% pruned.LTP$tip.label)
silva_gcku_full <- silva_gcku_full[!duplicated(silva_gcku_full$NodeNames),]
rownames(silva_gcku_full) <- silva_gcku_full$NodeNames
silva_gcku_full$Ku <- as.character(silva_gcku_full$Ku)

#Plot tree
setwd("~/gcku/Shareable/Figs")
pdf(file="Ku_SILVA_TREE.pdf",width=5,height=5)
ku <- as.numeric(as.character(silva_gcku_full$Ku))
names(ku) <- rownames(silva_gcku_full)
ku <- ku[pruned.LTP$tip.label]
plot(pruned.LTP,show.tip.label=FALSE,type="fan")
tiplabels(pch = 21, cex=0.5,col=rgb(1-ku,0,ku,0.25),bg=rgb(1-ku,0,ku,0.25))
dev.off()


