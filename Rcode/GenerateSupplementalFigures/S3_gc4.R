# JLW - 2019


library(dplyr)
library(data.table)
library(vioplot)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/gcku/Shareable/Data")
load(gc4,file="gc4.RData")
load("GCKU.RData")
gc_ku <- merge.easy(gc_ku,gc4,key="Accession2")

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="GC4KU_Violin.pdf",width=5,height=5)
vioplot(gc_ku$GC4[gc_ku$Ku==0 & !is.na(gc_ku$GC4)],gc_ku$GC4[gc_ku$Ku==1 & !is.na(gc_ku$GC4)],col="gray",names=c("No Ku","Ku"))
title(ylab="GC4 Content")
dev.off()


