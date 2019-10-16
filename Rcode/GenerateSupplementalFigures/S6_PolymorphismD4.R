# JLW - 2019
# Use D4 polymorphism data (from ATGC database) to estimate mutational bias, calculate expected GC
# Polymorphism GC content calculated using GC_Polymorphism_D4.R 

library(data.table)
library(ggplot2)
library(dplyr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/gcku/Shareable/Data")
load("GCKU.RData")
gc_ku <- gc_ku %>% subset(select=c(Accession2,Ku))

setwd("~/gcku/Shareable/Data")
load("PolyGC4BiasRate_GC4expected.RData")
poly_df <- poly_df %>% subset(Divergence<=0.01) %>% subset(nInformative>4)
poly_df$GC4poly <- 1/(1+poly_df$m)

#Attach Ku data
poly_df$Seq <- unlist(lapply(strsplit(poly_df$Seq,split="[|]"),"[",6))
poly_df$Accession2 <- paste("GCF",gsub("^.*GCF_","",poly_df$Seq),sep="_")
poly_df <- merge.easy(poly_df,gc_ku,key="Accession2") %>% subset(!is.na(Ku))

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("ATGC_GC4polyvsGCseqvsKu.pdf",width=7,height=5)
ggplot(poly_df,aes(y=GCseq,x=GC4poly,color=factor(Ku))) + geom_point(alpha=0.25) + 
  stat_smooth(method="lm") + ylim(0,1) + xlim(0,1) + geom_abline(slope=1,intercept=0,lty=2) +
  xlab("Expected GC at Fourfold Degenerate Sites") + ylab("Actual GC")+ theme_bw()+
  theme(text = element_text(size=18))+  scale_color_discrete(name = "",labels=c("No Ku","With Ku"))
dev.off()

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("ATGC_GC4polyvsGC4seqvsKu.pdf",width=7,height=5)
ggplot(poly_df,aes(y=GC4seq,x=GC4poly,color=factor(Ku))) + geom_point(alpha=0.25) + 
  stat_smooth(method="lm") + ylim(0,1) + xlim(0,1) + geom_abline(slope=1,intercept=0,lty=2) +
  xlab("Expected GC at Fourfold Degenerate Sites") + ylab("Actual GC at Fourfold Degenerate Sites")+ theme_bw()+
  theme(text = element_text(size=18))+ scale_color_discrete(name = "",labels=c("No Ku","With Ku"))
dev.off()



