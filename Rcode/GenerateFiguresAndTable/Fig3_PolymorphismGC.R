# JLW - 2019
# Use polymorphism data (from ATGC database) to estimate mutational bias, calculate expected GC
# Polymorphism GC content calculated using GC_Polymorphism.R 

library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/gcku/Shareable/Data")
load("GCKU.RData")
gc_ku <- gc_ku %>% subset(select=c(Accession2,Ku))

###################
# GC at All Sites #
###################

#load polymorphism data
setwd("~/gcku/Shareable/Data")
load("PolyGCBiasRate.RData")
load("GCKU.RData")

#Subset only to pairs of sequence that have recently diverged, but still have at least 5 informative sites
poly_df <- poly_df %>% subset(Divergence<=0.01) %>% subset(nInformative>4)

#Calculate expected GC content (eqn. from Long et al. 2018, Nat. Eco. Evo.)
poly_df$GCpoly <- 1/(1+poly_df$m)

#Total informative Sites
nrow(poly_df)
sum(poly_df$nInformative)

#Attach Ku data
poly_df$Seq <- unlist(lapply(strsplit(poly_df$Seq,split="[|]"),"[",6))
poly_df$Accession2 <- paste("GCF",gsub("^.*GCF_","",poly_df$Seq),sep="_")
poly_df <- merge.easy(poly_df,gc_ku,key="Accession2") %>% subset(!is.na(Ku))

table(poly_df$Ku)

p1 <- ggplot(poly_df,aes(y=GCseq,x=GCpoly,color=factor(Ku))) + geom_point(alpha=0.25) + 
  stat_smooth(method="lm") + ylim(0,1) + xlim(0,1) + geom_abline(slope=1,intercept=0,lty=2) +
  xlab("Expected GC") + ylab("Actual GC")+ theme_bw()+
  theme(text = element_text(size=18),legend.position = "none")

setwd("~/gcku/Shareable/Figs")
pdf("ATGC_GCpolyvsGCseqvsKu.pdf",width=7,height=5)
p1
dev.off()

###################################
# GC at Fourfold Degenerate Sites #
###################################

#load polymorphism data
setwd("~/gcku/Shareable/Data")
load("PolyGCBiasRate.RData")
poly_df <- poly_df %>% subset(Divergence<=0.01)
poly_df$GCpoly <- 1/(1+poly_df$m)
polyAll_df <- poly_df %>% subset(select=c(Seq,GCpoly,nInformative))
# Add GC4 data
setwd("~/gcku/Shareable/Data")
load("PolyGC4BiasRate_GC4expected.RData")
poly_df <- poly_df %>% subset(Divergence<=0.01) %>% subset(select=-c(nInformative))
poly_df$GC4poly <- 1/(1+poly_df$m)
poly_df <- merge.easy(poly_df,polyAll_df,key="Seq") %>% subset(nInformative>4)

#Attach Ku data
poly_df$Seq <- unlist(lapply(strsplit(poly_df$Seq,split="[|]"),"[",6))
poly_df$Accession2 <- paste("GCF",gsub("^.*GCF_","",poly_df$Seq),sep="_")
poly_df <- merge.easy(poly_df,gc_ku,key="Accession2") %>% subset(!is.na(Ku))

p2 <- ggplot(poly_df,aes(y=GC4seq,x=GCpoly,color=factor(Ku))) + geom_point(alpha=0.25) + 
  stat_smooth(method="lm") + ylim(0,1) + xlim(0,1) + geom_abline(slope=1,intercept=0,lty=2) +
  xlab("Expected GC") + ylab("Actual GC at Fourfold Degenerate Sites")+ theme_bw()+
  theme(text = element_text(size=18),legend.position = c(0.875,0.12))+ scale_color_discrete(name = "",labels=c("No Ku","With Ku"))

setwd("~/gcku/Shareable/Figs")
pdf("ATGC_GCpolyvsGC4seqvsKu_GC4vAllexp.pdf",width=7,height=5)
p2
dev.off()

setwd("~/gcku/Shareable/Figs")
pdf("Fig3_polymorphism.pdf",width=12,height=5)
ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"),hjust=-3.5,vjust=2)
dev.off()
