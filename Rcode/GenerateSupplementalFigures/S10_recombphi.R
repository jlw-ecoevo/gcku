# JLW -2019

library(data.table)
library(ggplot2)
library(dplyr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}


setwd("~/gcku/Shareable/Data")
load("ATGC_GC.Rdata")
load("PhiRecombTest.RData")
load("GCKU.RData")
names(atgc_acc)[2] <- "ATGC" #I named a column ACTG. Nobody's perfect.
names(atgccog_gc)[1] <- "ATGC"

###########################
# Between-Genome Analysis #
###########################

# % Genes Recombining In an ATGC
phi_percent <- phi_df %>% group_by(ATGC) %>% summarize(PercentRecomb=sum(PvalCorrectedBH<0.05)/length(PvalCorrectedBH))

# Connect with Ku and Genomic GC
gc_ku <- gc_ku %>% subset(select=c(Accession2,GC,Ku))
gc_ku$Accession <- gsub("[.].*","",gc_ku$Accession2)
atgc_gcku <- merge.easy(atgc_acc,gc_ku,key="Accession") %>% 
  group_by(ATGC) %>% summarize(GC=mean(GC,na.rm=T),Ku=mean(Ku,na.rm=T))
phi_gcku <- merge.easy(phi_percent,atgc_gcku,key="ATGC")

hist(phi_gcku$Ku)

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="RecombKu.pdf",width=5,height=5)
ggplot(phi_gcku,aes(x=Ku,y=PercentRecomb)) + geom_point() + stat_smooth(method="lm") +
  theme_bw() + xlab("Ku Freq. in ATGC") + ylab("% Recombining Genes") + theme(text = element_text(size=18))
dev.off()

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="RecombGC.pdf",width=5,height=5)
ggplot(phi_gcku,aes(x=PercentRecomb,y=GC)) + geom_point() + stat_smooth(method="lm") +
  theme_bw() + xlab("% Recombining Genes") + ylab("Mean GC of Genomes in ATGC") + theme(text = element_text(size=18))
dev.off()
