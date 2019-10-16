# JLW -2019
# Fig 4 - Recombination along the genome
# The data to run this analysis is already in /Data but if you would like to generate the recombination data yourself use the scripts in /OtherScripts (Unpack_ATGCDatabase_forPhiPack.R, PhiPack_Run.R, PhiOut_Process.R, ATGCCOGpair_GC.R)


library(data.table)
library(ggplot2)
library(dplyr)
library(vioplot)
library(ggpubr)

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


# Connect with Ku and Genomic GC
gc_ku <- gc_ku %>% subset(select=c(Accession2,GC,Ku))
gc_ku$Accession <- gsub("[.].*","",gc_ku$Accession2)
atgc_gcku <- merge.easy(atgc_acc,gc_ku,key="Accession") %>% 
  group_by(ATGC) %>% summarize(GC=mean(GC,na.rm=T),Ku=mean(Ku,na.rm=T))

atgccog_gc$ATGC.COG <- paste(atgccog_gc$ATGC,atgccog_gc$COG,sep=".")
phi_df$ATGC.COG <- paste(phi_df$ATGC,phi_df$COG,sep=".")

phi_gc  <- merge.easy(phi_df,atgccog_gc,key="ATGC.COG")

pgc <- phi_gc %>% group_by(ATGC.x) %>% summarise(GCs=mean(GC[PvalCorrectedBH<0.05]),
                                               GCns=mean(GC[PvalCorrectedBH>=0.05]))
names(pgc)[1] <- "ATGC"

#Remove ATGCs w/out significant recombining genes
pgc <- pgc %>% subset(!is.na(GCs))

#add recomb, ku
pgc <- merge.easy(pgc,atgc_gcku,key="ATGC") 

p1 <- ggplot(pgc,aes(x=GCns,y=GCs,color=Ku)) + geom_point() + 
  geom_abline(intercept=0,slope=1,lty=2) + theme_bw() +
  xlab("Mean GC Non-Recombining Genes") + ylab("Mean GC Recombining Genes") +
  theme(text = element_text(size=14)) + scale_color_continuous(name="Ku Freq.\n in Cluster")

setwd("~/gcku/Shareable/Figs/")
pdf(file="GCRecomb_VsKu.pdf",width=5,height=7)
p1
dev.off()

pgc$Diff <- pgc$GCs-pgc$GCns
vioplot(pgc$Diff,col="grey")
abline(h=0,lty=2)

pgc$Ku2 <- "No Ku"
pgc$Ku2[pgc$Ku>0] <- "Ku"
pgc$Ku2 <- as.factor(pgc$Ku2)
p2 <- ggplot(pgc,aes(y=Diff,x=Ku2)) + geom_violin(fill="gray") + geom_boxplot(width=0.1) + 
  theme_bw() + xlab("") + ylab("(Mean GC Recombining Genes) - (Mean GC Non-Recombining Genes)") +
  theme(text = element_text(size=14)) + geom_abline(slope=0,intercept=0,lty=2)

setwd("~/gcku/Shareable/Figs/")
pdf(file="GCRecomb_VsKu_Violin.pdf",width=5,height=7)
p2
dev.off()

setwd("~/gcku/Shareable/Figs/")
pdf(file="Fig4_recomb.pdf",width=10,height=7)
ggarrange(p1,p2,ncol=2,widths=c(2,1),labels=c("(a)","(b)"),hjust=-3,vjust=2)
dev.off()

t.test(pgc$GCs,pgc$GCns,paired = T)

length(pgc$Diff[pgc$Ku==0])
length(pgc$Diff[pgc$Ku>0])
