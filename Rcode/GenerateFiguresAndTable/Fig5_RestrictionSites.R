# JLW - 2019
# Code to generate panels in figure 4



library(ggplot2)
library(seqinr) 
library(ggpubr)

# Panel A - GC content of RE recognition sequences
# REs summarized nicely in FlankGCRfast_flankdist1.txt (original enzyme list in REBASEGenomeEnzymes.RData but not with geenomic GC for comparison)
setwd("~/gcku/Shareable/Data")
x <- read.delim("~/gcku/Data/FlankGCRfast_flankdist1.txt",header = F) %>% subset(!is.na(V4)) 
x$ID <- paste(x$V1,x$V3,sep="_")
names(x) <- c("File","GCseq","RecSeq","Flank1","nSites","ID")

#Calculate GC of rec seqs and average for each genome
x$RecGC <- unlist(lapply(strsplit(as.character(x$RecSeq),split=""),GC))
y <- x %>% group_by(File) %>% summarize(SeqGC=mean(GCseq),RecGC=mean(RecGC))

p1 <- ggplot(y,aes(x=SeqGC,y=RecGC)) + geom_point(alpha=1) + stat_smooth(method="lm",col="red") + 
  ylim(0,1) + xlim(0,1) + geom_abline(intercept=0,slope=1,lty=2) + 
  xlab("Genomic GC") + ylab("Mean Recognition GC") +theme_bw() + 
  theme(text = element_text(size=18))

setwd("~/gcku/Shareable/Figs")
pdf("RRecGCVsGenomicGC.pdf",width=7,height=5)
p1
dev.off()


#Panel B - GC content near restriction sites (data generated using Summarize_Restriction_Site_FLanks.R)

setwd("~/gcku/Shareable/Data")
load("GCRMFlanks.RData")

p2 <- ggplot(distdiff_df,aes(x=Dist,y=GCdiff)) + geom_point(size=1) +
  geom_errorbar((aes(ymin=CIlower,ymax=CIupper)),alpha=0.25) +# ylim(-0.025,0.1)+ 
  xlab("Distance from Recognition Site (bp)")+ theme_bw() + geom_abline(slope=0,intercept=0,lty=2) + 
  ylab(expression(paste("Excess GC (",GC[Flank]-GC[Null],")"))) + geom_line(lwd=0.25) + 
  theme(text = element_text(size=18)) + xlim(0,200) + stat_smooth(col="red")

setwd("~/gcku/Shareable/Figs")
pdf("GCFlank_AT75_WithNull.pdf",width=7,height=5,pointsize = 10)
p2
dev.off()


setwd("~/gcku/Shareable/Figs/")
pdf(file="Fig5_rm.pdf",width=14,height=5)
ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"),hjust=-4,vjust=2)
dev.off()

