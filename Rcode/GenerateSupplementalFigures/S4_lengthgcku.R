
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

set.seed(20011)

setwd("~/gcku/Shareable/Data")
load("SILVAGCKu.RData")
load("refseq_all_genome_lengths.RData")
load("taxonomy.RData")

silva_gcku_full <- merge.easy(silva_gcku_full,gl_df,key="Accession2")
silva_gcku_full <- silva_gcku_full %>% subset(!is.na(GC) & !is.na(GenomeLength))
silva_gcku_full$GenomeLength <- as.numeric(silva_gcku_full$GenomeLength)
x <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% summarize(Ku=mean(as.numeric(Ku)))
hist(x$Ku)
sum(x$Ku <1 & x$Ku>0)/nrow(x)
x <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% summarize(VGC=var(GC))
max(x$VGC,na.rm=T)
hist(x$VGC,breaks=100)
silva_gcku <- silva_gcku %>% group_by(SpeciesSILVA) %>% sample_n(1)
silva_gcku_full <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% sample_n(1)

p1 <- ggplot(silva_gcku_full[silva_gcku_full$Ku==1,],aes(x=GenomeLength,y=GC)) + 
  geom_point(color=rgb(0,0,1,0.25))  + scale_x_log10(limits=c(10^5.8,10^7.2)) + 
  geom_density_2d(color="black") + stat_smooth(method="lm",color="blue")  + ylim(c(0.2,0.75))+ theme_bw()+
  theme(text = element_text(size=20)) + xlab("Genome Length")
p2 <- ggplot(silva_gcku_full[silva_gcku_full$Ku==0,],aes(x=GenomeLength,y=GC)) + 
  geom_point(color=rgb(1,0,0,0.25)) + scale_x_log10(limits=c(10^5.8,10^7.2)) + 
  geom_density_2d(color="black") + stat_smooth(method="lm",color="red") + ylim(c(0.2,0.75))+ theme_bw()+
  theme(text = element_text(size=20)) + xlab("Genome Length")

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("GCKU_GenomeLength_spp_linear.pdf",width=8,height=5)
ggarrange(p1,p2,ncol=2,labels=c("(Ku)","(No Ku)"),hjust=-1)
dev.off()


tax_df <- tax_df %>% subset(select=c(Species, Phylum, Family)) %>% unique() %>% subset (!(is.na(Species) | is.na(Phylum)))
silva_gcku_full$Species <- silva_gcku_full$SpeciesSILVA
sgckufp <- merge.easy(silva_gcku_full,tax_df,key="Species")

sgcku_proteo <- sgckufp %>% subset(Phylum=="Proteobacteria")
sgcku_actino <- sgckufp %>% subset(Phylum=="Actinobacteria")
sgcku_firmi <- sgckufp %>% subset(Phylum=="Firmicutes")

p_actino<- ggplot(sgcku_actino,aes(x=GenomeLength,y=GC,color=factor(Ku))) + 
  geom_point(aes(fill=factor(Ku)),shape=21,alpha=0.25,size=1)  + scale_x_log10() +
  geom_density_2d(lty=2,alpha=0.7) + stat_smooth(method="lm")  + theme_bw()+
  guides(color=FALSE,fill=FALSE) + ylab("GC")+ 
  theme(text = element_text(size=15)) + xlab("Genome Length")
p_firmi<- ggplot(sgcku_firmi,aes(x=GenomeLength,y=GC,color=factor(Ku))) + 
  geom_point(aes(fill=factor(Ku)),shape=21,alpha=0.25,size=1)  + scale_x_log10() +
  geom_density_2d(lty=2,alpha=0.7) + stat_smooth(method="lm")  + theme_bw()+
  guides(color=FALSE,fill=FALSE) + ylab("GC")+ 
  theme(text = element_text(size=15)) + xlab("Genome Length")
p_proteo <- ggplot(sgcku_proteo,aes(x=GenomeLength,y=GC,color=as.factor(Ku))) + 
  geom_point(aes(fill=factor(Ku)),shape=21,alpha=0.25,size=1)  + scale_x_log10() +
  geom_density_2d(lty=2,alpha=0.7) + stat_smooth(method="lm") + ylab("GC")+ theme_bw()+
  theme(text = element_text(size=15)) + xlab("Genome Length")+ 
  scale_color_discrete(name = "Ku",labels=c("FALSE","TRUE"))+ guides(fill=FALSE)
setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("GCKU_GenomeLength_spp_Phyla.pdf",width=12,height=5)
ggarrange(p_actino,p_firmi,p_proteo,ncol=3,widths=c(3,3,4),labels=c("Actinobacteria","Firmicutes","Proteobacteria"))
dev.off()
