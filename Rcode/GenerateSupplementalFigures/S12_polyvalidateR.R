
library(data.table)
library(ggplot2)
library(dplyr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/gcku/Shareable/Data")
long <- read.csv("long2018.csv")
load("PolyGC4BiasRate_GC4expected.RData")
load("GCKU.RData")
gc_ku <- gc_ku %>% subset(select=c(Accession2,Ku))

poly_df <- poly_df %>% subset(Divergence<=0.01) %>% subset(nInformative>4)
poly_df$GC4poly <- 1/(1+poly_df$m)

poly_df$Seq <- unlist(lapply(strsplit(poly_df$Seq,split="[|]"),"[",6))
poly_df$Accession2 <- paste("GCF",gsub("^.*GCF_","",poly_df$Seq),sep="_")
poly_df <- merge.easy(poly_df,gc_ku,key="Accession2") %>% subset(!is.na(Ku))

poly_df$Organism <- gsub("_"," ",gsub("\\..*","",poly_df$Seq))
long$Organism <- as.character(long$Organism)
long$Expected <- 1/(1+long$m)

long$Species <- paste(unlist(lapply(strsplit(long$Organism,split=" "),"[",1)),
      unlist(lapply(strsplit(long$Organism,split=" "),"[",2)))
long$Species <- gsub("Burkhoderia","Burkholderia",long$Species)

pind <- lapply(long$Species,grep,poly_df$Organism)
pind_genus <- lapply(long$Genus,grep,poly_df$Organism)

long_validate <- matrix(nrow=length(unlist(pind))+length(unlist(pind_genus)),ncol=4)
counter <- 1
counter_genus <- 1
for(i in 1:length(pind)){
  if(length(pind[[i]])>0){
    for(j in 1:length(pind[[i]])){
      long_validate[counter_genus+counter,1] <- long$Expected[i]
      long_validate[counter_genus+counter,2] <- poly_df$GC4poly[unlist(pind)[counter]]
      long_validate[counter_genus+counter,3] <- i
      counter <- counter+1
    } 
  } else if(length(pind_genus[[i]])>0){
    for(j in 1:length(pind_genus[[i]])){
      long_validate[counter_genus+counter,1] <- long$Expected[i]
      long_validate[counter_genus+counter,2] <- poly_df$GC4poly[unlist(pind_genus)[counter_genus]]
      long_validate[counter_genus+counter,3] <- i
      long_validate[counter_genus+counter,4] <- 1
      counter_genus <- counter_genus+1
    }
  }
}
long_validate <- long_validate[!is.na(long_validate[,1]),]

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("ExpectedGC_Long2018vsATGC.pdf",width=7,height=5)
ggplot(data=as.data.frame(long_validate),aes(x=long_validate[,1],y=long_validate[,2],group=long_validate[,1])) + 
  geom_boxplot() + xlim(0,1) + ylim(0,1) + xlab("Expected GC (Long 2018 Estimate)") + 
  ylab("Expected GC (ATGC GC4 Estimate)") + geom_abline(slope=1,intercept=0,lty=2) 
dev.off()

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("ExpectedGC_Long2018vsATGC_JustSpp.pdf",width=7,height=5)
ggplot(data=as.data.frame(long_validate[is.na(long_validate[,4]),]),
       aes(x=long_validate[is.na(long_validate[,4]),1],
           y=long_validate[is.na(long_validate[,4]),2],
           group=long_validate[is.na(long_validate[,4]),1])) + 
  geom_boxplot() + xlim(0,1) + ylim(0,1) + xlab("Expected GC (Long 2018 Estimate)") + 
  ylab("Expected GC (ATGC GC4 Estimate)") + geom_abline(slope=1,intercept=0,lty=2) 
dev.off()







