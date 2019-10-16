
library(dplyr)
library(data.table)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/gcku/Shareable/Data")
load("GCKU.RData")
hrr <- read.csv("HRR_VosAndDidelot2009.csv")
hr <- read.csv("capsulerecomb.csv")
gc_ku <- gc_ku %>% subset(select=c(Ku,Species,Accession2)) #%>% unique()

hr$Organism <- paste(hr$Genus,hr$Species)
# Some spp names are synonyms:
hr$Organism[!hr$Organism %in% gc_ku$Species]
hr$Organism[hr$Organism=="Borrelia burgdorferi"] <- "Borreliella burgdorferi"
hr$Organism[hr$Organism=="Enterobacter aerogenes"] <- "Klebsiella aerogenes"
hr$Organism[hr$Organism=="Haemophilus ducreyi"] <- "[Haemophilus] ducreyi"
hr$Organism[hr$Organism=="Mycobacterium abscessus"] <- "Mycobacteroides abscessus"
hr$Organism[hr$Organism=="Propionibacterium acnes"] <- "Cutibacterium acnes"
#Calculate recomb %
hr$HRCHI <- hr$CHI/hr$PHI.families
hr$HRPHI <- hr$PHI/hr$PHI.families
#
hr <-  hr %>% subset(select=c(Organism,HRCHI,HRPHI))
names(hr)[1] <- "Species"
kuhr <- merge.easy(gc_ku,hr,key="Species")
kuhr <- kuhr %>% group_by(Species) %>% summarise_all(mean,na.rm=T)

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("Ku_vs_HR_PHI.pdf",width=7,height=5)
plot(kuhr$Ku,kuhr$HRPHI,xlab="Ku Frequency in Species",ylab="PHI",cex.axis=1.5,cex.lab=1.5)
lm_phi <- lm(HRPHI~Ku,data=kuhr)
abline(lm_phi,col="red")
dev.off()

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("Ku_vs_HR_CHI.pdf",width=7,height=5)
plot(kuhr$Ku,kuhr$HRCHI,xlab="Ku Frequency in Species",ylab="CHI",cex.axis=1.5,cex.lab=1.5)
lm_phi <- lm(HRPHI~Ku,data=kuhr)
abline(lm_phi,col="red")
dev.off()

hrr$Species <- as.character(hrr$Species)
hrr$Species[hrr$Species=="Clostridium difficile"] <- "Clostridioides difficile"
hrr$Species[2]
hrr$Species[2] <- "Candidatus Pelagibacter ubique"
hrr$Species[hrr$Species=="Haemophilus parasuis"] <- "Glaesserella parasuis"
hrr$Species[hrr$Species=="Bacillus weihenstephanensis"] <- "Bacillus mycoides"
hrr$Species[hrr$Species=="Microcoleus chthonoplastes"] <- "Coleofasciculus chthonoplastes"
hrr$Species[34]
hrr$Species[34] <- "Escherichia coli"
hrr$Species[!hrr$Species %in% gc_ku$Species]

hrrku <- merge.easy(hrr,gc_ku,key="Species")
hrrku <- hrrku %>% group_by(Species) %>% summarise_all(mean)

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("Ku_vs_HR.pdf",width=7,height=5)
plot(r.m~Ku,data=hrrku,log="y",ylab="r/m",xlab="Ku Frequency in Species",cex.axis=1.5,cex.lab=1.5)
abline(lm(log10(r.m)~Ku,data=hrrku),col="red")
dev.off()

