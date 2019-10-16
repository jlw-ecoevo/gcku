# JLW - 2019
# Summarizes the information about sites flanking restriction sites found in Bases_Flanking_restriction_Sites.R
# Compares to Null distribution (see Generate_Null_RestrictionEnzymes.R and Bases_Flanking_restriction_Sites.R)
# *** The product of this script ("GCRMFlanks.RData") is already available in the /Data folder - if you want to run this script yourself you will first need to run the scripts noted above

library(dplyr)
library(data.table)
library(ggplot2)
library(seqinr)
library(ggpubr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

#species names
setwd("~/gcku/Shareable/Data")
load("gcku_acc2spp.RData")
names(acc_spp_df) <- c("Accession","Species","Organism")

#Compile actual flanks
addFlank <- function(x,flankdist){
  xf <- read.delim(paste("~/gcku/Shareable/Data/FlankGCRfast_flankdist",flankdist,
                         ".txt",sep=""),header = F) %>% subset(!is.na(V4)) 
  xf$ID <- paste(xf$V1,xf$V3,sep="_")
  xf <- xf %>% subset(select=c(ID,V4))
  xf[[2]] <- as.numeric(as.character(xf[[2]]))
  names(xf)[2] <- paste("Flank",flankdist,sep="")
  return(merge.easy(x,xf,key="ID"))
}
x <- read.delim("~/gcku/Shareable/Data/FlankGCRfast_flankdist1.txt",header = F) %>% subset(!is.na(V4)) 
x$Accession <- x$V1 %>% gsub(pattern=".*/",replace="") %>% 
  strsplit(split="_") %>% lapply("[",2) %>% unlist() %>% paste0("GCF_",.)
x <- merge.easy(x,acc_spp_df,key="Accession")
x$ID <- paste(x$V1,x$V3,sep="_")
x$V4 <- as.numeric(as.character(x$V4))
names(x)[2:6] <- c("File","GCseq","RecSeq","Flank1","nSites")
for(i in 2:200){
  print(i)
  x <- addFlank(x,i)
}
x$LenRec <- nchar(x$RecSeq%>% as.character())
x$GCrec <- unlist(lapply(strsplit(as.character(x$RecSeq),split=""),GC))
xfull <- x


#Compile null flanks
addFlankN <- function(x,flankdist){
  xf <- read.delim(paste("~/gcku/Shareable/Data/FlankGCRfast_NULLAltSeq_flankdist",flankdist,
                         ".txt",sep=""),header = F) %>% subset(!is.na(V4)) 
  xf$ID <- paste(xf$V1,xf$V3,sep="_")
  xf <- xf %>% subset(select=c(ID,V4))
  xf[[2]] <- as.numeric(as.character(xf[[2]]))
  names(xf)[2] <- paste("Flank",flankdist,sep="")
  return(merge.easy(x,xf,key="ID"))
}
xN <- read.delim("~/gcku/Shareable/Data/FlankGCRfast_NULLAltSeq_flankdist1.txt",header = F) %>% subset(!is.na(V4)) 
xN$ID <- paste(xN$V1,xN$V3,sep="_")
xN$V4 <- as.numeric(as.character(xN$V4))
names(xN) <- c("File","GCseq","RecSeq","Flank1","nSites","ID")
for(i in 2:200){
  print(i)
  xN <- addFlankN(xN,i)
}
xN$LenRec <- nchar(xN$RecSeq%>% as.character())
xN$GCrec <- unlist(lapply(strsplit(as.character(xN$RecSeq),split=""),GC))
xNfull <- xN

bootMeanDist <- function(x,n=1000){
  boot_vec <- numeric(n)
  for(i in 1:n){
    boot_vec[i] <- mean(sample(x,length(x),replace=T),na.rm=T)
  }
  return(boot_vec)
}

summarizeFlank <- function(x,flankdist){
  x <- x[which(rowSums(is.na(x))==0),]
  flnk <- paste("Flank",flankdist,sep="")
  flnkN <- paste("NFlank",flankdist,sep="")
  bootflnk <- bootMeanDist(as.numeric(as.character(x[[flnk]]))-as.numeric(as.character(x[[flnkN]])))
  distdiff_df <- data.frame(Dist=flankdist,
                            GCdiff=mean(as.numeric(as.character(x[[flnk]]))-as.numeric(as.character(x[[flnkN]])),na.rm=T),
                            GCnf=mean(as.numeric(as.character(x[[flnkN]])),na.rm=T),
                            GCf=mean(as.numeric(as.character(x[[flnk]])),na.rm=T),
                            CIlower=quantile(bootflnk,0.025,na.rm=T),
                            CIupper=quantile(bootflnk,0.975,na.rm=T),
                            Qlower=quantile(as.numeric(as.character(x[[flnk]]))-as.numeric(as.character(x[[flnkN]])),0.025,na.rm=T),
                            Qupper=quantile(as.numeric(as.character(x[[flnk]]))-as.numeric(as.character(x[[flnkN]])),0.975,na.rm=T))
  return(distdiff_df)
}


x <- xfull %>% subset(GCrec<=0.25)
xN <- xNfull %>% subset(GCrec<=0.25)
xN$GCseq <- xN$GCseq %>% as.character() %>% as.numeric()
x$ID2 <- paste0(x$File,x$GCseq)
xN$ID2 <- paste0(xN$File,xN$GCseq)
names(xN)[grep("Flank",names(xN))] <- paste0("N",names(xN)[grep("Flank",names(xN))] )
# x <- x %>% group_by(ID2) %>% sample_n(1)
# xN <- xN %>% group_by(ID2) %>% sample_n(1)
xa <- merge.easy(x,xN,key="ID2")
distdiff_df <- data.frame()
for(i in 1:200){
  print(i)
  distdiff_df <- rbind(distdiff_df,summarizeFlank(xa,i))
}

setwd("~/gcku/Shareable/Data")
save(distdiff_df,file="GCRMFlanks.RData")

