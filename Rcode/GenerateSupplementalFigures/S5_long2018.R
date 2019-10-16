
library(lattice)
library(ggplot2)

setwd("~/gcku/Shareable/Data")
long <- read.csv("~/gcku/Data/long2018.csv")
long$expectGC <- 1/(1+long$m)

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="MutationBias_ExpectedGC_Long2018.pdf",width=5,height=5)
vioplot(long$expectGC[long$Ku==0],long$expectGC[long$Ku==1],col="gray",names=c("No Ku","Ku"))
title(ylab="Expected GC Content")
stripchart(expectGC~Ku, vertical = TRUE, data = long, 
           method = "jitter", add = TRUE, pch = 21, col = 'black',bg='blue',cex=1)
dev.off()

