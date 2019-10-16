# JLW - 2019

library(dplyr)

setwd("~/gcku/Shareable/Data")
load("PTKU.RData")

pt_ku$Ku <- pt_ku$Ku %>% as.logical() %>% as.numeric()
x <- cor(pt_ku)
y <- x[1:15,"Ku"]
x <- x[1:15,"GC"]

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf("CorGCKu.pdf",width=7,height=5)
par(mar=c(5,5,2,2))
plot(x,y,xlab="Correlation with GC Content",
     ylab="Correlation with Ku Incidence",xlim=c(-.5,0.6),ylim=c(-.5,0.6),
     cex.lab=1.5,cex.axis=1.5)
abline(a=0,b=1,lty=2)
dev.off()
