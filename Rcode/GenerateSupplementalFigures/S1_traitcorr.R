# JLW - 2019

library(corrplot)

setwd("~/gcku/Shareable/Data")
load("PTKU.RData")

pt_ku$Ku <- as.numeric(as.logical(as.character(pt_ku$Ku)))
M <- cor(pt_ku[,1:15])
corrplot(M,method="pie")

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="CorrTrait.pdf",width=6,height=6)
corrplot(M,method="pie")
dev.off()

