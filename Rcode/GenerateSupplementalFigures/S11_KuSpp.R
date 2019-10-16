

library(dplyr)
library(ggplot2)


setwd("~/gcku/Shareable/Data")
load("GCKU.RData")

gcku_spp <- gc_ku %>% group_by(Species) %>% summarize(Ku=mean(Ku),N=length(Accession2))

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="Ku_Conservation.pdf",width=5,height=5)
ggplot(gcku_spp,aes(x=N,y=Ku)) + geom_point(alpha=0.025,size=5) + scale_x_log10() +
  theme_bw() + xlab("Number of Genomes for Each Species in Dataset") + 
  ylab("Percent of Genomes with Ku in a Species") + theme(text = element_text(size=14)) 
dev.off()
