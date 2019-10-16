# JLW - 2019
# PCA of select ProTraits (Brbic et al. 2016, NAR) traits w/ GC content and Ku presence/absence overlays


library(dplyr)
library(data.table)
library(readr)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(cluster)
library(gridExtra)
library(RColorBrewer)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}

setwd("~/gcku/Shareable/Data")
load("GCKU.RData")
load("PTKU.RData")

#PCA on trait data
pt_ku <- pt_ku %>% as.data.frame(stringsAsFactors=F)
pca_traits <- prcomp(pt_ku[,which(!names(pt_ku) %in% c("Ku","GC"))],scale=T)
autoplot(pca_traits)

#Weighting
sort(pca_traits$rotation[,1])
sort(pca_traits$rotation[,2])

#Plot two PCs with Ku overlay
lx <- 0.9
ly <- 0.25
xy <- data.frame(list(PCA1=pca_traits$x[,1],
                      PCA2=pca_traits$x[,2],
                      Ku=pt_ku$Ku,
                      GC=pt_ku$GC))
features_traits <- data.frame(list(PCA1=pca_traits$rotation[,1],
                                   PCA2=pca_traits$rotation[,2],
                                   Feature=rownames(pca_traits$rotation)))
#http://rforpublichealth.blogspot.hk/2014/02/ggplot2-cheatsheet-for-visualizing.html
#placeholder plot - prints nothing at all
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
#scatterplot of x and y variables
scatter <- autoplot(pca_traits, data = pt_ku, colour = 'Ku',
                    loadings = TRUE, loadings.colour = 'black',loadings.label.colour="black",
                    loadings.label.repel=T,loadings.label.alpha=1,loading.label.lineheight=20,
                    loadings.label = TRUE, loadings.label.size = 3,alpha=0.5,
                    label.show.legend=F) + theme_bw() + 
  theme(legend.position=c(lx,ly),legend.justification=c(1,1),legend.box.background = element_rect(colour = "black")) 
#marginal density of x - plot on top
plot_top <- ggplot(xy, aes(PCA1, fill=factor(Ku))) + 
  geom_density(alpha=.75) + theme_bw() + 
  theme(legend.position = "none") + xlab("") 
#marginal density of y - plot on the right
plot_right <- ggplot(xy, aes(PCA2, fill=factor(Ku))) + 
  geom_density(alpha=.75) + 
  coord_flip() + theme_bw() + 
  theme(legend.position = "southeast",axis.text.x = element_text(angle = -90, hjust = 1))+ 
  xlab("")

setwd("~/gcku/Shareable/Figs")
pdf("KU_PCA.pdf",width=8,height=8)
grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


#Plot PCA w/ GC overlay
scatter <- autoplot(pca_traits, data = pt_ku, colour = 'GC',
                    loadings = F, loadings.colour = 'black',loadings.label.colour="black",
                    loadings.label.repel=T,loadings.label.alpha=1,loading.label.lineheight=20,
                    loadings.label = F, loadings.label.size = 3,alpha=0.5,
                    label.show.legend=F) + theme_bw() + scale_color_distiller(palette = "RdPu",direction=1) + 
  theme(legend.position=c(lx,ly),legend.justification=c(1,1),legend.box.background = element_rect(colour = "black")) 
scatter

#Linear model of GC content w/ respect to traits
lm_gc <- lm(GC~.,data=pt_ku)
summary(lm_gc)
setwd("~/gcku/Shareable/Data")
write.csv(summary(lm_gc)$coefficients,file="gc_lm.csv")

setwd("~/gcku/Shareable/Figs")
pdf("GCA_PCA.pdf",width=8,height=8)
scatter
dev.off()


