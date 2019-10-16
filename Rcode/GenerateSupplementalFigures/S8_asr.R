
library(dplyr)
library(data.table)
library(corHMM)
library(ape)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/gcku/Shareable/Data")
load("GCKU.RData")
load("taxonomy.RData")
LTP <- read.tree("LTPs128_SSU_tree.newick") 
load("SILVAGCKu.RData")

tax_df <- tax_df %>% subset(select=c(Species, Genus, Family, Order, Class, Phylum)) %>% 
  unique() %>% subset(!(is.na(Species) | is.na(Genus)))
gc_ku <- merge.easy(gc_ku,tax_df,key="Species")
lgcku <- gc_ku %>% subset(Ku==1 & GC < 0.4)

table(lgcku$Genus)
table(lgcku$Family)
table(lgcku$Order)
table(lgcku$Class)
table(lgcku$Phylum)

sum(lgcku$Family=="Bacillaceae",na.rm=T)/nrow(lgcku)
sum(lgcku$Family=="Flavobacteriaceae",na.rm=T)/nrow(lgcku)
sum(lgcku$Family=="Sphingobacteriaceae",na.rm=T)/nrow(lgcku)

bspp <- gc_ku %>% subset(Family=="Bacillaceae")
bspp$GC %>% hist()

sum(bspp$GC<0.5)/nrow(bspp)
sum(bspp$GC<0.4)/nrow(bspp)

################################
# Look at Tree for Bacillaceae #
################################

#Line up tree w/ Ku data
silva_gcku_full <- silva_gcku_full %>% subset(!is.na(Ku)) %>% subset(Accession2 %in% bspp$Accession2)
species <- silva_gcku_full$NodeNames
pruned.LTP <- drop.tip(LTP,LTP$tip.label[!(LTP$tip.label %in% species)])
silva_gcku_full <- silva_gcku_full %>% subset(NodeNames %in% pruned.LTP$tip.label)
silva_gcku_full <- silva_gcku_full[!duplicated(silva_gcku_full$NodeNames),]
rownames(silva_gcku_full) <- silva_gcku_full$NodeNames
silva_gcku_full$Ku <- as.character(silva_gcku_full$Ku)

#Plot tree
setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="Ku_SILVA_TREE_Bacillaceae.pdf",width=5,height=5)
ku <- as.numeric(as.character(silva_gcku_full$Ku))
names(ku) <- rownames(silva_gcku_full)
ku <- ku[pruned.LTP$tip.label]
plot(pruned.LTP,show.tip.label=FALSE,type="fan")
tiplabels(pch = 21, cex=1,col=rgb(1-ku,0,ku,0.25),bg=rgb(1-ku,0,ku,0.25))
dev.off()

##################################
# Ancestral State Reconstruction #
##################################

ku_incidence <- silva_gcku_full %>% subset(select=c(NodeNames,Ku))

# 
# # This takes a long time to run
# one.cat <- corHMM(pruned.LTP,data=ku_incidence,rate.cat=1)
# two.cat <- corHMM(pruned.LTP,data=ku_incidence,rate.cat=2)
# setwd("~/gcku/Shareable/Data")
# save(one.cat,two.cat,file="ASR_Bacillaceae.RData")

setwd("~/gcku/Shareable/Data")
load("ASR_Bacillaceae.RData")

one.cat$AICc # 257.3119
two.cat$AICc #263.8347

setwd("~/gcku/Shareable/Figs/Supplemental")
pdf(file="GC_SILVA_ASRonerate_Bacillaceae.pdf",width=5,height=5)
plotRECON(one.cat$phy, one.cat$states,type="fan",show.tip.label=FALSE)
dev.off()
one.cat$loglik

















