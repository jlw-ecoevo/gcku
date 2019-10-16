
library(dplyr)
library(data.table)
library(readr)
library(ape)
library(phylolm)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}


setwd("~/gcku/Shareable/Data")
load("SILVAGCKu.RData")
load("refseq_all_genome_lengths.RData")
LTP <- read.tree("LTPs128_SSU_tree.newick")
load("taxonomy.RData")

#Add genome length column
silva_gcku_full <- merge.easy(silva_gcku_full,gl_df,key="Accession2")
silva_gcku_full <- silva_gcku_full %>% subset(!is.na(GC) & !is.na(GenomeLength))
silva_gcku_full$GenomeLength <- as.numeric(silva_gcku_full$GenomeLength)


#Ok to sample one genome per tip on SILVA tree  - Ku incidence and GC content basically same across genomes in spp
x <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% summarize(Ku=mean(as.numeric(Ku)))
hist(x$Ku)
sum(x$Ku <1 & x$Ku>0)/nrow(x)
x <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% summarize(VGC=var(GC))
max(x$VGC,na.rm=T)
hist(x$VGC,breaks=100)
#Do sampling
silva_gcku_full <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% sample_n(1)


####################
# Phylo Correction #
####################

#Logit transform GC column
silva_gcku_full$GCorig <- silva_gcku_full$GC
silva_gcku_full$GC <- log(silva_gcku_full$GC/(1-silva_gcku_full$GC))

#Uncorrected model
lm_gc <- lm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),data=silva_gcku_full)
summary(lm_gc)

#Prune SILVA tree - -only nodes for which we have genomes
#Make sure data rows align w/ node names
species <- silva_gcku_full$NodeNames
pruned.LTP <- drop.tip(LTP,LTP$tip.label[!(LTP$tip.label %in% species)])
sum(pruned.LTP$edge.length==0)
min(pruned.LTP$edge.length[pruned.LTP$edge.length>0])
pruned.LTP$edge.length <- pruned.LTP$edge.length + 1e-8
silva_gcku_full <- silva_gcku_full %>% subset(NodeNames %in% pruned.LTP$tip.label)
silva_gcku_full <- silva_gcku_full[!duplicated(silva_gcku_full$NodeNames),]
rownames(silva_gcku_full) <- silva_gcku_full$NodeNames
silva_gcku_full$Ku <- as.character(silva_gcku_full$Ku)

#Number of Ku/No-Ku taxa
table(ku=silva_gcku_full$Ku)

#Orenstein-Uhlenbeck Model (ALL)
phy_mod <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                   data=silva_gcku_full,phy=pruned.LTP,model="OUrandomRoot",
                   measurement_error=T)
summary(phy_mod)
phy_mod$aic

#Brownian-Motion Model (ALL)
phy_mod <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                   data=silva_gcku_full,phy=pruned.LTP,model="BM",
                   measurement_error=T)
summary(phy_mod)
phy_mod$aic

################
# Phylum Model #
################

#Add column for phyla
tax_df <- tax_df %>% subset(select=c(Species, Phylum, Family)) %>% unique() %>% subset (!(is.na(Species) | is.na(Phylum)))
silva_gcku_full$Species <- silva_gcku_full$SpeciesSILVA
sgckufp <- merge.easy(silva_gcku_full,tax_df,key="Species")

#Phyla
table(sgckufp$Phylum)[which(table(sgckufp$Phylum)>1000)] # 3 phyla w/ 1000+ spp
sgcku_proteo <- sgckufp %>% subset(Phylum=="Proteobacteria")
sgcku_actino <- sgckufp %>% subset(Phylum=="Actinobacteria")
sgcku_firmi <- sgckufp %>% subset(Phylum=="Firmicutes")

#Number of Ku/No-Ku taxa
table(sgcku_actino$Ku)
table(sgcku_proteo$Ku)
table(sgcku_firmi$Ku)


### Proteobacteria

rownames(sgcku_proteo) <- sgcku_proteo$NodeNames
proteo_tree <- drop.tip(pruned.LTP,
                        pruned.LTP$tip.label[!(pruned.LTP$tip.label %in% sgcku_proteo$NodeNames)])
phy_mod_proteo <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                   data=sgcku_proteo,phy=proteo_tree,model="BM",
                   measurement_error=T)
summary(phy_mod_proteo)
phy_mod_proteo$aic
#
phy_mod_proteo <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                          data=sgcku_proteo,phy=proteo_tree,model="OUrandomRoot",
                          measurement_error=T)
summary(phy_mod_proteo)
phy_mod_proteo$aic


### Actinobacteria

rownames(sgcku_actino) <- sgcku_actino$NodeNames
actino_tree <- drop.tip(pruned.LTP,
                        pruned.LTP$tip.label[!(pruned.LTP$tip.label %in% sgcku_actino$NodeNames)])
phy_mod_actino <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                          data=sgcku_actino,phy=actino_tree,model="BM",
                          measurement_error=T)
summary(phy_mod_actino)
phy_mod_actino$aic
#
phy_mod_actino <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                          data=sgcku_actino,phy=actino_tree,model="OUrandomRoot",
                          measurement_error=T)
summary(phy_mod_actino)
phy_mod_actino$aic

### Firmicutes

rownames(sgcku_firmi) <- sgcku_firmi$NodeNames
firmi_tree <- drop.tip(pruned.LTP,
                        pruned.LTP$tip.label[!(pruned.LTP$tip.label %in% sgcku_firmi$NodeNames)])
phy_mod_firmi <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                          data=sgcku_firmi,phy=firmi_tree,model="BM",
                          measurement_error=T)
summary(phy_mod_firmi)
phy_mod_firmi$aic
#
phy_mod_firmi <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                         data=sgcku_firmi,phy=firmi_tree,model="OUrandomRoot",
                         measurement_error=T)
summary(phy_mod_firmi)
phy_mod_firmi$aic

############################################
# Models Just w/ Genera Where Ku Conserved #
############################################

setwd("~/gcku/Shareable/Data")
load("SILVAGCKu.RData")
load("refseq_all_genome_lengths.RData")
LTP <- read.tree("LTPs128_SSU_tree.newick")
load("taxonomy.RData")

#Add genome length column
silva_gcku_full <- merge.easy(silva_gcku_full,gl_df,key="Accession2")
silva_gcku_full <- silva_gcku_full %>% subset(!is.na(GC) & !is.na(GenomeLength))
silva_gcku_full$GenomeLength <- as.numeric(silva_gcku_full$GenomeLength)

#Add column for Genus
silva_gcku_full$Species <- silva_gcku_full$SpeciesSILVA
tax_df <- tax_df %>% subset(select=c(Species, Genus)) %>% 
  unique() %>% subset (!(is.na(Species) | is.na(Genus))) %>%
  subset(Species %in% silva_gcku_full$Species)
silva_gcku_full<- merge.easy(silva_gcku_full,tax_df,key="Species")

#Restrict to Genera w/ multiple representatives and uniform Ku
x <- silva_gcku_full %>% group_by(Genus) %>% count()
y <- silva_gcku_full %>% group_by(Genus) %>% summarize(Kuf=mean(Ku))
silva_gcku_full <- silva_gcku_full %>% 
  subset(Genus %in% x$Genus[x$n>1]) %>%
  subset(Genus %in% y$Genus[y$Kuf > 0.975 | y$Kuf < 0.025])

#Ok to sample one genome per tip on SILVA tree  - Ku incidence and GC content basically same across genomes in spp
x <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% summarize(Ku=mean(as.numeric(Ku)))
hist(x$Ku)
sum(x$Ku <1 & x$Ku>0)/nrow(x)
x <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% summarize(VGC=var(GC))
max(x$VGC,na.rm=T)
hist(x$VGC,breaks=100)
#Do sampling
silva_gcku_full <- silva_gcku_full %>% group_by(SpeciesSILVA) %>% sample_n(1)

#Logit transform GC column
silva_gcku_full$GCorig <- silva_gcku_full$GC
silva_gcku_full$GC <- log(silva_gcku_full$GC/(1-silva_gcku_full$GC))

#Uncorrected model
lm_gc <- lm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),data=silva_gcku_full)
summary(lm_gc)

#Prune SILVA tree - -only nodes for which we have genomes
#Make sure data rows align w/ node names
species <- silva_gcku_full$NodeNames
pruned.LTP <- drop.tip(LTP,LTP$tip.label[!(LTP$tip.label %in% species)])
sum(pruned.LTP$edge.length==0)
min(pruned.LTP$edge.length[pruned.LTP$edge.length>0])
pruned.LTP$edge.length <- pruned.LTP$edge.length + 1e-8
silva_gcku_full <- silva_gcku_full %>% subset(NodeNames %in% pruned.LTP$tip.label)
silva_gcku_full <- silva_gcku_full[!duplicated(silva_gcku_full$NodeNames),]
rownames(silva_gcku_full) <- silva_gcku_full$NodeNames
silva_gcku_full$Ku <- as.character(silva_gcku_full$Ku)

#Number of Ku/No-Ku taxa
table(ku=silva_gcku_full$Ku)

#Orenstein-Uhlenbeck Model (ALL)
phy_mod <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                   data=silva_gcku_full,phy=pruned.LTP,model="OUrandomRoot",
                   measurement_error=T)
summary(phy_mod)
phy_mod$aic

#Brownian-Motion Model (ALL)
phy_mod <- phylolm(GC~Ku + log10(GenomeLength) + Ku:log10(GenomeLength),
                   data=silva_gcku_full,phy=pruned.LTP,model="BM",
                   measurement_error=T)
summary(phy_mod)
phy_mod$aic






