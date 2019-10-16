# JLW - 2019
# Put together diverse sources of data about our set of genomes from refseq (Ku presence, GC content, etc.)
# The generated datafiles are available in the /Data folder (shouldn't need to re-run this code)

### this takes a bit of time to run - all relevant data already available in /Data folder for downstream analysis

library(dplyr)
library(data.table)
library(readr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/gcku/Shareable/Data")
load("KuSearch_full.RData") # output from hmmsearch w/ pfam profile PF02735
load("genomic_gc.RData") # output from bbmap stats.sh program 
load("PT_CRISPR_FT.RData") # transformed ProTraits (Brbic et al. 2016, NAR) data from Weissman et al. 2019, ISME J


#Bonferroni correction for # genomes searched
e_cutoff <- 1e-2/(length(unique(ku_fdf$Accession2))) 
ku_fdf <- ku_fdf %>% subset(as.numeric(as.character(full.evalue))<e_cutoff)

#Fix up accession numbers, make binary Ku column
ku_fdf$Ku <- 1
ku_fdf <- ku_fdf %>% subset(select=c(Accession2,Ku))
ku_fdf$Accession2 <- paste("GCF_",unlist(lapply(strsplit(ku_fdf$Accession2,split="_"),"[",2)),sep="")
names(gc_df)[1] <- "Accession2"
gc_ku <- merge.easy(gc_df,ku_fdf,key="Accession2")
gc_ku$Ku[is.na(gc_ku$Ku)] <- 0

# #Query NCBI for Species and Organism names associated w/ each genome
# #Uncomment below to run (takes some time)
#
# library(rentrez) 
# accessions <- unique(gc_ku$Accession2)
# getSpecies <- function(accession){
#   id <- entrez_search(db="assembly",term=accession)
#   assembly <- entrez_summary(db="assembly",id=id$ids)
#   Sys.sleep(0.01)
#   return(c(accession,assembly$speciesname,assembly$organism))
# }
# acc_spp_list <- list()
# for(i in 1:length(accessions)){
#   print(i)
#   acc_spp_list[[i]] <- try(getSpecies(accessions[i]))
# }
# acc_spp_list[[1000]]
# acc_spp_df <- as.data.frame(do.call("rbind",acc_spp_list),stringsAsFactors=FALSE)
# setwd("~/gcku/Shareable/Data")
# save(acc_spp_df,file="gcku_acc2spp.RData")
setwd("~/gcku/Shareable/Data")
load("gcku_acc2spp.RData")
names(acc_spp_df) <- c("Accession2","Species","SppName")

#Save GC vs Ku Dataset
gc_ku <- merge.easy(gc_ku,acc_spp_df,key="Accession2")
setwd("~/gcku/Shareable/Data")
save(gc_ku,file="GCKU.RData")

# Cleanup trait data, select columns for relevant traits, merge with GC/Ku data
pt_crispr$Accession2 <- rownames(pt_crispr)
names(pt_crispr) <- gsub("Score.","",names(pt_crispr))
keywords <- c("=soil",
              "oxygen",
              "spore",
              "sporu",
              "nitrogen",
              "temperature",
              "terrestrial",
              "accession")
trait_ind <- grep(paste(keywords,collapse="|"),names(pt_crispr),ignore.case = T)
pt_crispr <- pt_crispr[,trait_ind]
gc_ku <- gc_ku %>% subset(select=c(Accession2,Ku,GC))
gc_ku$Ku <- factor(gc_ku$Ku>0)
pt_ku <- merge.easy(pt_crispr,gc_ku,key="Accession2")
pt_ku$Accession2 <- NULL

#Save Trait Dataset
setwd("~/gcku/Shareable/Data")
save(pt_ku,file="PTKU.RData")

#Load in SILVA living tree Species Names and Accessions
setwd("~/gcku/Shareable/Data")
load("SILVA_LivingTree_accession.RData")
names(SILVA_LivingTree_accession_df) <- c("AccessionSILVA","SpeciesSILVA","Accession2")
SILVA <- SILVA_LivingTree_accession_df %>% unique()

#Load SILVA Tree - Match nodenames to accessions
lines <- readLines("LTPs128_SSU_tree.newick")
lines <- unique(unlist(strsplit(lines,split="[(]")))
lines <- unique(unlist(strsplit(lines,split=")")))
lines <- unique(unlist(strsplit(lines,split=":")))
SILVA$NodeNames <- lapply(SILVA$AccessionSILVA,grep,x=lines,value=TRUE)
SILVA <- as.data.frame(lapply(SILVA,as.character),stringsAsFactors = FALSE)
#Some naming matching issues
SILVA$NodeNames[SILVA$AccessionSILVA=="Y18835"] <- "Hymenobacter_ocellatus__Y18835__Cytophagaceae__HymOcell"
SILVA$NodeNames[SILVA$AccessionSILVA=="X55598"] <- "Mycobacterium_aichiense__X55598__Mycobacteriaceae__MyoAichi"
SILVA$NodeNames[SILVA$AccessionSILVA=="Y18654"] <- "Lactobacillus_fornicalis__Y18654__Lactobacillaceae__LtbForni"
SILVA$NodeNames[SILVA$AccessionSILVA=="Y18097"] <- "Vagococcus_salmoninarum__Y18097__Enterococcaceae__VagSalmo"
SILVA$NodeNames[SILVA$AccessionSILVA=="M24483"] <- "Spiroplasma_poulsonii__M24483__Spiroplasmataceae__SrsPouls"
SILVA$NodeNames[SILVA$AccessionSILVA=="M59083"] <- "Acetitomaculum_ruminis__M59083__Lachnospiraceae__AcmRumin"
SILVA$NodeNames[SILVA$AccessionSILVA=="CP001336"] <- "Desulfitobacterium_hafniense__CP001336__Peptococcaceae__DscHafn2"
SILVA$NodeNames[SILVA$AccessionSILVA=="U88891"] <- "Desulfotomaculum_halophilum__U88891__Peptococcaceae__DtmHalop"
SILVA$NodeNames[SILVA$AccessionSILVA=="CP000924"] <- "Thermoanaerobacter_pseudethanolicus__CP000924__Thermoanaerobacteraceae__TmrPseu2"
SILVA$NodeNames[SILVA$AccessionSILVA=="AY596297"] <- "Haloarcula_marismortui__AY596297__Halobacteriaceae__HrcMaris"
SILVA <- rbind(SILVA,c(0,0,0,SILVA$Classification[SILVA$AccessionSILVA=="CP001336"],"Desulfitobacterium_hafniense__CP001336__Peptococcaceae__DscHafni",SILVA$CRISPR[SILVA$AccessionSILVA=="CP001336"]))
SILVA <- rbind(SILVA,c(0,0,0,SILVA$Classification[SILVA$AccessionSILVA=="CP000924"],"Thermoanaerobacter_pseudethanolicus__CP000924__Thermoanaerobacteraceae__TmrPseud",SILVA$CRISPR[SILVA$AccessionSILVA=="CP000924"]))
SILVA <- rbind(SILVA,c(0,0,0,SILVA$Classification[SILVA$AccessionSILVA=="AY596297"],"Haloarcula_marismortui__AY596297__Halobacteriaceae__HrcMari2",SILVA$CRISPR[SILVA$AccessionSILVA=="AY596297"]))

SILVA <- as.data.frame(lapply(SILVA,as.character),stringsAsFactors = FALSE)
setwd("~/gcku/Shareable/Data")
save(SILVA,file="SILVA1.RData")

#Oxygen data from NCBI
load("biosample_data_NCBI_processed.RData")
oxygen_data <- sample_data %>% subset(!is.na(`Oxygen Requirement`)|!is.na(rel_to_oxygen)|!is.na(`Relationship to oxygen`)) %>% 
  subset(select=c(Accession2,`Oxygen Requirement`,rel_to_oxygen,`Relationship to oxygen`))
oxygen_data$OxygenRequirement <- oxygen_data$rel_to_oxygen
oxygen_data$OxygenRequirement[is.na(oxygen_data$OxygenRequirement)] <- oxygen_data$`Oxygen Requirement`[is.na(oxygen_data$OxygenRequirement)]
oxygen_data$OxygenRequirement[is.na(oxygen_data$OxygenRequirement)] <- oxygen_data$`Relationship to oxygen`[is.na(oxygen_data$OxygenRequirement)]
oxygen_data <- oxygen_data %>% subset(select=c(Accession2,OxygenRequirement))
oxygen_data$OxygenRequirement <- gsub("bic","be",tolower(oxygen_data$OxygenRequirement))
oxygen_data$OxygenRequirement <- gsub("facilitative","facultative",oxygen_data$OxygenRequirement)
oxygen_data$OxygenRequirement <- gsub("facultatively","facultative",oxygen_data$OxygenRequirement)
oxygen_data$OxygenRequirement <- gsub("obligate ","",oxygen_data$OxygenRequirement)
oxygen_data <- oxygen_data %>% subset(!(OxygenRequirement %in% c("-","unknown","missing")))
unique(oxygen_data$OxygenRequirement)
table(oxygen_data$OxygenRequirement)

#Put SILVA tree data together w/ GC and Ku
gc_ku <- merge.easy(gc_ku,oxygen_data,key="Accession2")
gc_ku <- gc_ku[!is.na(gc_ku$GC),] 
gc_ku[is.na(gc_ku)] <- "NA"
silva_gcku_full <- merge.easy(SILVA,gc_ku,key="Accession2")
silva_gcku <- silva_gcku_full %>% subset(OxygenRequirement %in% c("aerobe","anaerobe"))
setwd("/home/jlweissman/gcku/Data")
save(silva_gcku,silva_gcku_full,file="SILVAGCKu.RData")

