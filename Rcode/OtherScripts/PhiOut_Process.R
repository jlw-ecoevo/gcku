# JLW - 2019

# Compiles p-values from PhiPack and does BH correction

# First run:
# cd ~/gcku/Shareable/ATGC/PhiOutPerm/
# grep -R "PHI (Permutation):" > ../phi_p_vals.txt


library(dplyr)

setwd("~/gcku/Data")
phi_p <- readLines("phi_p_vals.txt")
ATGC <- phi_p %>% gsub(pattern="[.].*",replace="")
COG <- phi_p %>% strsplit(split="[.]") %>% lapply("[",2) %>% unlist()
Pval <- phi_p %>% gsub(pattern=".*:",replace="") %>% 
  gsub(pattern="[()].*",replace="") %>% trimws() %>% as.numeric()


hist(Pval%>%log10())

phi_df <- data.frame(ATGC=ATGC,
                     COG=COG,
                     Pval=Pval,
                     PvalCorrectedBH=p.adjust(Pval,method = "BH"))


setwd("~/gcku/Shareable/Data")
save(phi_df,file="PhiRecombTest.RData")
