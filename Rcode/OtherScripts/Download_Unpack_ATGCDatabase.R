# JLW - 2019
# Download ATGC clusters w/ at least three members, unpack into format useable for polymorphism analysis 
# Note that because the unpacked db is >12gb it is not included in this repository

#Make folder to hold db
system("mkdir ~/gcku/Shareable/ATGC; mkdir ~/gcku/Shareable/ATGC/ATGC_nucleotidealignments ; mkdir ~/gcku/Shareable/ATGC/ATGC_unpacked")

#Download clusters w/ >=3 members
setwd("~/gcku/Shareable/Data")
atgc_metadata <- read.csv("atgcdata.csv") #From ATGC website - Number of members in each cluster
to_download <- as.character(atgc_metadata$ATGC.[atgc_metadata$numgenomes>=3])
download_folder <- "~/gcku/Shareable/ATGC/ATGC_nucleotidealignments"
download_links <- paste("http://dmk-brain.ecn.uiowa.edu/ATGC/data/",to_download,
                        "/",to_download,".indexnuclalns.zip",sep="")
commands <- paste("cd ",download_folder,"; wget -q ",download_links,sep="")
lapply(commands,system)
#Unzip
system("cd ~/gcku/Shareable/ATGC/ATGC_nucleotidealignments; unzip '*.zip'")

#################################################################################################################
# Unpack Clusters (concatenate core gene alignments for each genome in a cluster, make new fasta for each genome)

grabOrganism <- function(organism,cog_keep_ind,header_inds,headers,atgc_file,atgc_id){
  org_header_ind <- intersect(grep(organism,headers),cog_keep_ind)
  lines <- paste(paste(header_inds[org_header_ind]+1,header_inds[org_header_ind+1]-1,sep=","),collapse="p -ne ")
  lines <- paste(header_inds[org_header_ind[1]],lines,sep="p -ne ")
  filename <- paste("~/gcku/Shareable/ATGC/ATGC_unpacked/ATGC",
                    atgc_id,"/",organism,".fasta",sep="")
  # print(paste("cd ~/ATGC_test; sed -ne ",lines,"p ",atgc_file," > ",filename,sep=""))
  system(paste("cd ~/gcku/Shareable/ATGC/ATGC_nucleotidealignments; sed -ne ",
               lines,"p ",atgc_file," > ",filename,sep=""))
}

unpackATGC <- function(atgc_id){
  atgc_file <- paste("ATGC",atgc_id,".indexnuclalns.fa",sep="")
  foldername <- paste("~/gcku/Shareable/ATGC/ATGC_unpacked/ATGC",atgc_id,"/",sep="")
  system(paste("mkdir -p",foldername))
  system(paste('cd ~/gcku/Shareable/ATGC/ATGC_nucleotidealignments; grep "taxonId" ',atgc_file,
               ' > ',foldername,'headers.txt',sep=""))
  
  setwd(foldername)
  x <- readLines("headers.txt")
  organisms <- unique(unlist(lapply(strsplit(x,split="[|]"),"[",6)))
  cogs <- (unlist(lapply(strsplit(x,split="[|]"),"[",2)))
  # hist(table(cogs))
  cogs_core <- names(table(cogs))[table(cogs)==length(table(organisms))]
  # cog_keep_ind <- grep(paste(cogs_core,collapse="|"),x)
  cog_keep_ind <- unlist(lapply(cogs_core,grep,x=x))
  system(paste('cd ~/gcku/Shareable/ATGC/ATGC_nucleotidealignments/; grep -Fn "CogId" ',
               atgc_file,' | cut -d : -f 1 > ',foldername,'header_indices.txt',sep=""))
  header_inds <- as.numeric(readLines("header_indices.txt"))
  
  system(paste("mkdir -p",foldername))
  lapply(organisms,grabOrganism,cog_keep_ind=cog_keep_ind,header_inds=header_inds,headers=x,atgc_file=atgc_file,atgc_id=atgc_id)
  
  system(paste("wc -l ",foldername,"*.fasta > ",foldername,"lengths.txt",sep=""))
  l <- readLines(paste(foldername,"lengths.txt",sep=""))
  lengths <- unlist(as.numeric(lapply(strsplit(trimws(l),split=" "),"[",1)))
  lengths <- lengths[1:(length(lengths)-1)]
  files <- unlist(lapply(strsplit(trimws(l),split=" "),"[",2))
  files <- paste(files[which(lengths==max(lengths))],collapse=" ")
  
  system(paste("cd ",foldername,"; cat ",files," > ATGC",atgc_id,".fasta",sep=""))
  system(paste("cd ",foldername,"; snp-sites -o ATGC",atgc_id,"snp.fasta ATGC",atgc_id,".fasta",sep=""))
}

setwd("~/gcku/Shareable/Data/")
atgc_metadata <- read.csv("atgcdata.csv")
to_unpack <- gsub("ATGC","",as.character(atgc_metadata$ATGC.[atgc_metadata$numgenomes>=3 & atgc_metadata$numgenomes<5]))
lapply(to_unpack,unpackATGC)
