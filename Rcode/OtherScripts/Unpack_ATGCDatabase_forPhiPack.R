# JLW - 2019

#Unpack ATGC database into single gene alignments (one for each cluster-gene pair)
# argument input to easily start separate jobs (one per ATGC cluster) on compute cluster

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (for file number)", call.=FALSE)
} else {
  file_ind <- as.numeric(as.character(args[1]))
}

getGene <- function(file){
  #Find all sequence starts w/ descriptions
  header_cmd <- paste0('grep -n ">CogId" ',file,' > ',file,'.header')
  system(header_cmd)
  #List of ATGC/COG combos
  headers <- system(paste0('awk -F "|" \'{ print $2 }\' ',file,'.header'),intern=T) 
  unique_headers <- headers %>% unique()
  header_lines <- system(paste0('awk -F ":" \'{ print $1 }\' ',file,'.header'),intern=T)
  header_ind <- c(header_lines,
                  as.numeric(system(paste0('wc -l <',file),intern=T))+1) %>% as.numeric
  #Extract each ATGC/GOC alignment
  for(h in unique_headers){
    gene_ind <- which(headers==h)
    start_align <- header_ind[gene_ind[1]]
    end_align <- header_ind[gene_ind[length(gene_ind)]+1]-1
    write_cmd <- paste0('sed -n ',
                        format(start_align,scientific=FALSE),
                        ',',
                        format(end_align,scientific=FALSE),
                        'p ',
                        file,
                        ' > '
                        ,h,
                        '.fasta')
    system(write_cmd)
  }
}

#system("cd ~/gcku/Shareable/ATGC/ATGC_nucleotidealignments;; ls *.fa > ../files.txt")
files <- readLines("~/gcku/Shareable/ATGC/files.txt") %>% paste0("~/gcku/Shareable/ATGC/ATGC_nucleotidealignments;",.)
setwd("~/gcku/Shareable/ATGC/GeneAlignments")
getGene(files[file_ind])






