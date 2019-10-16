# JLW - 2019

#Will run PhiPack on a set of gene alignments (stepsize argument controls size for each job)
# Need phipack to be installed at ~/gcku/Shareable/ATGC/PhiPack

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (for file number)", call.=FALSE)
} else {
  file_ind <- as.numeric(as.character(args[1]))
  stepsize <- as.numeric(as.character(args[2]))
}

files <- readLines("~/gcku/Shareable/ATGC/gene_alignment_files.txt")
file_ind_todo <- ((file_ind-1)*stepsize + 1):min((file_ind*stepsize),length(files))
file <- files[file_ind_todo]

phi_cmd <- paste0("cd ~/gcku/Shareable/ATGC/PhiPack; ./Phi -p 10000 -v TRUE -f ",
                  "~/gcku/Shareable/ATGC/GeneAlignments/",file,
                  " > ",
                  "~/gcku/Shareable/ATGC/PhiOutPerm/",gsub(".fasta",".out",file))

lapply(phi_cmd,system)