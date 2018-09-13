library(Biostrings)
library(pbapply)
library(msa)
library(dplyr)

source("./functions/read_uniprot.R")

seq2bio <- function(x)
  AAString(paste0(x, collapse = ""))

remove_uncultured <- function(x) 
  ifelse(x[[1]] == "uncultured", x[[2]], x[[1]])


unify_gammas <- function(x)
  ifelse(x %in% c("gamma", "Gammaproteobacteria", "Gamma-proteobacterium"), 
         "Gamma-proteobacterium", x[[1]])

csga <- read_uniprot("./data/uniprot-csgA.txt", "signal")
csgc <- read_uniprot("./data/uniprot-csgC.txt", "signal")

ecoli <- csga[unname(sapply(csga, function(i) grepl("Escherichia", attr(i, "OS")) & grepl("coli", attr(i, "OS")) ))]

#table(lengths(ecoli))
#29 102 110 149 151 152 
#1   1   4   1  97  26
# only 149 to 152 are csgAs
ecoli <- ecoli[lengths(ecoli) < 200 & lengths(ecoli) > 148]

# PubMed:1357528 and PubMed:1677357 - signal ppetide 20 aa long

library(seqinr)
write.fasta(sapply(ecoli, function(i) i[attr(i, "SP"):length(i)]), 
            names = names(ecoli), file.out = "ecoli_csga.fasta")

#salmonella <- csga[unname(sapply(csga, function(i) grepl("Salmonell", attr(i, "OS"))))]

