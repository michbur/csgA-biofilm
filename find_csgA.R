library(biogram)
library(seqinr)
library(dplyr)

seqs <- read_fasta("/home/michal/Dropbox/dropbox-amylogram/E. coli biofilm screen/Sequences/csg/csg_all.fasta")
seq_k12 <- read_fasta("data/CsgA-K12-gene.fasta")

find_sequence <- function(pattern, seqs, max.mismatch = 60) {
  seq2bio <- function(x)
    Biostrings::AAString(paste0(x, collapse = ""))
  
  pattern_bs <- seq2bio(pattern)
  
  lapply(seqs, function(ith_seq) {

    matching_seq <- Biostrings::matchPattern(pattern_bs, seq2bio(ith_seq), max.mismatch = max.mismatch, 
                                             min.mismatch = 0, with.indels = TRUE, algorithm = "auto")
    
    ranges <- as.matrix(slot(matching_seq, "ranges"))
    seq_nuc <- ith_seq[ranges[1]:(ranges[1] + ranges[2] - 1)]
    seq_prot <- seqinr::translate(seq_nuc)
    list(nuc = seq_nuc,
         prot = seq_prot)
  })
}

csga_seqs <- find_sequence(seq_k12[[1]], seqs)
