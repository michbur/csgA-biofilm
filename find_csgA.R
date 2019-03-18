library(biogram)
library(seqinr)
library(dplyr)

seqs <- read_fasta("data/csg_all.fasta")
seq_k12_csga <- read_fasta("data/CsgA-K12-gene.fasta")
seq_k12_csgb <- read_fasta("data/CsgB-K12-gene.fasta")
seq_k12_csgc <- read_fasta("data/CsgC-K12-gene.fasta")
seq_k12_csgd <- read_fasta("data/CsgD-K12-gene.fasta")
seq_k12_csge <- read_fasta("data/CsgE-K12-gene.fasta")
seq_k12_csgf <- read_fasta("data/CsgF-K12-gene.fasta")
seq_k12_csgg <- read_fasta("data/CsgG-K12-gene.fasta")

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


csga_seqs <- find_sequence(seq_k12_csga[[1]], seqs)
csgb_seqs <- find_sequence(seq_k12_csgb[[1]], seqs)
csgc_seqs <- find_sequence(seq_k12_csgc[[1]], seqs)
csgd_seqs <- find_sequence(seq_k12_csgd[[1]], seqs)
csge_seqs <- find_sequence(seq_k12_csge[[1]], seqs)
csgf_seqs <- find_sequence(seq_k12_csgf[[1]], seqs)
csgg_seqs <- find_sequence(seq_k12_csgg[[1]], seqs)


biogram::write_fasta(lapply(csga_seqs, function(i) i[["prot"]][i[["prot"]] != "*"]), file = "data/CsgA_prot.fasta")
biogram::write_fasta(lapply(csgb_seqs, function(i) i[["prot"]]), file = "data/CsgB_prot.fasta")
biogram::write_fasta(lapply(csgc_seqs, function(i) i[["prot"]]), file = "data/CsgC_prot.fasta")
biogram::write_fasta(lapply(csgd_seqs, function(i) i[["prot"]]), file = "data/CsgD_prot.fasta")
biogram::write_fasta(lapply(csge_seqs, function(i) i[["prot"]]), file = "data/CsgE_prot.fasta")
biogram::write_fasta(lapply(csgf_seqs, function(i) i[["prot"]]), file = "data/CsgF_prot.fasta")
biogram::write_fasta(lapply(csgg_seqs, function(i) i[["prot"]]), file = "data/CsgG_prot.fasta")

biogram::write_fasta(lapply(csga_seqs, function(i) i[["nuc"]]), file = "data/CsgA_nuc.fasta")
biogram::write_fasta(lapply(csgb_seqs, function(i) i[["nuc"]]), file = "data/CsgB_nuc.fasta")
biogram::write_fasta(lapply(csgc_seqs, function(i) i[["nuc"]]), file = "data/CsgC_nuc.fasta")
biogram::write_fasta(lapply(csgd_seqs, function(i) i[["nuc"]]), file = "data/CsgD_nuc.fasta")
biogram::write_fasta(lapply(csge_seqs, function(i) i[["nuc"]]), file = "data/CsgE_nuc.fasta")
biogram::write_fasta(lapply(csgf_seqs, function(i) i[["nuc"]]), file = "data/CsgF_nuc.fasta")
biogram::write_fasta(lapply(csgg_seqs, function(i) i[["nuc"]]), file = "data/CsgG_nuc.fasta")
