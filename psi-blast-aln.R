library(ggplot2)
library(dplyr)
library(ggseqlogo)

# important regions ---------------------------------

region_borders <- list(R1 = 43L:65, 
                       #R2 = 66L:87,
                       R3 = 88L:110,
                       #R4 = 111L:132,
                       R5 = 133L:151)

alph_color_list <- list(list(col = "#8dd3c7", aa = c("G", "A", "S", "T")),
                        list(col = "#ffffb3", aa = c("C", "V", "I", "L", "P", "F", "Y", "M", "W")),
                        list(col = "#bebada", aa = c("N", "H", "Q")),
                        list(col = "#fb8072", aa = c("D", "E")),
                        list(col = "#80b1d3", aa = c("K", "R")),
                        list(col = "white", aa = c("-")))


all_lines <- readLines("/home/michal/Dropbox/dropbox-amylogram/PSI-blast/CsgA_muscle.fas")

all_lines <- readLines("/home/jarek/Dropbox/amyloids/PSI-blast/CsgA_muscle.fas")

prot_id <- cumsum(grepl("^>", all_lines))

all_prots <- split(all_lines, prot_id)

aln_dat <- lapply(all_prots, function(ith_prot) {
  strsplit(paste0(ith_prot[-1], collapse = ""), "")[[1]]
}) %>% 
  do.call(rbind, .)

real_positions <- cumsum(aln_dat[1, ] != "-")

all_names <- sapply(all_prots, function(ith_prot) ith_prot[[1]])
aln_dat[, real_positions %in% region_borders[[1]]]


ggplot() + geom_logo( apply(aln_dat[, real_positions %in% region_borders[[2]]], 1, paste0, collapse = ""), col_scheme = cs2) + theme_logo()

aln_dat[, real_positions %in% region_borders[[1]]]

library(AmyloGram)
AmyloGram_model[["enc"]]
# use with make_col_scheme

am <- AmyloGram_model[["enc"]]
cs2 = make_col_scheme(chars=c("G", 
                              "K", "P", "R", 
                              "I", "L", "V", 
                              "F", "W", "Y", 
                              "A", "C", "H", "M", 
                              "D", "E", "N", "Q", "S", "T")), 
                      groups = c("gr1", "gr2", "gr3", "gr4", "gr5", "gr6")

# Generate sequence logo
ggseqlogo(seqs_dna$MA0001.1, col_scheme=cs2)





# functions - csgA(region) CsgA(1)

csgA <- function(x){
region_borders <- list(R1 = 43L:65, 
                       #R2 = 66L:87,
                       R3 = 88L:110,
                       #R4 = 111L:132,
                       R5 = 133L:151)

# all_lines <- readLines("/home/michal/Dropbox/dropbox-amylogram/PSI-blast/CsgA_muscle.fas")

all_lines <- readLines("/home/jarek/Dropbox/amyloids/PSI-blast/CsgA_muscle.fas")

prot_id <- cumsum(grepl("^>", all_lines))

all_prots <- split(all_lines, prot_id)

aln_dat <- lapply(all_prots, function(ith_prot) {
  strsplit(paste0(ith_prot[-1], collapse = ""), "")[[1]]
}) %>% 
  do.call(rbind, .)

real_positions <- cumsum(aln_dat[1, ] != "-")

all_names <- sapply(all_prots, function(ith_prot) ith_prot[[1]])
aln_dat[, real_positions %in% region_borders[[1]]]

ggplot() + geom_logo( apply(aln_dat[, real_positions %in% region_borders[[x]]], 1, paste0, collapse = "")) + theme_logo()
}



csgB <- function(x){
  region_borders <- list(R1 = 45L:66, 
                         R2 = 67L:88,
                         R3 = 89L:110,
                         R4 = 111L:132,
                         R5 = 133L:154
                         )

  # all_lines <- readLines("/home/michal/Dropbox/dropbox-amylogram/PSI-blast/CsgB_muscle.fas")
  
  all_lines <- readLines("/home/jarek/Dropbox/amyloids/PSI-blast/CsgB_muscle.fas")
  
  prot_id <- cumsum(grepl("^>", all_lines))
  
  all_prots <- split(all_lines, prot_id)
  
  aln_dat <- lapply(all_prots, function(ith_prot) {
    strsplit(paste0(ith_prot[-1], collapse = ""), "")[[1]]
  }) %>% 
    do.call(rbind, .)
  
  real_positions <- cumsum(aln_dat[1, ] != "-")
  
  all_names <- sapply(all_prots, function(ith_prot) ith_prot[[1]])
  aln_dat[, real_positions %in% region_borders[[1]]]
  
  
  ggplot() + geom_logo( apply(aln_dat[, real_positions %in% region_borders[[x]]], 1, paste0, collapse = "")) + theme_logo()
}
