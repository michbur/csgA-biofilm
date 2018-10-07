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

prot_id <- cumsum(grepl("^>", all_lines))

all_prots <- split(all_lines, prot_id)

aln_dat <- lapply(all_prots, function(ith_prot) {
  strsplit(paste0(ith_prot[-1], collapse = ""), "")[[1]]
}) %>% 
  do.call(rbind, .)

real_positions <- cumsum(aln_dat[1, ] != "-")

all_names <- sapply(all_prots, function(ith_prot) ith_prot[[1]])
aln_dat[, real_positions %in% region_borders[[1]]]


ggplot() + geom_logo( apply(aln_dat[, real_positions %in% region_borders[[1]]], 1, paste0, collapse = "") ) + theme_logo()

aln_dat[, real_positions %in% region_borders[[1]]]

library(AmyloGram)
AmyloGram_model[["enc"]]
# use with make_col_scheme
