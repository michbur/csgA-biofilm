library(msa)
library(dplyr)
library(ggplot2)

# add reference gene K12
# confirm border of regions in region_border object

mySeqs <- readAAStringSet("./data/CsgA_prot.fasta")

aln <- msa(mySeqs, method = "Muscle")
aln_m <- msaConvert(aln, "seqinr::alignment") %>% 
  as.vector %>% 
  getElement("seq") %>% 
  strsplit("") %>% 
  do.call(rbind, .)

diff_pos <- which(apply(aln_m, 2, function(i) length(unique(i))) > 1)

region_borders <- list(R1 = 43L:65, 
                       #R2 = 66L:87,
                       R3 = 88L:110,
                       #R4 = 111L:132,
                       R5 = 133L:151) %>% 
  unlist

lapply(diff_pos, function(ith_row)
  data.frame(table(aln_m[, ith_row]), 
             region = ith_row %in% region_borders,
             pos = paste0("Position: ", ith_row), stringsAsFactors = FALSE) %>% 
    mutate(Freq = Freq/sum(Freq))) %>% 
  bind_rows() %>% 
  mutate(pos = factor(pos, levels = unique(pos))) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = region)) +
  geom_col() +
  scale_x_discrete("Amino acid") +
  scale_y_continuous("Frequency") +
  scale_fill_manual("Amyloid-specific region", values = c("#d8b365", "#5ab4ac")) +
  facet_wrap(~ pos, scales = "free_x", nrow = 3) +
  theme_bw() 
