library(dplyr)
library(reshape2)
library(ggplot2)

all_lines <- readLines("ecoli_csga.aln")

prot_id <- cumsum(grepl("^>", all_lines))

all_prots <- split(all_lines, prot_id)

lapply(all_prots, function(ith_prot) {
  list(name = sub("^>", "", ith_prot[1]),
             aa = strsplit(paste0(ith_prot[-1], collapse = ""), "")[[1]])
})

raw_pim <- read.table("pim.txt", header = FALSE, stringsAsFactors = FALSE) %>% 
  select(-V1) %>% 

colnames(raw_pim) <- c("strain", raw_pim[[1]])

m_pim <- melt(raw_pim) %>% 
  arrange(strain, variable) %>% 
  right_join(combn(x = unique(raw_pim[["strain"]]), m = 2, simplify = TRUE) %>%
               t %>%
               data.frame(),
             by = c("strain" = "X1", "variable" = "X2")) %>% 
  mutate(strain = factor(strain, levels = unique(raw_pim[["strain"]])),
         variable = factor(variable, levels = rev(unique(raw_pim[["strain"]]))))



p1 <- ggplot(m_pim, aes(x = strain, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient("PID", low = "#fc8d59", high = "#99d594") +
  scale_x_discrete("Strain") +
  scale_y_discrete("Strain") +
  theme_bw() +
  theme(axis.text = element_blank())


p2 <- ggplot(m_pim, aes(x = strain, y = variable, fill = cut(value, c(0, 50, 90, 95, 99, 100), include.lowest = TRUE, right = TRUE))) +
  geom_tile(color = "black") +
  scale_fill_discrete("PID") +
  scale_x_discrete("Strain") +
  scale_y_discrete("Strain") +
  theme_bw() +
  theme(axis.text = element_blank())

save(p1, p2, file = "plots.RData")
