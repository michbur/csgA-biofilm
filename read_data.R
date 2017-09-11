library(dplyr)
library(readxl)
library(reshape2)

dat_vs <- read_xlsx("/home/michal/Dropbox/BTU/staż_BTU/Wyniki_excel/11-09-17 VS-MB.xlsx") %>% 
  melt(variable.name = "replicate") %>% 
  mutate(replicate = factor(replicate, levels = rev(levels(replicate))))
