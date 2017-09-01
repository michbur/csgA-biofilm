library(dplyr)

raw_dat <- readLines("data/raw_CV1.txt")

gsub(x = raw_dat, pattern = '\\"', replacement = "",) %>% 
  strsplit(split = ",") %>% 
  do.call(rbind, .)

strsplit(raw_dat[-1], ',"') %>% 
  lapply(gsub, pattern = '\\"', replacement = "") %>% 
  lapply(function(i) unlist(strsplit(i, "[:digit:],"))) 


  lapply(gsub, pattern = "7Â", replacement = "37Â") 


