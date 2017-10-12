library(qpcR)
library(dplyr)
library(reshape2)
library(ggplot2)

dat <- read.csv2("/home/michal/Dropbox/BTU/CFX96_Real-Time_System_qpcr/170917_csgd_csga_fumc_all_strain -  Quantification Amplification Results_FAM.csv") %>% 
  select(-X)

fits <- modlist(dat, model = l4)

#sapply(fits, function(x) Cy0(x))

mdat <-melt(dat, id.vars = "Cycle") %>% 
  inner_join(data.frame(variable = gsub(pattern = "*", replacement = "", x = sapply(fits, function(i) i[["names"]]), fixed = TRUE),
                        model = sapply(fits, function(i) is.null(i[["MODEL"]]))))

ggplot(mdat, aes(x = Cycle, y = value, color = model, group = variable)) +
  geom_line() 

sapply(fits[!(sapply(fits, function(i) is.null(i[["MODEL"]])))], Cy0)

