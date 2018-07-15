library(dplyr)

if(Sys.info()[["nodename"]] == "amyloid")
  seq_path <- "/home/michal/Dropbox/artykuły/PSI-blast/Bacteria (taxid:2)/csgA/iteracja_5"


all_lines <- readLines(paste0(seq_path, "/Sequences producing significant alignments with E-value BETTER than threshold.txt"))

prot_name_id <- grepl("^>", all_lines)

removed_species <- all_lines[prot_name_id] %>% 
  strsplit("[", fixed = TRUE) %>% 
  sapply(first) %>% 
  strsplit(" ")

protein_descriptions <- lapply(removed_species, function(i) tolower(i[-1]))

names(protein_descriptions) <- gsub(">", "", sapply(removed_species, first))

description_points <- read.table("./data/CsgA-points.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)

all_points <- lapply(1L:nrow(description_points), function(ith_description_name_id)
  sapply(protein_descriptions, function(ith_description) {
    ifelse(any(grepl(description_points[ith_description_name_id, "name"], ith_description)),
           description_points[ith_description_name_id, "points"], 0)
  })
) %>% 
  do.call(cbind, .) %>% 
  rowSums() %>% 
  sort(decreasing = TRUE)
