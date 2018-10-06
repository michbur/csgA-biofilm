library(dplyr)
library(rentrez)
library(pbapply)

if(Sys.info()[["nodename"]] %in% c("amyloid", "lori"))
  seq_path <- "/home/michal/Dropbox/dropbox-amylogram/PSI-blast/Bacteria (taxid:2)/csgA/iteracja_5"

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

download_res <- pblapply(names(which(all_points > 5)), function(single_term) try({
  #WP_003858043.1
  prot_id <- entrez_search(db = "protein", term = paste0(single_term, "[Accession]"), 
                           config = httr::config(http_version = 2))
  genomes_links <- entrez_link(dbfrom='protein', id=prot_id[["ids"]], db = "nuccore", config = httr::config(http_version = 2))
  
  genomes <- lapply(genomes_links[["links"]][["protein_nuccore_wp"]], function(ith_id)
    entrez_fetch(db = "nuccore", id = ith_id, rettype = "gb", 
                 retmode = "text", config = httr::config(http_version = 2))) %>% 
    unlist
  
  genomes_path <- paste0("/home/michal/Dropbox/dropbox-amylogram-comp/all-genomes-from-NCBI/", single_term, ".gbk")
  cat(genomes, file = genomes_path)
  single_term
}))
