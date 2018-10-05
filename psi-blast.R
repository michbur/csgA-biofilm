library(dplyr)
library(rentrez)
library(XML)
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

#grep("WP_010429577.1", all_lines[prot_name_id])
res <- pblapply(names(which(all_points > 5)), function(single_term) try({
  #WP_003858043.1
  prot_id <- entrez_search(db = "protein", term = paste0(single_term, "[Accession]"), 
                           config = httr::config(http_version = 2))
  genomes_links <- entrez_link(dbfrom='protein', id=prot_id[["ids"]], db = "nuccore", config = httr::config(http_version = 2))
  
  genomes <- lapply(genomes_links[["links"]][["protein_nuccore_wp"]], function(ith_id)
                    entrez_fetch(db = "nuccore", id = ith_id, rettype = "gb", 
                          retmode = "text", config = httr::config(http_version = 2))) %>% 
    unlist
  
  dir.create(single_term)
  genomes_path <- paste0("./", single_term, "/genomes.gbk")
  cat(genomes, file = genomes_path)
  system(paste0("awk -v n=1 '/^$/{close(\"out-genome\"n);n++;next} {print > \"", single_term,"/out-genome\"n}' ", genomes_path))
  file.remove(genomes_path)
  
  single_genomes <- list.files(paste0("./", single_term), full.names = TRUE)
  
  genome_proteins <- lapply(single_genomes, function(single_genome) try({
    
    genome_lines <- readLines(single_genome)
    protein_line_id <- grep(single_term, genome_lines)
    
    # get the name of the source organism
    
    source_name <- genome_lines[grep("DEFINITION  ", genome_lines):(grep("ACCESSION  ", genome_lines) - 1)] %>% 
      gsub(pattern = "DEFINITION  ", replacement = "") %>% 
      gsub(pattern = "[ ]{2,}", replacement = "") %>% 
      paste0(collapse = " ")
    
    
    # identify line defining the coordinates of CsgA cds
    cds_line_id <- protein_line_id
    
    
    while(!grepl("CDS", genome_lines[cds_line_id])) {
      cds_line_id <- cds_line_id - 1
    }
    
    # convert start coordinate to numeric
    csgA_start <- strsplit(genome_lines[cds_line_id], split = "..", fixed = TRUE)[[1]][[1]] %>% 
      gsub(pattern = "[^0-9]", replacement = "") %>% 
      as.numeric
    
    all_cds_id <- grep("     CDS             ", genome_lines)
    
    all_cds_start <- genome_lines[all_cds_id] %>% 
      gsub(pattern = "[^0-9\\.]", replacement = "") %>% 
      strsplit(split = "..", fixed = TRUE) %>% 
      sapply(first) %>% 
      as.numeric
    
    # choose only CDSs that start at least 5000 bp before CsgA
    # and at most 7000 bp after CsgA
    line_id_df <- data.frame(start = all_cds_start, line_id = all_cds_id) %>% 
      filter(start > (csgA_start - 5000) & start < (csgA_start + 7000))
    
    # CDS of the interest - near CsgA
    coi <- genome_lines[min(line_id_df[["line_id"]]):max(line_id_df[["line_id"]])]
    
    # probable ends of cds in the gbk file
    putative_cds_id_end <- (which(!grepl("                     ", coi)) - 1)[-1]
    
    all_cds_in_coi <- lapply(grep("     CDS             ", coi), function(i) {
      end_id <- which.max(i - putative_cds_id_end[putative_cds_id_end > i])
      if(length(end_id) > 0) {
        i:(putative_cds_id_end[putative_cds_id_end > i])[end_id]
      } else {
        NULL
      }
    })
    
    lapply(all_cds_in_coi[lengths(all_cds_in_coi) != 0], function(single_cds_id) try({
      single_cds <- coi[single_cds_id]
      
      single_prot_seq <- single_cds[grep("translation=\"", single_cds):length(single_cds)] %>% 
        gsub(pattern = '[ "=/a-z]', replacement = "")  %>% 
        strsplit("") %>% 
        unlist
      
      single_prot_id <- single_cds[grep("/protein_id=\"", single_cds)] %>% 
        gsub(pattern = '[ "=/]', replacement = "") %>% 
        gsub(pattern = "protein_id", replacement = "")
      
      single_prot_name <- single_cds[grep("/product=\"", single_cds)] %>% 
        gsub(pattern = '["=/]', replacement = "") %>% 
        gsub(pattern = "product", replacement = "") %>% 
        gsub(pattern = "[ ]{2,}", replacement = "")
      
      
      list(definition = source_name,
           id = single_prot_id,
           name = single_prot_name,
           seq = single_prot_seq)  
    }, silent = TRUE)
    )
    
  }, silent = TRUE)) 
  
  unlink(single_term, recursive = TRUE)
  
  genome_proteins
}, silent = TRUE)) 

save(res, file = "/home/michal/Dropbox/dropbox-amylogram/PSI-blast/NCBI-Csg.RData")
