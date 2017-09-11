read_uniprot <- function (connection, ft_names, kwds = NULL)  {
  all_lines <- readLines(connection)
  prot_ids <- grep("\\<ID   ", all_lines)
  fts <- lapply(ft_names, function(single_ft) {
    signalHsmm:::get_ft(all_lines, prot_ids, single_ft)
  })

  all_seqs <- cbind(id = prot_ids, signalHsmm:::get_add_id(all_lines))
  
  list_prots <- lapply(1L:nrow(all_seqs), function(i) {
    start_seq <- all_seqs[i, "seq_start"]
    end_seq <- all_seqs[i, "seq_end"]
    ith_seq <- strsplit(gsub(" ", "", paste0(all_lines[start_seq:end_seq], 
                                             collapse = "")), "")[[1]]
    class(ith_seq) <- "SeqFastaAA"
    aa_name <- strsplit(all_lines[all_seqs[i, "id"]], "   ")[[1]][2]
    attr(ith_seq, "name") <- aa_name
    attr(ith_seq, "Annot") <- paste0(">", aa_name)
    attr(ith_seq, "class") <- "SeqFastaAA"
    os <- all_lines[all_seqs[i, "os_start"]:all_seqs[i, 
                                                       "os_end"]]
    attr(ith_seq, "OS") <- paste0(vapply(strsplit(os, "OS   "), 
                                         function(single_os) single_os[2], "a"), collapse = " ")
    oc <- all_lines[all_seqs[i, "oc_start"]:all_seqs[i, 
                                                       "oc_end"]]
    attr(ith_seq, "OC") <- paste0(vapply(strsplit(oc, "OC   "), 
                                         function(single_oc) single_oc[2], "a"), collapse = " ")
    browser()
    ith_seq
  })
  
  names(list_prots) <- sapply(list_prots, function(i) attr(i, "name"))
  list_prots
}