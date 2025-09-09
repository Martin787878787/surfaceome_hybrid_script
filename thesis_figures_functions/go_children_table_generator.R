go_children_table_generator <- function(my_terms) {
  
  library(GO.db) # to extract GO terms and (IMPORTANT!) all their children
  
  GO_descendants_df     <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/GO_term_children/GO_term_children.csv", header = TRUE , check.names = FALSE)
  
  query_log <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/GO_term_children/query_log.csv", header = TRUE)
  query_log <- rbind(query_log, data.frame(dummy = c("__new query__", my_terms)))
  write.csv(query_log, "/Users/mgesell/Desktop/currentR/git/shs_resources/GO_term_children/query_log.csv", row.names = FALSE)
  
  # Function: get all children of a GO term, *auto-detects* ontology
  get_all_children <- function(go_id) {
    ontology <- Ontology(go_id)[[1]]  # auto-detect for this term
    children_env <- switch(
      toupper(ontology),
      BP = GOBPCHILDREN,
      MF = GOMFCHILDREN,
      CC = GOCCCHILDREN,
      stop(paste("Invalid ontology for", go_id))
    )
    rec_children <- function(term) {
      kids <- unlist(mget(term, children_env, ifnotfound=NA))
      kids <- kids[!is.na(kids)]
      if (length(kids) == 0) return(NULL)
      c(kids, unlist(lapply(kids, rec_children)))
    }
    unique(rec_children(go_id))
  }
  
  # Prepare descendants list for all your terms
  descendants_list <- lapply(my_terms, function(go) unique(c(go, get_all_children(go))))
  
  # Prepare for data frame output
  max_len <- max(sapply(descendants_list, length))
  descendants_mat <- sapply(descendants_list, function(x) { length(x) <- max_len; x })
  GO_descendants_df <- as.data.frame(descendants_mat, stringsAsFactors=FALSE)
  
  # Use GO term names as column names
  colnames(GO_descendants_df) <- Term(GOTERM[my_terms])
  
  # Write to csv
  write.csv(GO_descendants_df, "/Users/mgesell/Desktop/currentR/git/shs_resources/GO_term_children/GO_term_children.csv", row.names=FALSE)
  return(GO_descendants_df)
  
}