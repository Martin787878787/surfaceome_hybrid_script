go_children_table_generator2 <- function(my_terms, only_new = TRUE) {
  library(GO.db)
  
  out_file <- "/Users/mgesell/Desktop/currentR/git/shs_resources/GO_term_children/GO_term_children.csv"
  
  # Read the existing table or create an empty one (first row is special!)
  if (file.exists(out_file)) {
    GO_descendants_df <- read.csv(out_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    GO_descendants_df <- as.data.frame(matrix(nrow = 0, ncol = 0))
  }
  
  # The first row of each column contains its GO ID (match to my_terms)
  if (nrow(GO_descendants_df) >= 1) {
    existing_ids <- as.character(GO_descendants_df[1, ])  # it's a named vector
  } else {
    existing_ids <- character(0)
  }
  
  # Decide which terms are new
  if (only_new) {
    terms_to_search <- my_terms[ !my_terms %in% existing_ids ]
    if (length(terms_to_search) == 0) {
      message("All terms already present. Nothing to do.")
      return(GO_descendants_df)
    }
  } else {
    terms_to_search <- my_terms
  }
  
  # Helper to get all children, including self, for each term
  get_all_children <- function(go_id) {
    ontology <- Ontology(go_id)[[1]]
    children_env <- switch(
      toupper(ontology),
      BP = GOBPCHILDREN,
      MF = GOMFCHILDREN,
      CC = GOCCCHILDREN,
      stop(paste("Invalid ontology for", go_id))
    )
    rec_children <- function(term) {
      kids <- unlist(mget(term, children_env, ifnotfound = NA))
      kids <- kids[!is.na(kids)]
      if (length(kids) == 0) return(NULL)
      c(kids, unlist(lapply(kids, rec_children)))
    }
    unique(c(go_id, rec_children(go_id)))
  }
  
  # For each new term, get descendants
  descendants_list <- lapply(terms_to_search, get_all_children)
  
  # Prepare and name each column: first row = GO ID, column name = Term(GOTERM[go_id])
  term_names <- as.character(Term(GOTERM[terms_to_search]))
  term_names[ is.na(term_names) ] <- terms_to_search[ is.na(term_names) ]  # fallback to ID
  
  # Each column: first row is GO ID, remaining are descendants
  max_length <- max(sapply(descendants_list, length)) + 1  # +1 for the ID in first row
  make_col <- function(go_id, go_descs) {
    col <- c(go_id, setdiff(go_descs, go_id))
    length(col) <- max_length
    col
  }
  descendants_mat <- mapply(
    FUN = make_col,
    go_id = terms_to_search,
    go_descs = descendants_list,
    SIMPLIFY = FALSE
  )
  descendants_df_new <- as.data.frame(descendants_mat, stringsAsFactors = FALSE)
  colnames(descendants_df_new) <- term_names
  
  # Bind new columns to the existing data, padding rows with NA
  new_nrows <- max(nrow(GO_descendants_df), nrow(descendants_df_new))
  extend_df <- function(df, n_rows) {
    if (nrow(df) < n_rows) {
      toadd <- matrix(NA, ncol = ncol(df), nrow = n_rows - nrow(df))
      colnames(toadd) <- colnames(df)
      df <- rbind(df, toadd)
    }
    df
  }
  if (ncol(GO_descendants_df) == 0) {
    final_df <- extend_df(descendants_df_new, new_nrows)
  } else {
    GO_descendants_df <- extend_df(GO_descendants_df, new_nrows)
    descendants_df_new <- extend_df(descendants_df_new, new_nrows)
    final_df <- cbind(GO_descendants_df, descendants_df_new)
  }
  
  write.csv(final_df, out_file, row.names = FALSE)
  return(final_df)
}
