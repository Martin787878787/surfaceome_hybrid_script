# belongs to shs_downstream
go_gost <- function(query_list, max_term_size, set) {
  # function requires ensembl
  # ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  
  ####### backbone dummy
  result <- gost(query = gconvert(query_list, #setdiff(data_CSC_w_ENSG_table$Act_0, data_CSC_w_ENSG_table$Act_2h),  # this works best, uniprot identifiers worst, primary_gene_names good
                                  numeric_ns = "",
                                  mthreshold = Inf,
                                  filter_na = TRUE) %>% 
                   pull(target), 
                 organism = "hsapiens",
                 ordered_query = FALSE,
                 multi_query = FALSE,
                 significant = TRUE, # <<- filtering outside of function!s
                 exclude_iea = FALSE,
                 measure_underrepresentation = FALSE,
                 evcodes = TRUE,  # TRUE to show which protein are in intersection
                 user_threshold = 0.05,
                 correction_method = "g_SCS",
                 domain_scope = "annotated",
                 # custom_bg = background_upUlrich,
                 numeric_ns = "",
                 sources = NULL, #When sources is set to NULL (the default), the analysis includes all Gene Ontology branches (GO:BP, GO:MF, GO:CC), KEGG pathways, Reactome (REAC), transcription factors (TF), miRNA targets (MIRNA), Human Phenotype Ontology (HP), CORUM protein complexes, WikiPathways (WP), and tissue-specific annotations (HPA).
                 as_short_link = FALSE,
                 highlight = TRUE)
  # subset and trim result
  result_result  <- result$result  %>% 
    filter(term_size < max_term_size) %>%
    mutate(comparison = paste(set)) %>%
    dplyr::select(term_name, comparison, term_size, p_value, recall, intersection) %>%
    arrange(desc(recall)) 
  
  # translate Ensembl IDs into gene name ans protein names
  mapping <- getBM(
    attributes = c(
      "ensembl_gene_id",       # Ensembl Gene ID
      "external_gene_name",    # Gene Name
      "uniprot_gn_id",         # UniProt Accession (Identifier)
      "uniprot_gn_symbol"      # UniProt Entry Name
    ),
    filters = "ensembl_gene_id",
    values  = unique(unlist(strsplit(result_result$intersection, ","))),
    mart    = ensembl
  )
  
  # Define the lookup function ----- 
  lookup_mapping <- function(ens_ids, mapping, column_to_retrieve) {
    # Always match with uniprot_gn_id
    uniprot_ids <- mapping$uniprot_gn_id[match(ens_ids, mapping$ensembl_gene_id)]
    # Retrieve the requested column values
    retrieved_values <- mapping[[column_to_retrieve]][match(ens_ids, mapping$ensembl_gene_id)]
    # Combine the uniprot_gn_id with the retrieved values
    combined_values <- paste(retrieved_values, sep = ",")
    # Remove any NA entries and collapse into a single string
    paste(combined_values[!is.na(uniprot_ids) & !is.na(retrieved_values)], collapse = ",")
  }
  
  # Define the apply_mapping function
  apply_mapping <- function(data, mapping, column_to_retrieve, new_column_name) {
    # Split the intersection column
    data$intersection_split <- strsplit(data$intersection, ",")
    
    # Apply the lookup function
    data[[new_column_name]] <- sapply(
      data$intersection_split, 
      lookup_mapping, 
      mapping = mapping, 
      column_to_retrieve = column_to_retrieve
    )
    
    # Remove the temporary split column
    data$intersection_split <- NULL
    
    return(data)
  }
  
  # specify which columns to translate
  result_result <- apply_mapping(result_result, mapping, column_to_retrieve = "external_gene_name", new_column_name = "intersection_gene")
  result_result <- apply_mapping(result_result, mapping, column_to_retrieve = "uniprot_gn_id"     , new_column_name = "intersection_up")
  
  
  return(result_result) 
}