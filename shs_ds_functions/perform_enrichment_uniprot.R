# belongs to shs_downstream
perform_enrichment_uniprot <- function(query_list, analysis_type = c("GO", "GO_slim", "hallmark"), set) {
  analysis_type <- match.arg(analysis_type) # analysis_type <- "GO"
  
  # Connect to Ensembl BioMart
  ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Convert UniProt IDs to Ensembl gene IDs
  ensembl_ids_query_info <- getBM(attributes = c("uniprotswissprot", "ensembl_gene_id", "hgnc_symbol"),
                   filters    = "uniprotswissprot",
                   values     = query_list,
                   mart       = ensembl_human)  
  
  ensembl_ids_query <- ensembl_ids_query_info%>%
    filter(ensembl_gene_id != "") %>%
    pull(ensembl_gene_id) %>%
    unique()

  if (length(ensembl_ids_query) == 0) {
    stop("No Ensembl IDs found for the provided UniProt IDs.")
  }
  
  ###### define background human protein coding gene background
  swissprot_accessions <- proteome %>% pull(entry) # Your downloaded list
  background_universe_human_protein_coding_genes <- getBM(
    attributes = c("uniprotswissprot", "ensembl_gene_id"),
    filters    = "uniprotswissprot",
    values = swissprot_accessions,
    mart = ensembl_human
  ) %>%
  filter(ensembl_gene_id != "") %>%
  pull(ensembl_gene_id) %>%
  unique()
  paste("fyi proteome ref contains proteins:      ", length(proteome %>% pull(entry)))
  paste("ensembl backgroudn contains annotations: ", length(background_universe_human_protein_coding_genes))
  
  # Perform enrichment based on selected analysis type
  if (analysis_type == "GO") {  # .....................................................................................................
    enrich_res <- enrichGO(gene    = ensembl_ids_query,
                           OrgDb   = org.Hs.eg.db,
                           keyType = "ENSEMBL",
                           ont     = "ALL",  # One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                           universe      = background_universe_human_protein_coding_genes,
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.2 ,
                           minGSSize     = 2,
                           maxGSSize     = 100)
    
  } else if (analysis_type == "GO_slim") {  # ..........................................................................................
    # obtain goslim terms for human protein coding gene background
    term2gene <- getBM(attributes = c("goslim_goa_accession", "goslim_goa_description", "ensembl_gene_id"),   # retrieve GO_slim info 
                       filters    = "ensembl_gene_id",
                       values     = background_universe_human_protein_coding_genes,
                       mart       = ensembl_human) %>%
      filter(goslim_goa_accession != "") %>%
      dplyr::rename(term = goslim_goa_accession, description = goslim_goa_description, gene = ensembl_gene_id)
    
    if (nrow(term2gene) == 0) {
      stop("No GO slim terms found for the provided genes.")
    }

    enrich_res <- enricher(gene          = ensembl_ids_query,
                           TERM2GENE     = term2gene %>% dplyr::select(term, gene) ,
                           universe      = background_universe_human_protein_coding_genes,
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.2 ,
                           minGSSize     = 2,
                           maxGSSize     = 100)
    
    enrich_res@result <- enrich_res@result %>%
      mutate(ONTOLOGY = paste(analysis_type)) %>%   # label analysis type (add collumn to allow easy rbind of results)
      dplyr::select(-Description) %>%
      left_join(                                    # map the term description (lost in enricher())
        term2gene %>% 
          dplyr::select(term, description) %>%
          distinct(term, .keep_all = TRUE),    
        by = c("ID" = "term")  
      ) %>%
      dplyr::rename(Description = description) %>%
      dplyr::select(ONTOLOGY, ID, Description, everything())
    
  } else if (analysis_type == "hallmark") { # ..........................................................................................
    # Fetch Hallmark gene sets
    hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
    # Prepare term2gene (data frame format)
    term2gene <- hallmark_sets %>%
      dplyr::select(gs_name, gs_description, ensembl_gene) %>%
      distinct() %>%
      dplyr::rename(term = gs_name, description = gs_description, gene = ensembl_gene)
    
    if (nrow(term2gene) == 0) {
      stop("No Hallmark terms found.")
    }
    
    # Run enrichment with standardized parameters
    enrich_res <- enricher(
      gene          = ensembl_ids_query,
      TERM2GENE     = term2gene %>% dplyr::select(term, gene),
      universe      = background_universe_human_protein_coding_genes,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 2,
      maxGSSize     = 100
    )
    # Add ontology label for consistency
    enrich_res@result <- enrich_res@result %>%
      mutate(ONTOLOGY = paste(analysis_type)) %>%   # label analysis type (add collumn to allow easy rbind of results)
      dplyr::select(-Description) %>%
      left_join(                                    # map the term description (lost in enricher())
        term2gene %>% 
          dplyr::select(term, description) %>%
          distinct(term, .keep_all = TRUE),  
        by = c("ID" = "term") 
      ) %>%
      dplyr::rename(Description = description) %>%
      dplyr::select(ONTOLOGY, ID, Description, everything())
    
    
  }
  
  # add column to distinguish code iterations
  enrich_res@result <- enrich_res@result %>%  
    mutate(comparison    = paste(set),
           analysis_type = paste(analysis_type))
  
  ## mapping of hgnc gene info
  # Create mapping vector from annotation dataframe
  id_to_symbol <- setNames(
    ensembl_ids_query_info$hgnc_symbol,
    ensembl_ids_query_info$ensembl_gene_id
  )
  # Translate IDs in matched column
  enrich_res@result$geneID_hgnc <- sapply(strsplit(enrich_res@result$geneID, "/"), function(ids) {
    translated <- ifelse(
      ids %in% names(id_to_symbol),
      id_to_symbol[ids],
      ids  # Keep original ID if not found
    )
    paste(translated, collapse = ", ")
  })
  
  
  # asdf <- enrich_res@result
  # View(asdf)
  
  return(enrich_res)
}
