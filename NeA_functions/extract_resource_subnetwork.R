extract_resource_subnetwork <- function(ppi_network, pageRank_retain_all_changes){
  
  
  # load the network resource (full network based on specified parameters)
  network_interactions <- read.table(paste(result_directory, "intermediate/", "network_interactions_", ppi_network, ".tsv", sep = "" ), header = TRUE)
  
  # load the proteins that survived the pageRank filtering (depending on pageRank_retain_all_changes this will have produced different lists of proteins)
  if (pageRank_retain_all_changes == FALSE) {
    propagated_proteins <- read.table( paste0(result_directory, "intermediate/propagated_proteins.tsv")                   ,  header = TRUE)
  } else if (pageRank_retain_all_changes == TRUE) {
    propagated_proteins <- read.table( paste0(result_directory, "intermediate/propagated_proteins_plus_changes_input.tsv"),  header = TRUE)
  }
  
  # filter for pageRank survivor proteins
  network_interactions_pageRank_filtered <- network_interactions %>%
    filter(protein1_entry %in% propagated_proteins$x) %>%
    distinct()
  
  ## translate entry name format to desired output (such that in cytoscape or wherever downstream readable names appear)
  # load uniprot
  up_human <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/_NeA_uniprot_string_ensembl_key.csv", header = TRUE) %>%
    select(-string, -reviewed) %>%
    mutate(entry_name = gsub("_HUMAN", "", entry_name))
  
  network_interactions_pageRank_filtered <- network_interactions_pageRank_filtered %>% 
    left_join(up_human %>% select(entry, entry_name, gene), 
              by = c("protein1_entry" = "entry")) %>%
    rename("protein1_entry_name" = "entry_name", "protein1_gene" = "gene") %>% 
    left_join(up_human %>% select(entry, entry_name, gene), 
              by = c("protein2_entry" = "entry")) %>%
    rename("protein2_entry_name" = "entry_name", "protein2_gene" = "gene") %>%
    distinct() %>%
    select(protein1_entry, protein2_entry, protein1_entry_name, protein2_entry_name, protein1_gene, protein2_gene)
  
   # export result network (load into cytoscape via file > import > import from file > then select source and target column)
  write.table(network_interactions_pageRank_filtered  , paste0(result_directory, "intermediate/network_interactions_pageRank_filtered.tsv")     , quote = F , sep = "\t", row.names = F)
  
  # write.table(network_interactions_pageRank_filtered %>% select(protein1_entry, protein2_entry)          , paste0(result_directory, "intermediate/network_interactions_pageRank_filtered_entry.tsv")     , quote = F , sep = "\t", row.names = F)
  # write.table(network_interactions_pageRank_filtered %>% select(protein1_entry_name, protein2_entry_name), paste0(result_directory, "intermediate/network_interactions_pageRank_filtered_entry_name.tsv"), quote = F , sep = "\t", row.names = F)
  # write.table(network_interactions_pageRank_filtered %>% select(protein1_gene, protein2_gene)            , paste0(result_directory, "intermediate/network_interactions_pageRank_filtered_gene.tsv")      , quote = F , sep = "\t", row.names = F)
  
}
