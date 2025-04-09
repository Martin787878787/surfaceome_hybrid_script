filter_string_network <- function(string_network, string_score_to_use,  string_score_cutoff,
                                  changes,
                                  network_output_filename){
  
  # load data -----------------------------------------------------------------------------------------------------------------------------------------------------------
  if (string_network == "uniprot_unreviewed") {  
    string_input <-   fread("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/string_net_unrev.csv")  # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_unrev.csv

  } else if (string_network == "uniprot_swissprot") {
    string_input <-  fread("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/string_net_upsp.csv")   # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_upsp.csv
  } 
  
  ## string df filtering ------------------------------------------------------------------------------------------------------------------------------------------------
  # filter score
  string_output <- string_input %>%
    dplyr::select(protein1_entry, protein2_entry, protein1_entry_name, protein2_entry_name, !!sym(string_score_to_use) ) %>%
    filter(!!sym(string_score_to_use) >= string_score_cutoff) %>%
    select(-!!sym(string_score_to_use))
  
  # filter ppi of interst (ppioi) ---------------   = first interaction shell of pois = ppioi +1 (all ppis that involve at lease 1 of pois)
  string_output <- string_output %>%
    filter(protein1_entry %in% changes | protein2_entry %in% changes)
  
  # filter ppi NA / AB-BA / self-loop -----------
  string_output <- string_output %>% 
    ##   filter out NAs in mapping ------
    filter(!is.na(protein1_entry), !is.na(protein2_entry)) %>%
    ##   filter out self-loops  A = A --
    filter(protein1_entry != protein2_entry) %>%
    # # kick AB BA duplicates ----------
    mutate(protein_min = pmin(protein1_entry, protein2_entry),
           protein_max = pmax(protein1_entry, protein2_entry)
    ) %>%
    distinct(protein_min, protein_max, .keep_all = TRUE) %>%
    dplyr::select(-protein_min, -protein_max) 
    
    ##   reduce to columns of interest
  string_output <- string_output[, c("protein1_entry", "protein2_entry")]
  
  # export -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  write.table(string_output, paste(result_directory, "intermediate/", network_output_filename, ".tsv", sep = "" ), quote = F , sep = "\t", row.names = F)
  
}