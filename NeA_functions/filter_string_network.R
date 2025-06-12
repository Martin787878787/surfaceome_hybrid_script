filter_string_network <- function(string_network, string_score_cutoff, 
                                  changes, 
                                  network_output_filename){

 ### new uniprot to string code 2025-05-05 Martin (see old code below) #########################################################################################################################
 ###############################################################################################################################################################################################

  ## step 0: check if required string resource exists - if not download ---------------------------------------------------------------------------------------------------------------------------------------------------------
  file_path <- paste0("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/_upsp_string_",  string_network, "_cutoff", string_score_cutoff, ".csv"  ) # Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/_upsp_string_
  
  # Check if the file exists
  if (!file.exists(file_path)) {
    message(
      " ------------------------------------------------------------------------------------------------------------------------------ \n",
      "   --- ATTENTION - local string resource with specified parameters does not exist. will download it. will take 5-10 minutes  ---- \n",
      "   ------------------------------------------------------------------------------------------------------------------------------ "
    )
    ## step 0.1: source human proteome. used to query string & key to translate ENSP column to entry format ---------------------------------------------------------------------------------
    library(queryup)
    up_human <- get_uniprot_data(
      query = list(
        organism_id = "9606",
        reviewed = "true"
      ),
      columns = c("accession", "id", "gene_primary", "xref_string", "reviewed")   # "xref_ensembl",
    )$content %>%
      rename("entry" = "Entry", "entry_name" = "Entry Name", "gene" = "Gene Names (primary)", "string" = "STRING", "reviewed" = "Reviewed") %>%
      mutate(string = gsub(";", "", string)) 

    write.csv(up_human,  "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/_NeA_uniprot_string_ensembl_key.csv", row.names = FALSE)  # /Users/mgesell/Desktop/currentR/git/shs_resources/_NeA_uniprot_string_ensembl_key.csv
    #
    up_human <- read.csv("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/_NeA_uniprot_string_ensembl_key.csv", header = TRUE)     # /Users/mgesell/Desktop/currentR/git/shs_resources/_NeA_uniprot_string_ensembl_key.csv
    paste( # FYI
      "Out of", length(unique(up_human$entry)), "proteins, only",
      length(unique(up_human$string)), "string annotated.",
      length(unique(up_human$entry)) - length(unique(up_human$string)), "missing (",
      round((length(unique(up_human$entry)) - length(unique(up_human$string))) / length(unique(up_human$entry)) * 100, 1), "% missing)"
    )
    
    ## step 0.2: string - download, filter and translate to entry format  -----------------------------------------------------------------------------------
    library(STRINGdb)
    # Initialize STRINGdb instance with human species and score threshold 400
    string_db <- STRINGdb$new(
      version = "12.0",
      species = 9606,                          # NCBI taxonomy ID for Homo sapiens
      score_threshold = string_score_cutoff,   # Includes interactions_string_up_human with combined score ≥ ...
      network_type    = string_network         # Includes functional and physical interactions_string_up_human
    )
    # Map UniProt IDs to string identifiers
    string_up_human <- string_db$map(
      data.frame(uniprot_id = up_human$entry),
      "uniprot_id",
      removeUnmappedRows = TRUE  # Remove unmappable IDs
    )
    # Get ALL interactions_string_up_human between string_up_human proteins and their partners
    # and filter
    interactions_string_up_human <- string_db$get_interactions(string_up_human$STRING_id) %>%
      filter(!is.na(from), !is.na(to)) %>%  # filter out NAs in mapping ------
    filter(from != to) %>%                # filter out self-loops  A = A --
      mutate(protein_min = pmin(from, to),
             protein_max = pmax(from, to)
      ) %>%
      distinct(protein_min, protein_max, .keep_all = TRUE) %>%   # kick AB BA duplicates ----------
    dplyr::select(-protein_min, -protein_max) %>%
      distinct() # Remove duplicates
    
    # Map 'from' and 'to' columns in interactions_string_up_human to UniProt entry
    interactions_string_up_human <- interactions_string_up_human %>%
      left_join(up_human, by = c("from" = "string")) %>%
      rename(protein1_entry  = entry) %>%
      left_join(up_human, by = c("to" = "string")) %>%
      rename(protein2_entry  = entry)
    
    upsp_string <- interactions_string_up_human %>%
      dplyr::select(protein1_entry, protein2_entry) %>%
      distinct()
    write.csv(upsp_string, paste0("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/_upsp_string_", string_network, "_cutoff", string_score_cutoff, ".csv"),  # /Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/_upsp_string_
              row.names = FALSE)
    
  } # _____________________________________________________________________________________________________________________________________________________________________________________
  #########################################################################################################################################################################################
  #########################################################################################################################################################################################

  
  # read string resource file with specified parameters
  upsp_string <- read.csv(file_path, header = TRUE) %>% 
    filter(!is.na(protein1_entry) & !is.na(protein2_entry))
  ## step 2: subset string ---------------------------------------------------------------------------------------------------------------------------------------------------------
  string_output <- upsp_string %>% 
    filter(protein1_entry %in% changes | protein2_entry %in% changes) %>% # at least one match with changes list (point is to add +1 interaction layer later)
    distinct()
  
  ## step 3: export -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  write.table(string_output, paste(result_directory, "intermediate/", network_output_filename, ".tsv", sep = "" ), quote = F , sep = "\t", row.names = F)
  
  ## step 4: reporting  --------------------------------------------------------------------------------------------------------------------------------------------------------------
  message(paste0("_____ selected string resource (stringency filters) contains ppi info for ", length(unique(c(upsp_string$protein1_entry, upsp_string$protein2_entry))) , " unique proteins ______"))
  message(paste0("_____ string mapping of ",  length(changes) , " proteins yielded ", nrow(string_output), " unique PPIs of ",  length(c(unique(string_output$protein1_entry, string_output$protein2_entry))), " unique proteins ______"))
}





# §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
# §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
# previous code §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

# filter_string_network <- function(string_network, string_score_to_use,  string_score_cutoff,
#                                   changes,
#                                   network_output_filename){
  # # load data -----------------------------------------------------------------------------------------------------------------------------------------------------------
  # if (string_network == "uniprot_unreviewed") {  
  #   string_input <-   fread("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/string_net_unrev.csv")  # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_unrev.csv
  # 
  # } else if (string_network == "uniprot_swissprot") {
  #   string_input <-  fread("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/string_net_upsp.csv")   # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_upsp.csv
  # } 
  # 
  # ## string df filtering ------------------------------------------------------------------------------------------------------------------------------------------------
  # # filter score
  # string_output <- string_input %>%
  #   dplyr::select(protein1_entry, protein2_entry, protein1_entry_name, protein2_entry_name, !!sym(string_score_to_use) ) %>%
  #   filter(!!sym(string_score_to_use) >= string_score_cutoff) %>%
  #   select(-!!sym(string_score_to_use))
  # 
  # # filter ppi of interst (ppioi) ---------------   = first interaction shell of pois = ppioi +1 (all ppis that involve at lease 1 of pois)
  # string_output <- string_output %>%
  #   filter(protein1_entry %in% changes | protein2_entry %in% changes)
  # 
  # # filter ppi NA / AB-BA / self-loop -----------
  # string_output <- string_output %>% 
  #   ##   filter out NAs in mapping ------
  #   filter(!is.na(protein1_entry), !is.na(protein2_entry)) %>%
  #   ##   filter out self-loops  A = A --
  #   filter(protein1_entry != protein2_entry) %>%
  #   # # kick AB BA duplicates ----------
  #   mutate(protein_min = pmin(protein1_entry, protein2_entry),
  #          protein_max = pmax(protein1_entry, protein2_entry)
  #   ) %>%
  #   distinct(protein_min, protein_max, .keep_all = TRUE) %>%
  #   dplyr::select(-protein_min, -protein_max) 
  #   
  #   ##   reduce to columns of interest
  # string_output <- string_output[, c("protein1_entry", "protein2_entry")]
# }