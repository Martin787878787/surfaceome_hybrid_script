filter_biogrid_network <- function(biogrid_input_file, changes, network_output_filename) {
  
  ##### biogrid = https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.242/BIOGRID-ALL-4.4.242.mitab.zip  ###################################################################################################################################################
  ## biogrid data type = ppi dataset (protein-protein relationship)

  # #code to get biogrid into right format
  # biogrid_human  <- read_protti(biogrid_input_file) %>%
  #   filter(taxid_interactor_a == "taxid:9606" & taxid_interactor_b == "taxid:9606" ) %>%
  #   # allow only phyiscal / direct interactions (= complex partners or direct binding partner)
  #   filter(interaction_types %in% c("psi-mi:\"MI:0407\"(direct interaction)", "psi-mi:\"MI:0915\"(physical association)")) %>%                 # unique(biogrid_human$interaction_types) to display categories
  #   # exclude: "imaging technique" (biochemical)" "(enzymatic study)" "(genetic interference)" "(unspecified method)" "(protein complementation assay)"
  #   filter(interaction_detection_method %in% c("psi-mi:\"MI:0018\"(two hybrid)", "psi-mi:\"MI:0055\"(fluorescent resonance energy transfer)",  # unique(biogrid_human$interaction_detection_method)
  #                                   "psi-mi:\"MI:1313\"(bioid)"     , "psi-mi:\"MI:0047\"(far western blotting)",
  #                                   "psi-mi:\"MI:0004\"(affinity chromatography technology)", "psi-mi:\"MI:0030\"(cross-linking study)",
  #                                   "psi-mi:\"MI:0096\"(pull down)", "psi-mi:\"MI:0114\"(x-ray crystallography)")
  #          ) %>%
  #   dplyr::select(interaction_identifiers, publication_identifiers, interaction_detection_method, interaction_types,  alt_i_ds_interactor_a, alt_i_ds_interactor_b) %>%
  #   # part below takes long time command takes a long time
  #   mutate(protein1_entry = str_extract(alt_i_ds_interactor_a, "(?<=swiss-prot:)[^|]+(?=\\|refseq)"),
  #          protein2_entry = str_extract(alt_i_ds_interactor_b, "(?<=swiss-prot:)[^|]+(?=\\|refseq)")) %>%
  #   dplyr::select(-c(alt_i_ds_interactor_a, alt_i_ds_interactor_b)) %>%
  #   select(protein1_entry, protein2_entry) %>%
  #   # ...
  #   filter(!is.na(protein1_entry), !is.na(protein2_entry)) %>%  # filter out NAs in mapping ------
  #   filter(protein1_entry != protein2_entry)                    # filter out self-loops  A = A --
  # 
  #  write.csv(biogrid_human, "/Users/mgesell/PhD/local_resources/PPIs_and_complexes/_biogrid_human_filtered_method_and_type.csv", row.names = FALSE)
  #  # ______________________________________________________________________________________________________________________________________________________

  # determine presence of PPI in dataset
  biogrid_physical <- read.csv("/Users/mgesell/PhD/local_resources/PPIs_and_complexes/_biogrid_human_filtered_method_and_type.csv", header = TRUE)
  
  biogrid_physical_inverted_dummy <- biogrid_physical %>%
    mutate(dummy1 = protein1_entry,
           dummy2 = protein2_entry,
           protein1_entry = dummy2, # invert annotation
           protein2_entry = dummy1) %>%
    select(-dummy1, -dummy2)   

  biogrid_phys_query_subset <- rbind(biogrid_physical_inverted_dummy, biogrid_physical) %>%
    distinct() %>%
    group_by(protein1_entry) %>%
    mutate(
      interactors_biogrid_phys_global = paste(unique(protein2_entry), collapse = ", "),
      node_degree_biogrid_global      = n_distinct(protein2_entry)
    ) %>%
    ungroup() %>% 
    ## filtering
    # check which ppis are present - values can range 0-2  (caluclate number query_entries that occur per protein1_entry and protein2_entry column)
    mutate(ppi = rowSums(across(starts_with("protein"), ~.x %in% changes))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi == 1) %>% # filter for at least one interactors identified 
    # filter(ppi == 2) %>% # filter for both interactors identified 
    dplyr::select(-ppi) %>%
    distinct() %>%
    group_by(protein1_entry) %>%
    mutate(
      interactors_biogrid_phys_changes = paste(unique(protein2_entry), collapse = ", "),
      node_degree_biogrid_changes      = n_distinct(protein2_entry)
    ) %>%
    ungroup() %>% 
    mutate(protein_min = pmin(protein1_entry, protein2_entry),
           protein_max = pmax(protein1_entry, protein2_entry)
    ) %>%
    distinct(protein_min, protein_max, .keep_all = TRUE) %>%   # kick AB BA duplicates ----------
    dplyr::select(-protein_min, -protein_max) %>%
    distinct() # Remove duplicates
  
  boxplot(
    biogrid_phys_query_subset$node_degree_biogrid_global,
    biogrid_phys_query_subset$node_degree_biogrid_changes,
    names = c("Full", "Changes Subnetwork "),
    main = "Phyiscal filtered Biogrid: Node Degrees"
  )  
  
  ## export -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  write.table(biogrid_phys_query_subset %>% select(protein1_entry, protein2_entry)
              , paste(result_directory, "intermediate/", network_output_filename, ".tsv", sep = "" ), quote = F , sep = "\t", row.names = F)
  
  ## step 4: reporting  --------------------------------------------------------------------------------------------------------------------------------------------------------------
  message(paste0("_____ using physical biogrid network of ", length(  unique( c(biogrid_phys_query_subset$protein1_entry, biogrid_phys_query_subset$protein2_entry) )  )  , " unique proteins and ", nrow(biogrid_phys_query_subset), " ppis ________ "))

  
}