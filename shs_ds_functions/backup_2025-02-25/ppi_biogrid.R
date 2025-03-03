# belongs to shs_downstream
ppi_biogrid <- function(query_list, mode, set) {
  #################################################################################################################################################################################################################
  ##### biogrid = https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.242/BIOGRID-ALL-4.4.242.mitab.zip  ###################################################################################################################################################
  ## biogrid data type = ppi dataset (protein-protein relationship)
  # setdiff(query_list, tcr_chains)
  
  # #code to get biogrid into right format
  # biogrid_human  <- read_protti("/Users/mgesell/PhD/local_resources/PPIs_and_complexes/BIOGRID-ALL-4.4.242_20250218.mitab.txt") %>%
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
  #   mutate(interactor_a = str_extract(alt_i_ds_interactor_a, "(?<=swiss-prot:)[^|]+(?=\\|refseq)"),
  #          interactor_b = str_extract(alt_i_ds_interactor_b, "(?<=swiss-prot:)[^|]+(?=\\|refseq)")) %>%
  #   select(-c(alt_i_ds_interactor_a, alt_i_ds_interactor_b))
  #  write.csv(biogrid_human, "/Users/mgesell/PhD/local_resources/PPIs_and_complexes/_biogrid_human_mg")
  
  # determine presence of PPI in dataset
  biogrid_phys_query_subset <- read_protti("/Users/mgesell/PhD/local_resources/PPIs_and_complexes/_biogrid_human_mg", head = TRUE) %>%
    dplyr::select(-v1) %>%
    # check which ppis are present - values can range 0-2  (caluclate number query_entries that occur per interactor_a and interactor_b column)
    mutate(ppi = rowSums(across(starts_with("interactor_"), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi == 2) %>% # filter for both interactors identified 
    dplyr::select(-ppi)
  
  # to get for each protein all interactions swap and rowbind interactor_a and interactor_b column  
  biogrid_phys_query_subset <- biogrid_phys_query_subset %>%
    mutate(dummy        = interactor_a,   # store info
           interactor_a = interactor_b,   # overwrite (swap info)
           interactor_b = dummy) %>%      # overwrite (swap info)
    dplyr::select(-dummy) %>%
    rbind(biogrid_phys_query_subset) %>%  # append original matrix
    distinct()                            # not necessary but make sure
  
  # determine node degree
  biogrid_phys_query_subset <- biogrid_phys_query_subset %>% 
    group_by(interactor_a) %>%
    summarize(interactors_biogrid_phys = paste(unique(interactor_b), collapse = ", ")) %>%
    right_join(biogrid_phys_query_subset, by = "interactor_a") %>%
    mutate(node_degree_biogrid =  str_count(interactors_biogrid_phys, ",") + 1)  # each connected protein separated by , --> number of connections = count of , +1 (last has no ,)
  
  # trim df down to unique info
  biogrid_phys_query_subset_trimmed <- biogrid_phys_query_subset %>%
    dplyr::select(interactor_a, interactors_biogrid_phys, node_degree_biogrid) %>%
    distinct() %>%
    arrange(desc(node_degree_biogrid)) %>%
    mutate(comparison = paste(set)) %>%
    rename("node_bg" = "interactor_a")
  
  # investigate query connectivity
  boxplot(biogrid_phys_query_subset_trimmed$node_degree_biogrid)
  summary(biogrid_phys_query_subset_trimmed$node_degree_biogrid)
  paste0(length(query_list), " proteins were queried")                    
  paste0(nrow(biogrid_phys_query_subset_trimmed), " proteins that have PPI according to database")
  paste0(round(nrow(biogrid_phys_query_subset_trimmed)/length(query_list)*100, 1), " % of queried proteins have PPI partner in query") 
  paste0(round(nrow(biogrid_phys_query_subset_trimmed)/length(setdiff(query_list, tcr_chains))*100, 1), " % of queried proteins have PPI partner in query (TCR chains excluded)") 
  paste0(median(biogrid_phys_query_subset_trimmed$node_degree_biogrid), " = median node degree)")
 
   return(biogrid_phys_query_subset_trimmed)
}