# belongs to shs_downstream
network_complex_ppi <- function(query_list, complex_db, ppi_db, set) {
  
  # # combine complex databases ------------------------------------------------------------------------------------------------------------------------------------------------------
  # # load corum complexes --
  # ppi_corum <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_human_MG20250305.csv", head=TRUE) %>%
  #   dplyr::select(subunits_uniprot_id) %>%
  #   # get complex format (e.g. Q92828;Q13227;O15379;O75376;O60907;Q9BZK7) to ppi format interactor_a - interactor_b
  #   separate_rows(subunits_uniprot_id, sep = ";") %>%
  #   group_by(complex_id) %>%
  #   mutate(
  #     interactor_a = first(subunits_uniprot_id),
  #     interactor_b = subunits_uniprot_id
  #   ) %>%
  #   filter(interactor_a != interactor_b) %>%
  #   ungroup() %>%
  #   dplyr::select(complex_id, interactor_a, interactor_b) %>%
  #   rename("interaction_identifiers" = "complex_id")
  # 
  # write.csv(ppi_corum, "/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_human_ppi_MG20250306.csv", row.names = FALSE)
  ppi_corum <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_human_ppi_MG20250306.csv", header = TRUE) %>%
    # kick AB BA duplicates
    mutate(AB = paste0(interactor_a, interactor_b),
           BA = paste0(interactor_b, interactor_a)) %>%
    distinct(AB, BA, .keep_all = TRUE) %>%
    dplyr::select(-AB, -BA) 
  
  
  ppi_complex_portal_expanded <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_ppi_MG20250306.csv", header = TRUE)  %>%
    # kick AB BA duplicates
    mutate(AB = paste0(interactor_a, interactor_b),
           BA = paste0(interactor_b, interactor_a)) %>%
    distinct(AB, BA, .keep_all = TRUE) %>%
    dplyr::select(-AB, -BA) 
  
  ## combine corum and complex-portal complexes ====================================================================================================
  ppi_complex_resources <- rbind(ppi_corum, ppi_complex_portal_expanded)
  # fyi ...
  cat(
    "Total complex PPIs:     ", ppi_complex_resources %>% nrow(), "\n",
    "Unique PPIs:            ", 
    ppi_complex_resources %>% dplyr::select(-interaction_identifiers) %>% distinct() %>% nrow(), 
    " (", 
    round((1 - ppi_complex_resources %>% dplyr::select(-interaction_identifiers) %>% distinct() %>% nrow() / 
             ppi_complex_resources %>% nrow()) * 100, 1), 
    "% overlap)\n",
    "   Corum PPIs:          ", ppi_complex_resources %>% filter(grepl("corum_", interaction_identifiers)) %>% distinct() %>% nrow(), "\n",
    "   Complex Portal PPIs: ", ppi_complex_resources %>% filter(grepl("complexportal_", interaction_identifiers)) %>% distinct() %>% nrow(),
    "\n",
    sep = ""   )
  # ================================================================================================================================================
  
  ############################################################################################################################################################################################################################################
  ############################################################################################################################################################################################################################################
  ## combine ppi databases -------------------------------------------------------
  # load string ppi
  ppi_string <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/_string_net_mediumconfidence_physicalMG_upsp20250306.csv", head = TRUE) %>%
    dplyr::select(interaction_identifiers, interactor_a, interactor_b)   %>%
    # kick AB BA duplicates
    mutate(AB = paste0(interactor_a, interactor_b),
           BA = paste0(interactor_b, interactor_a)) %>%
    distinct(AB, BA, .keep_all = TRUE) %>%
    dplyr::select(-AB, -BA) 
  
  # load biogrid ppi
  ppi_biogrid <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/biogrid/_biogrid_human_ppi_global_mg20250306.csv", head = TRUE) %>%
    dplyr::select(interaction_identifiers, interactor_a, interactor_b) %>%
    mutate(interaction_identifiers = gsub(":", "_", interaction_identifiers))   %>%
    # kick AB BA duplicates
    mutate(AB = paste0(interactor_a, interactor_b),
           BA = paste0(interactor_b, interactor_a)) %>%
    distinct(AB, BA, .keep_all = TRUE) %>%
    dplyr::select(-AB, -BA) 
  
  # combine string and biogrid =====================================================================================================================
  ppi_network_resources <- rbind(ppi_string, ppi_biogrid) 
  # fyi ...
  cat(
    "Total network PPIs:      ", ppi_network_resources %>% nrow(), "\n",
    "Unique PPIs:             ", 
    ppi_network_resources %>% dplyr::select(-interaction_identifiers) %>% distinct() %>% nrow(), 
    " (", 
    round((1 - ppi_network_resources %>% dplyr::select(-interaction_identifiers) %>% distinct() %>% nrow() / 
             ppi_network_resources %>% nrow()) * 100, 1), 
    "% overlap)\n",
    "   STRING PPIs:        ", ppi_network_resources %>% filter(grepl("string_", interaction_identifiers)) %>% distinct() %>% nrow(), "\n",
    "   BioGRID PPIs:       ", ppi_network_resources %>% filter(grepl("biogrid_", interaction_identifiers)) %>% distinct() %>% nrow(),
    "\n",
    sep = ""
  )
  # ================================================================================================================================================
  
  ############################################################################################################################################################################################################################################
  ############################################################################################################################################################################################################################################
  ## combine complex and ppi levels -------------------------------------------------------
  ppi_meta <- rbind(ppi_complex_resources, ppi_network_resources)
  cat("corum count:         ", ppi_meta %>% filter(grepl("corum", interaction_identifiers)) %>% distinct() %>% nrow() , "\n",
      "complexportal count: ", ppi_meta %>% filter(grepl("complexportal", interaction_identifiers)) %>% distinct() %>% nrow(), "\n",
      "string count:        ", ppi_meta %>% filter(grepl("string", interaction_identifiers)) %>% distinct() %>% nrow(), "\n",
      "biogrid count:       ", ppi_meta %>% filter(grepl("biogrid", interaction_identifiers)) %>% distinct() %>% nrow(), "\n",
      sep = "")
  # fyi
  barplot(sapply(c("corum", "complexportal",  "biogrid", "string"), function(x) sum(grepl(x, ppi_meta$interaction_identifiers, fixed = TRUE))), 
          main = "PPI Counts by Source",           xlab = "Source",           ylab = "ppi count",           col = c("#fdbb84", "#e34a33", "#9ecae1", "#3182bd")   ) 
  barplot(log10(sapply(c("corum", "complexportal",  "biogrid", "string"), function(x) sum(grepl(x, ppi_meta$interaction_identifiers, fixed = TRUE)))), 
          main = "PPI Counts by Source",           xlab = "Source",           ylab = "log10(ppi count)",    col = c("#fdbb84", "#e34a33", "#9ecae1", "#3182bd")   ) 
  
  # workaround to get wide dataframe (because pivot wider transformation not feasible for dataframe of this size)
  ppi_meta_w <- list(ppi_string  %>% rename("string" = "interaction_identifiers") , 
                     ppi_biogrid %>% rename("biogrid" = "interaction_identifiers") ,
                     ppi_corum   %>% rename("corum" = "interaction_identifiers"), 
                     ppi_complex_portal_expanded %>% rename("complex_portal" = "interaction_identifiers") ) %>%
    reduce(full_join, by = c("interactor_a", "interactor_b")) %>%
    dplyr::select(interactor_a, interactor_b, string, biogrid, corum, complex_portal) # order columns
  
  ############################################################################################################################################################################################################################################
  ############################################################################################################################################################################################################################################
  # testing
  query_list = unique(data_LUX_prot_diff_v31_meta_full %>% pull(entry))
  set = "overall" 
  # testing end
  
  ## extract ppi info for query subset  -------------------------------------------------------
  ppi_meta_qery    <- ppi_meta %>%
    mutate(ppi_query_subset = rowSums(across(starts_with("interactor_"), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi_query_subset == 2) %>% # filter for interactor_a and interactor_b (sum = 2) identified 
    dplyr::select(-ppi_query_subset) 
  
  ppi_meta_w_qery  <- ppi_meta_w %>%
    mutate(ppi_query_subset = rowSums(across(starts_with("interactor_"), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi_query_subset == 2) %>% # filter for interactor_a and interactor_b (sum = 2) identified 
    dplyr::select(-ppi_query_subset) %>%
    mutate()
  
  # total proteins - determine the network   (i) based on these the overall network is constructed
  ppi_condition_meta <- ppi_meta_qery %>% dplyr::select(interactor_a, interactor_b) %>% distinct()
  
  
  comp_cp_LUX
  comp_cor
  
  ppi_str
  ppi_bg
  
  
  
  
  subsets = c("nCD4", "nnCD4", "nCD8", "nnCD8")
  # subset 1    (i) enriched --> color black (ii) recall = 1 --> complex edges green
  # subset 2    (i) enriched --> color black (ii) recall = 1 --> complex edges green
  for (i in data_subsets) {
    # apppend info to w dataframe for this subset
    # ppi detected, corum_recall, complexportal_recall
  }
  
  
  comp_cor 
  
  comp_cp_LUX
  
  
  
  
  
  
  
  
  # check whether complex recall = 1 (or >= 75% ???) - only then set complex_recall to TRUE (== highlight in network)
  
  # set up filter criterium?? e.g. must be in 2/4 databases or such 
  # usefull ???? 
  
  
  
  
}

