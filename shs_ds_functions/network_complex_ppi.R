# belongs to shs_downstream
network_complex_ppi <- function(query_list, complex_db, ppi_db, set) {
  
  # # combine complex databases ------------------------------------------------------------------------------------------------------------------------------------------------------
  # # load corum complexes --
  # ppi_corum <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_human_MG20250305.csv", head=TRUE) %>%
  #   dplyr::select(subunits_uniprot_id) %>%
  #   mutate(complex_id = paste0("corum_", row_number()))  %>%  # Add a unique ID for each complex
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
  
  ## load biogrid complexes --------------------------------------------------------------------------------------------------------------------------------------------------------
  # !!! CAVEAT i assume [] indicates "one or the other" in complex. so far did not account for this and just considered as ppis also withing bracket ... !!!
  # complex_portal_ <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_human_MG20250305.csv", head=TRUE)  %>%
  #   mutate(subunits_uniprot_ids_cp = gsub("\\([0-9]+\\)", ""     , identifiers_and_stoichiometry_of_molecules_in_complex),
  #          subunits_uniprot_ids_cp = gsub("\\[|\\]"         , ""     , subunits_uniprot_ids_cp),  # destroy sub-complex structure
  #          subunits_uniprot_ids_cp = gsub("\\|" , ","    , subunits_uniprot_ids_cp),  # destroy sub-complex structure
  #          subunits_uniprot_ids_cp = gsub("-PRO_"       , ",PRO_",  subunits_uniprot_ids_cp)  # replace - (which i think indicate protein RNA complex) with , to indicate complex members
  #   )   %>%
  #   dplyr::select(number_complex_ac, subunits_uniprot_ids_cp, complex_size)
  # 
  # # some complexes are described as complexes of complexes - in order to resolve this relationship to protein only annotation ...
  # # Function to expand complex identifiers
  # expand_complex <- function(complex_id, df) {
  #   # Locate the row where number_complex_ac matches the complex_id
  #   row <- df %>% filter(number_complex_ac == complex_id)
  #   # If the complex_id is not found, return the original complex_id
  #   if (nrow(row) == 0) {
  #     return(complex_id)
  #   }
  #   # Return the protein identifiers for that complex
  #   return(row$subunits_uniprot_ids_cp)
  # }
  # # Recursive replacement function
  # replace_complexes <- function(data) {
  #   new_data <- data %>%
  #     mutate(subunits_uniprot_ids_cp = sapply(subunits_uniprot_ids_cp, function(ids) {
  #       while (str_detect(ids, "CPX-")) {
  #         # Extract the first complex identifier found
  #         complex_id <- str_extract(ids, "CPX-\\d+")
  #         # Check if a complex_id was actually found
  #         if (is.na(complex_id)) {
  #           break # Exit the while loop if no complex ID is found
  #         }
  #         # Expand the complex identifier
  #         expanded_ids <- expand_complex(complex_id, data)
  #         # Replace the complex identifier with the expanded IDs
  #         ids <- str_replace(ids, complex_id, expanded_ids)
  #       }
  #       return(ids)
  #     }))
  #   return(new_data)
  # }
  # 
  # # Apply the recursive replacement
  # complex_portal_expanded <- replace_complexes(complex_portal_) %>% # Add a unique ID for each complex
  #   # get complex format (e.g. Q92828,Q13227,O15379,O75376,O60907,Q9BZK7) to ppi format interactor_a - interactor_b
  #   separate_rows(subunits_uniprot_ids_cp, sep = ",")
  # 
  #  # get pro mapping (proteoforms are saved in complex portal - in order to match them with dataset translate to uniprot reviewed format)
  #  human_unreviewd_isoforms_mapping <- read_protti("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/resources/resource_data_uniprot/human_UNREVIEWED_isoform_mapping_20250306.tsv", header = TRUE) %>%
  #    # Extract all PRO_ values and collapse them into a single string per row
  #    mutate(across(c("chain", "propeptide", "initiator_methionine",
  #                    "post_translational_modification", "lipidation",
  #                    "cross_link", "disulfide_bond", "glycosylation",
  #                    "modified_residue", "peptide", "transit_peptide",
  #                    "signal_peptide"), # specify your columns here
  #                  ~ str_extract_all(., "PRO_\\d+", simplify = TRUE) %>%
  #                    apply(1, paste, collapse = " ") %>%
  #                    str_trim()
  #    )) %>%
  #    mutate(all_pro_values = paste(chain, propeptide, initiator_methionine,
  #                                  post_translational_modification, lipidation,
  #                                  cross_link, disulfide_bond, glycosylation,
  #                                  modified_residue, peptide, transit_peptide,
  #                                  signal_peptide, sep = " ")) %>%
  # 
  #    dplyr::select(entry, entry_name, reviewed, all_pro_values) %>%
  #    mutate(all_pro_values = str_squish(all_pro_values)) %>%  # clean up unnecessary spaces (e.g., trailing spaces or double spaces)
  #    filter(all_pro_values != "") %>%                         # remove empty rows
  #    separate_longer_delim(all_pro_values, delim = " ") # get each PRO_ in to separate colum (sometimes multiple isoforms per protein --> to match them all assign individually (row wise) to entry)
  # 
  #  # translate PRO isoforms (uniprot>ptm>chain ) to uniprot format
  #  complex_portal_expanded_pro <- complex_portal_expanded %>% filter(grepl("PRO_", subunits_uniprot_ids_cp)) %>%
  #    left_join(human_unreviewd_isoforms_mapping, by = c("subunits_uniprot_ids_cp" = "all_pro_values") )
  #  # after transaltion trim to original format for merging
  #  table(complex_portal_expanded_pro$reviewed)
  #  # after translation trim to original format for merging
  #  complex_portal_expanded_pro <- complex_portal_expanded_pro %>%
  #    mutate(subunits_uniprot_ids_cp = entry) %>%
  #    dplyr::select(number_complex_ac, subunits_uniprot_ids_cp, complex_size)
  # 
  #  # get data in ppi format
  #  ppi_complex_portal_expanded <- rbind(complex_portal_expanded_pro, complex_portal_expanded %>% filter(!grepl("PRO_", subunits_uniprot_ids_cp)) ) %>%
  #    # kick isoforms ASDF123-1 ASDF123-3
  #    mutate(subunits_uniprot_ids_cp = str_remove(subunits_uniprot_ids_cp, "-.*$")) %>%
  #    distinct()  %>% # ensure unique values (might be residue from isoforms)
  #    group_by(number_complex_ac) %>%
  #    rename("interaction_identifiers" = "number_complex_ac", "subunits_uniprot_id" = "subunits_uniprot_ids_cp") %>%
  #    mutate(
  #      interactor_a = first(subunits_uniprot_id),
  #      interactor_b = subunits_uniprot_id       ,
  #      interaction_identifiers = gsub("CPX-", "complexportal_", interaction_identifiers)) %>%
  #    filter(interactor_a != interactor_b) %>%
  #    ungroup() %>%
  #    dplyr::select(interaction_identifiers, interactor_a, interactor_b)
  # #
  #  write.csv(ppi_complex_portal_expanded, "/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_ppi_MG20250306.csv", row.names = FALSE)
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
  summary(ppi_meta$interaction_identifiers)
  cat("corum count:         ", sum(ppi_meta$interaction_identifiers == "corum"), "\n",
      "string count:        ", sum(ppi_meta$interaction_identifiers == "string"), "\n",
      "complexportal count: ", sum(ppi_meta$interaction_identifiers == "complexportal"), "\n",
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
  ## extract ppi info for query subset  -------------------------------------------------------
  ppi_meta_qery    <- ppi_meta %>%
    mutate(ppi_query_subset = rowSums(across(starts_with("interactor_"), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi_query_subset == 2) %>% # filter for interactor_a and interactor_b (sum = 2) identified 
    dplyr::select(-ppi_query_subset) 
  
  ppi_meta_w_qery  <- ppi_meta_w %>%
    mutate(ppi_query_subset = rowSums(across(starts_with("interactor_"), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi_query_subset == 2) %>% # filter for interactor_a and interactor_b (sum = 2) identified 
    dplyr::select(-ppi_query_subset) 
  
  
  # check whether complex recall = 1 (or >= 75% ???) - only then set complex_recall to TRUE (== highlight in network)
  
  # set up filter criterium?? e.g. must be in 2/4 databases or such 
  # usefull ???? 
  
  
  
  
}

