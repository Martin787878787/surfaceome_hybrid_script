# belongs to shs_downstream
ppi_biogrid <- function(query_list, mode, set) {
  #################################################################################################################################################################################################################
  ##### biogrid = https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.242/BIOGRID-ALL-4.4.242.mitab.zip  ###################################################################################################################################################
  ## biogrid data type = ppi_query_subset dataset (protein-protein relationship)
  # setdiff(query_list, tcr_chains)
  
  # # #code to get biogrid into right format
  # biogrid_human  <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/biogrid/BIOGRID-ALL-4.4.242_20250218.mitab.txt") %>%
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
  #   dplyr::select(-c(alt_i_ds_interactor_a, alt_i_ds_interactor_b)) %>%
  #   mutate(interaction_detection_method = str_extract(interaction_detection_method, "\\(.*\\)") %>%
  #            str_remove("^\\(") %>%
  #            str_remove("\\)$"))
  #  write.csv(biogrid_human, "/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/biogrid/_biogrid_human_mg", row.names = FALSE)
  biogrid_phys_ <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/biogrid/_biogrid_human_mg", head = TRUE) 
  
  # determine presence of PPI in dataset
  
  # to get for each protein all interactions swap and rowbind interactor_a and interactor_b column  
  biogrid_phys <- biogrid_phys_ %>%
    # kick AB BA duplicates
    mutate(AB = paste0(interactor_a, interactor_b),
           BA = paste0(interactor_b, interactor_a)) %>%
    distinct(AB, BA, .keep_all = TRUE) %>%
    dplyr::select(-AB, -BA) %>%
    # controlled introduction of AB BA
    mutate(dummy        = interactor_a,   # store info
           interactor_a = interactor_b,   # overwrite (swap info)
           interactor_b = dummy) %>%      # overwrite (swap info)
    dplyr::select(-dummy) %>%
    rbind(biogrid_phys_) %>%  # append original matrix
    distinct()                            # not necessary but make sure
  
  # determine all interactors
  biogrid_phys_ppi_global <- biogrid_phys %>% 
    group_by(interactor_a) %>%
    mutate(interactors_biogrid_phys = paste(unique(interactor_b), collapse = ", "),
           node_degree_biogrid_global = length(unique(interactor_b))) %>%
    ungroup() %>%
    # determine all interactors
    # biogrid_phys_ppi_global <- biogrid_phys_ppi_global %>% 
    group_by(interactor_a, interactor_b) %>% # here group per pair the evidence - multiple layers might go for same 
    mutate(evidence_global       = paste(interaction_detection_method, collapse = ", "),
           evidence_count_global = length(interaction_detection_method)) %>%
    ungroup() %>%
    distinct()
  
  # write.csv(biogrid_phys_ppi_global %>% dplyr::select(interaction_identifiers, interactor_a, interactor_b, evidence_global, evidence_count_global),
  #           "/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/biogrid/_biogrid_human_ppi_global_mg20250306.csv", row.names = FALSE)
  
  # same as above just specifically for
  biogrid_phys_query_subset <- biogrid_phys_ppi_global %>%
    mutate(ppi_query_subset = rowSums(across(starts_with("interactor_"), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi_query_subset == 2) %>% # filter for both interactors identified 
    dplyr::select(-ppi_query_subset) %>%
    # -""- for query list
    group_by(interactor_a) %>%
    mutate(interactors_biogrid_query_subset = paste(unique(interactor_b), collapse = ", "),
           node_degree_biogrid_query_subset = length(unique(interactor_b))) %>%
    ungroup() %>%
    group_by(interactor_a, interactor_b) %>% # here group per pair the evidence - multiple layers might go for same 
    mutate(evidence_query_subset       = paste(interaction_detection_method, collapse = ", "),
           evidence_count_query_subset = length(interaction_detection_method)) %>%
    ungroup() %>%
    distinct()
  
  # boxplot node degrees  
  combined_df1 <- bind_rows(
    biogrid_phys_query_subset  %>% dplyr::select(interactor_a, node_degree_biogrid_global)       %>% mutate(source = "query_subset_global") %>% rename("node_degree" = "node_degree_biogrid_global")      ,
    biogrid_phys_query_subset  %>% dplyr::select(interactor_a, node_degree_biogrid_query_subset) %>% mutate(source = "query_subset_intern") %>% rename("node_degree" = "node_degree_biogrid_query_subset"),
    biogrid_phys_ppi_global    %>% dplyr::select(interactor_a, node_degree_biogrid_global)       %>% mutate(source = "global_biogrid")      %>% rename("node_degree" = "node_degree_biogrid_global")
  ) %>% distinct()
  
  # Create the boxplot --------------------------------------------------------------------------------------------------------
  ggplot(combined_df1, aes(x = source, y = node_degree, fill = source)) +
    geom_boxplot() +
    scale_fill_manual(values = c("query_subset_global" = "#74a9cf", 
                                 "query_subset_intern" = "#034e7b", 
                                 "global_string" = "#969696")) +
    stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                 hjust = -0.5, vjust = -1, color = "black", size = 5) +
    labs(x = "Data Source", y = "Node Degree", title = "Boxplot of Node Degrees") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    coord_cartesian(clip = "off")
  #  --------------------------------------------------------------------------------------------------------
  # boxplot evidence count
  combined_df2 <- bind_rows(
    biogrid_phys_query_subset  %>% dplyr::select(interactor_a, evidence_count_global)       %>% mutate(source = "query_subset_global") %>% rename("evidence_count" = "evidence_count_global")      ,
    biogrid_phys_query_subset  %>% dplyr::select(interactor_a, evidence_count_query_subset) %>% mutate(source = "query_subset_intern") %>% rename("evidence_count" = "evidence_count_query_subset"),
    biogrid_phys_ppi_global    %>% dplyr::select(interactor_a, evidence_count_global)       %>% mutate(source = "global_biogrid")      %>% rename("evidence_count" = "evidence_count_global")
  ) %>% distinct()
  
  # Create the boxplot --------------------------------------------------------------------------------------------------------
  ggplot(combined_df2, aes(x = source, y = evidence_count, fill = source)) +
    geom_boxplot() +
    scale_fill_manual(values = c("query_subset_global" = "#74a9cf", 
                                 "query_subset_intern" = "#034e7b", 
                                 "global_string" = "#969696")) +
    stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                 hjust = -0.5, vjust = -1, color = "black", size = 5) +
    labs(x = "Data Source", y = "Node Degree", title = "Number of Evidence Levels (Experimental Evidence)") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    coord_cartesian(clip = "off")
  #  --------------------------------------------------------------------------------------------------------
  
  # trim df down to unique info
  biogrid_phys_query_subset_trimmed <- biogrid_phys_query_subset %>%
    dplyr::select(interactor_a, node_degree_biogrid_query_subset, node_degree_biogrid_global) %>%
    distinct() %>%
    arrange(desc(node_degree_biogrid_query_subset)) %>%
    mutate(comparison = paste(set)) %>%
    rename("node_bg" = "interactor_a")
  
  # investigate query connectivity
  boxplot(biogrid_phys_query_subset_trimmed$node_degree_biogrid_query_subset)
  summary(biogrid_phys_query_subset_trimmed$node_degree_biogrid_query_subset)
  paste0(length(query_list), " proteins were queried")                    
  paste0(nrow(biogrid_phys_query_subset_trimmed), " proteins that have PPI according to database")
  paste0(round(nrow(biogrid_phys_query_subset_trimmed)/length(query_list)*100, 1), " % of queried proteins have PPI partner in query") 
  paste0(round(nrow(biogrid_phys_query_subset_trimmed)/length(setdiff(query_list, tcr_chains))*100, 1), " % of queried proteins have PPI partner in query (TCR chains excluded)") 
  paste0(median(biogrid_phys_query_subset_trimmed$node_degree_biogrid_query_subset), " = median node degree)")
  
  return(biogrid_phys_query_subset_trimmed)
}