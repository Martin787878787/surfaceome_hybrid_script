# belongs to shs_downstream
ppi_string <- function(query_list, mode, set) {
  # # extract and safe experimental string
  # string <- read.csv("~/PhD/local_resources/string/string_net_upsp.csv",  header = TRUE)  %>%
  #   dplyr::select(protein1_entry, protein2_entry, protein1_entry_name, protein2_entry_name, reviewed_protein1, experimental, combined_score) %>%
  #   filter(reviewed_protein1 == "reviewed" & experimental > 0) %>%
  #   dplyr::select(-reviewed_protein1)
  # write.csv(string, "~/PhD/local_resources/string/string_net_physicalMG_upsp.csv")

  # determine presence of PPI in dataset
  string_phys_query_subset <- read_protti("~/PhD/local_resources/string/string_net_physicalMG_upsp.csv", head = TRUE) %>%
    mutate(ppi = rowSums(across(c(protein1_entry, protein2_entry), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi == 2) %>% # filter for both interactors identified 
    dplyr::select(-ppi)
  # boxplot(string_phys_query_subset$experimental) # fyi
  
  # to get for each protein all interactions swap and rowbind interactor_a and interactor_b column  
  string_phys_query_subset <- string_phys_query_subset %>%
    mutate(dummy               = protein1_entry     , # store info
           dummy2              = protein1_entry_name, # store info
           protein1_entry      = protein2_entry     , # overwrite (swap info)
           protein1_entry_name = protein2_entry_name, # overwrite (swap info)
           protein2_entry      = dummy              , # overwrite (swap info)
           protein2_entry_name = dummy2) %>%          # overwrite (swap info)
    dplyr::select(-dummy, -dummy2) %>%
    rbind(string_phys_query_subset) %>%               # append original matrix
    distinct()    # not necessary but make sure
  
  # determine node degree
  string_phys_query_subset <- string_phys_query_subset %>% 
    group_by(protein1_entry) %>%
    summarize(interactors_string_phys = paste(unique(protein2_entry), collapse = ", ")) %>%
    right_join(string_phys_query_subset, by = "protein1_entry") %>%
    mutate(node_degree_string_phys =  str_count(interactors_string_phys, ",") + 1)  # each connected protein separated by , --> number of connections = count of , +1 (last has no ,)
  
  # trim df down to unique info
  string_phys_query_subset_trimmed <- string_phys_query_subset %>%
    dplyr::select(protein1_entry, interactors_string_phys, node_degree_string_phys) %>%
    distinct() %>%
    arrange(desc(node_degree_string_phys)) %>%
    mutate(comparison = paste(set)) %>%
    rename("node_str" = "protein1_entry")
  
  # investigate query connectivity
  boxplot(string_phys_query_subset_trimmed$node_degree_string_phys)
  summary(string_phys_query_subset_trimmed$node_degree_string_phys)
  paste0(length(query_list), " proteins were queried")                    
  paste0(nrow(string_phys_query_subset_trimmed), " proteins that have PPI according to database")
  paste0(round(nrow(string_phys_query_subset_trimmed)/length(query_list)*100, 1), " % of queried proteins have PPI partner in query") 
  paste0(round(nrow(string_phys_query_subset_trimmed)/length(setdiff(query_list, tcr_chains))*100, 1), " % of queried proteins have PPI partner in query (TCR chains excluded)") 
  paste0(median(string_phys_query_subset_trimmed$node_degree_string_phys), " = median node degree)")
  
  #__________________________________________________________________________________________________________________________________________________________________________________________________________________________
  
  return(string_phys_query_subset_trimmed)
  
}