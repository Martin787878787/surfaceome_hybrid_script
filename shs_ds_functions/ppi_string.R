# belongs to shs_downstream
ppi_string <- function(query_list, mode, set, string_phys_confidence) {
  # # extract and safe experimental string
  # string <- read.csv("~/PhD/local_resources/string/string_net_upsp.csv",  header = TRUE)  %>%
  #   dplyr::select(protein1_entry, protein2_entry, protein1_entry_name, protein2_entry_name, reviewed_protein1, experimental, combined_score) %>%
  #   filter(reviewed_protein1 == "reviewed" & experimental > 0) %>%
  #   dplyr::select(-reviewed_protein1)
  # write.csv(string, "~/PhD/local_resources/string/string_net_physicalMG_upsp.csv")
  
  # load network - see above for how the table was created (TABLE IS UNFILTERED - CONFIDENCE FILTERING HAPPENS BELOW)
  string_phys_ <- read_protti("~/PhD/local_resources/string/string_net_physicalMG_upsp.csv", head = TRUE) %>%
    dplyr::select(-v1)
  
  # define the stringency for physical interaction
  # not that this is a guesstimate of Martin based STRING recommendation to consider low / medium / high confidence (<0.4, >0.4, >0.7 combined_score)
  # summary(string_phys_$experimental)
    # Min.    1st Qu.     Median         Mean     3rd Qu.        Max. 
    # 42.0    >>71.0<<    >>100.0<<      152.1    >>166.0<<      999.0 
  # ggplot(data = string_phys_, aes(x = "", y = experimental)) +
  #   geom_boxplot() +
  #   stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
  #                vjust = -0.5)
  if (string_phys_confidence        == "low"   ) {        # string_phys_confidence = "low"
    string_phys_cutoff <- summary(string_phys_$experimental)["1st Qu."]   # >=  71.0
  } else if (string_phys_confidence == "medium") {        # string_phys_confidence = "medium"
    string_phys_cutoff <- summary(string_phys_$experimental)["Median"]    # >= 100.0
  } else if (string_phys_confidence == "high"  ) {        # string_phys_confidence = "high"
    string_phys_cutoff <- summary(string_phys_$experimental)["3rd Qu."]   # >= 166.0
  } else if (string_phys_confidence == "super_high"  ) {  # string_phys_confidence = "super_high"
    string_phys_cutoff <- 310   # >= 166.0
  }
  
 # appy string confidence filtering
 string_phys <- string_phys_ %>%
    filter(experimental >= string_phys_cutoff)
 paste0("Keeping ", round(nrow(string_phys)/nrow(string_phys_)*100,3), " % of string PPIs")

  # ensure to include A-B , B-A ppi relationship in string match frame (so no quer_list protein slips trough)
  string_phys <- string_phys %>%
    mutate(dummy               = protein1_entry     , # store info
           dummy2              = protein1_entry_name, # store info
           protein1_entry      = protein2_entry     , # overwrite (swap info)
           protein1_entry_name = protein2_entry_name, # overwrite (swap info)
           protein2_entry      = dummy              , # overwrite (swap info)
           protein2_entry_name = dummy2) %>%          # overwrite (swap info)
      dplyr::select(-dummy, -dummy2) %>%
      rbind(string_phys) %>%             # append original matrix
      distinct()   # make sure no duplications created
     
  # determine the phyiscal node degree of all string nodes
  string_phys_degree_global <- string_phys %>%
      group_by(protein1_entry) %>%
      summarize(interactors_string_phys_global = paste(unique(protein2_entry), collapse = ", "))  %>%
      mutate(node_degree_string_phys_global =  str_count(interactors_string_phys_global, ",") + 1)  # each connected protein separated by , --> number of connections = count of , +1 (last has no ,)
    
  # filter string for PPIs contained in query_list 
  string_phys_query_subset <- string_phys %>%
    mutate(ppi = rowSums(across(c(protein1_entry, protein2_entry), ~.x %in% query_list))) %>% #   0 neither, 1 one, 2 both interactors identified
    filter(ppi == 2) %>% # filter for both interactors identified 
    dplyr::select(-ppi) 
  # determine the phyiscal node degree WITHIN query list 
  node_degree_string_phys_query_subset <- string_phys_query_subset %>%
    group_by(protein1_entry) %>%
    summarize(interactors_string_phys_query_subset = paste(unique(protein2_entry), collapse = ", "))  %>%
    mutate(node_degree_string_phys_query_subset    =  str_count(interactors_string_phys_query_subset, ",") + 1)  # each connected protein separated by , --> number of connections = count of , +1 (last has no ,)
  
  # join 
  string_phys_query_subset <- string_phys_query_subset %>%
    left_join(string_phys_degree_global, by = "protein1_entry") %>%
    left_join(node_degree_string_phys_query_subset, by = "protein1_entry") %>%
    distinct()
  
  # boxplot node degrees  
  combined_df <- bind_rows(
    string_phys_query_subset  %>% dplyr::select(protein1_entry, node_degree_string_phys_global)       %>% mutate(source = "query_subset_global") %>% rename("node_degree" = "node_degree_string_phys_global")      ,
    string_phys_query_subset  %>% dplyr::select(protein1_entry, node_degree_string_phys_query_subset) %>% mutate(source = "query_subset_intern") %>% rename("node_degree" = "node_degree_string_phys_query_subset"),
    string_phys_degree_global %>% dplyr::select(protein1_entry, node_degree_string_phys_global)       %>% mutate(source = "global_string")       %>% rename("node_degree" = "node_degree_string_phys_global")
  ) %>% distinct()
  
  # Create the boxplot --------------------------------------------------------------------------------------------------------
  ggplot(combined_df, aes(x = source, y = node_degree, fill = source)) +
    geom_boxplot() +
    scale_fill_manual(values = c("query_subset_global" = "#74a9cf", 
                                 "query_subset_intern" = "#034e7b", 
                                 "global_string" = "#969696")) +
    stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                 hjust = -0.5, vjust = -1, color = "black", size = 5) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1, "cm")
    ) +
    labs(x = "Data Source", y = "Node Degree", title = "Boxplot of Node Degrees") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(clip = "off")

    # geom_boxplot() +
    # scale_fill_manual(values = c("query_subset_global" = "#74a9cf", 
    #                              "query_subset_intern" = "#034e7b", 
    #                              "global_string" = "#969696")) +
    # stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
    #              vjust = -0.5, color = "black", size = 5) +
    # theme_minimal() +
    # theme(
    #   axis.title.x = element_text(size = 14),
    #   axis.title.y = element_text(size = 14),
    #   axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    #   axis.text.y = element_text(size = 12),
    #   plot.title = element_text(size = 16),
    #   legend.title = element_text(size = 14),
    #   legend.text = element_text(size = 12),
    #   legend.key.size = unit(1, "cm")
    # ) +
    # labs(x = "Data Source", y = "Node Degree", title = "Boxplot of Node Degrees") +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #  --------------------------------------------------------------------------------------------------------
  
  # trim df down to unique info
  string_phys_query_subset_trimmed <- string_phys_query_subset %>%
    dplyr::select(protein1_entry, interactors_string_phys_global, node_degree_string_phys_global, interactors_string_phys_query_subset, node_degree_string_phys_query_subset) %>%
    distinct() %>%
    arrange(desc(node_degree_string_phys_global)) %>%
    mutate(comparison = paste(set)) %>%
    rename("node_str" = "protein1_entry")
  
  # investigate query connectivity
  paste0(length(unique(query_list)), " proteins were queried")                    
  paste0(length(unique(string_phys_query_subset_trimmed$node_str)), " proteins that have PPI according to database")
  paste0(round(length(unique(string_phys_query_subset_trimmed$node_str))/length(query_list)*100, 1), " % of queried proteins have PPI partner in query") 
  paste0(round(length(unique(string_phys_query_subset_trimmed$node_str))/length(setdiff(query_list, tcr_chains))*100, 1), " % of queried proteins have PPI partner in query (TCR chains excluded)") 
  paste0(median(string_phys_query_subset_trimmed$node_degree_string_phys_global),       " = median global node degree)")
  paste0(median(string_phys_query_subset_trimmed$node_degree_string_phys_query_subset), " = median query intrinsic node degree)")
  
  #__________________________________________________________________________________________________________________________________________________________________________________________________________________________
  
  return(string_phys_query_subset_trimmed)
  
}
