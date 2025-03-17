# belongs to shs_downstream
complex_corum <- function(query_list, set) {
  # 
  # not thate there are no spaces behind the ":"
  # corum_accepted_evidence <- c("MI:0004:affinity chromatography technology", "MI:0006:anti bait coimmunoprecipitation", "MI:0007:anti tag coimmunoprecipitation",
  #                              "MI:0018:two hybrid", "MI:0019:coimmunoprecipitation", "MI:0030:cross-linking study", "MI:0031:protein cross-linking with a bifunctional reagent",
  #                              "MI:0047:far western blotting", "MI:0055:fluorescent resonance energy transfer", "MI:0059:gst pull down", "MI:0096:pull down",
  #                              "MI:0107:surface plasmon resonance", "MI:0114:x-ray crystallography", "MI:0676:tandem affinity purification", "MI:0437:protein three hybrid"
  # )
  # # to display unique evidence levels execute this:
  #  # corum$purification_methods %>%
  #  #       strsplit(";") %>%
  #  #       unlist() %>%
  #  #       unique()
  # # load, filter and subset corum
  # corum <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/corum_humanComplexes_2025-01-21.txt") %>% # colnames(corum)
  #   filter(organism == "Human") %>% # not necessary
  #   mutate(complex_id = paste0("corum_", complex_id)) %>%
  #   dplyr::select(complex_id, pmid, complex_name, synonyms, subunits_uniprot_id, comment_complex, comment_members, comment_disease, comment_drug, purification_methods, functions_go_name) %>%
  #   filter(str_detect(purification_methods, paste(corum_accepted_evidence, collapse = "|"))) %>%
  #   mutate(complex_size = str_count(subunits_uniprot_id, ";")+1 ) %>%  # each component speparated by ; --> count occurences + 1
  #   mutate(complex_counter_cor = 0, # Initialize complex_counter_cor column in complex_portal with zeros
  #          matches = "")        # Initialize matches column
  # 
  # # get rid of isoform info
  # corum$subunits_uniprot_id <- gsub("-[^;]+(?=;|$)", "", corum$subunits_uniprot_id, perl = TRUE)
  # write.csv(corum, "/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_human_MG20250305.csv", row.names = FALSE)
  corum <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_human_MG20250305.csv", head = TRUE) 
  
  
  # Iterate through each entry in df2$entry
  for (entry in query_list) {
    # Use grepl to find matches in corum$columnX
    matches <- grepl(entry, corum$subunits_uniprot_id, fixed = TRUE)
    # Increment complex_counter_cor for matching rows
    corum$complex_counter_cor[matches] <- corum$complex_counter_cor[matches] + 1
    # Append the matched entry to the 'matches' column
    corum$matches[matches] <- ifelse(is.na(corum$matches[matches]) | corum$matches[matches] == "",
                                     entry,
                                     paste(corum$matches[matches], entry, sep = ", "))
  }
  
  corum_query_subset <- corum %>%
    dplyr::select(complex_id, complex_name, comment_complex, subunits_uniprot_id, matches, complex_size, complex_counter_cor) %>%
    filter(complex_counter_cor >= 2) %>%                          # single protein is no complex ...
    mutate(complex_recall_cor = complex_counter_cor/complex_size) %>% # calculate recall
    filter(complex_recall_cor  >= 0.5) %>%
    arrange(desc(complex_recall_cor)) %>%
    mutate(comparison = paste(set),
           subunits_uniprot_id = gsub(";", ", ", subunits_uniprot_id),
           subunits_uniprot_id = sapply(strsplit(subunits_uniprot_id, ", "), function(x) paste(sort(x), collapse = ", " ))) %>%  # alphabetical display of complex units (easier to compare with matches when ordered)) %>%
    rename("complex_name_cor" = "complex_name", "description_cor" = "comment_complex", "subunits_uniprot_id_cor" = "subunits_uniprot_id")
  
  #__________________________________________________________________________________________________________________________________________________________________________________________________________________________
  
  return(corum_query_subset)
  
}





