# belongs to shs_downstream
complex_compp <- function(query_list, set) {
  
  complex_portal <- read_protti("/Users/mgesell/PhD/local_resources/PPIs_and_complexes/ComplexPortal_human9606_20250218.tsv")   %>%
    filter(taxonomy_identifier == "9606") %>%
    dplyr::select(number_complex_ac, recommended_name, description, identifiers_and_stoichiometry_of_molecules_in_complex, experimental_evidence, complex_assembly) %>%
    mutate(complex_size = str_count(identifiers_and_stoichiometry_of_molecules_in_complex, "\\([0-9]+\\)"),   # each component marked by number which is given in () --> count "(..)" gives number components
           complex_counter_cp = 0, # Initialize complex_counter_cp column in complex_portal with zeros
           matches = "")        # Initialize matches column 
  
  # Iterate through each entry in query_list
  for (entry in query_list) {
    # Use grepl to find matches in complex_portal$identifiers_and_stoichiometry_of_molecules_in_complex
    matches <- grepl(entry, complex_portal$identifiers_and_stoichiometry_of_molecules_in_complex, fixed = TRUE)
    # Increment complex_counter_cp for matching rows
    complex_portal$complex_counter_cp[matches] <- complex_portal$complex_counter_cp[matches] + 1
    # Append the matched entry to the 'matches' column
    complex_portal$matches[matches] <- ifelse(complex_portal$matches[matches] == "", 
                                              entry, 
                                              paste(complex_portal$matches[matches], entry, sep = ", "))
  }
  
  complex_portal_query_subset <- complex_portal %>%
    mutate(subunits_uniprot_ids_cp = gsub("\\([0-9]+\\)", "", identifiers_and_stoichiometry_of_molecules_in_complex),
           subunits_uniprot_ids_cp = gsub("\\|", ", ", subunits_uniprot_ids_cp)) %>%
    dplyr::select(number_complex_ac, description, recommended_name, subunits_uniprot_ids_cp, matches, complex_size, complex_counter_cp) %>%
    filter(complex_counter_cp >= 2) %>%                          # single protein is no complex ...
    mutate(complex_recall_cp = complex_counter_cp/complex_size) %>% # calculate recall
    arrange(desc(complex_recall_cp))  %>%
    mutate(comparison = paste(set),
           subunits_uniprot_ids_cp = gsub(";", ", ", subunits_uniprot_ids_cp),
           subunits_uniprot_ids_cp = sapply(strsplit(subunits_uniprot_ids_cp, ", "), function(x) paste(sort(x), collapse = ", " ))) %>%  # alphabetical display of complex units (easier to compare with matches when ordered)
    rename("complex_name_cp" = "number_complex_ac", "description_cp" = "description")
  
  return(complex_portal_query_subset)
  
}