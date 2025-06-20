# purpose:   e.g. for extraction of condition (t cell subset) specific ppi network
extract_changes_internal_network <- function(changes, network_output_filename) {
  # read network (this one contains +1 interactor layer)
  network <- read.table(paste(result_directory, "intermediate/", network_output_filename, ".tsv", sep = "" ), sep = "\t", header =  TRUE)
  
  # filter for ppis that could be formed by changes list
  sub_network <- network %>%
    filter(protein1_entry %in% changes & protein2_entry %in% changes) 
  
  # export
  write.table(sub_network, paste(result_directory, "intermediate/", network_output_filename, "_changes_sub-network.tsv", sep = "" ), quote = F , sep = "\t", row.names = F)

}