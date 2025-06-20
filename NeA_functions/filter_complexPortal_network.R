filter_complexPortal_network <- function(complexPortal_input_file, complexPortal_mapping_table_file,
                                         changes,
                                         network_output_filename){
  
  ##### load data
  complexPortal_input <- fread(complexPortal_input_file)
  mapping_table       <- fread(complexPortal_mapping_table_file)
  
  ##### transform into an edge list
  list_of_complexes <- unique(complexPortal_input$`Recommended name`)
  
  # get a matrix of all complexes and their uniprot entries as binary interactions
  output <- matrix(ncol=3, nrow=0)
  

  for (i in 1:length(list_of_complexes)){
    current_complex <- list_of_complexes[i]  
    current_line <- subset(complexPortal_input, `Recommended name` == current_complex)
    
    members <- strsplit(current_line$`Identifiers (and stoichiometry) of molecules in complex`, "|", fixed = TRUE)
    members <- list2DF(members)
    colnames(members) <- "proteins"
    members <- gsub("\\s*\\([^\\)]+\\)","",as.character(members$proteins)) # remove brackets and the text within
    
    members <- members[members %in% mapping_table$Entry] # filter for valid uniport IDs
    
    
    if (length(members) > 1){
      
      combinations <- t(combn(members, 2))
      temp <- cbind(rep(current_complex, dim(combinations)[1]), combinations)
      
      output <- rbind(output, temp)
    } else if (length(members) == 1){
      temp <- cbind(current_complex, members, members)
      output <- rbind(output, temp)
    }
    
  }
  
  colnames(output) <- c("Complex", "protein1_entry", "protein2_entry")
  output <- as.data.frame(output)
  
  #### filter for self-loops
  index <- which(output$protein1_entry == output$protein2_entry)
  if (length(index) >0){
    output <- output[-index,]  
  }
  
  
  #### filter for AB = BA
  AB <- paste(output$protein1_entry, "_", output$protein2_entry, sep = "")
  BA <- paste(output$protein2_entry, "_", output$protein1_entry, sep = "")
  
  index <- which(AB == BA)
  if (length(index) >0){
    output <- output[-index,]  
  }
  
  output <- output[, c("protein1_entry", "protein2_entry")] %>% 
    filter(protein1_entry %in% changes | protein2_entry %in% changes) %>%# at least one match with changes list (point is to add +1 interaction layer later)
    distinct()
  
  #### export
  write.table(output, paste0(result_directory, "intermediate/", network_output_filename, ".tsv"), quote = F , sep = "\t", row.names = F)
  
}