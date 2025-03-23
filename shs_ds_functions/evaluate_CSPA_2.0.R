# belongs to shs_downstream
evaluate_CSPA_2.0 <- function(cell_type  ) {  #  cell_type = "Jurkat"   /   cell_type = "panT"  /  cell_type = ...

  #  CSPA_1.0             = read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_per_cell_type.csv"    , header = TRUE, sep = ",")
  cell_specific_CSPA   = read.csv(paste0("/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_2.0/", cell_type, "_CSPA.csv"), header = TRUE, sep = ",")
  
  # Upset Plot CSPA 1.0 vs. CSPA 1.0_cell_type vs. CSPA 2.0_cell_type
  colnames(cell_specific_CSPA[,c(5:7)])
  
  cell_specific_CSPA %>%
    dplyr::select(1, 5, 6, 7) %>%
    mutate(across(-1, as.numeric)) %>%
    column_to_rownames(var = names(.)[1]) %>%
    upset(nsets = 3, 
          mainbar.y.label = "Unique IDs Count",
          sets.x.label = "Set Size",
          text.scale = 1.2)

  if( sum(cell_specific_CSPA[, 6]) > 0 ){   # rescues cell lines not in CSPA_1.0 
    cell_specific_CSPA %>%
      dplyr::select(1, 5, 6) %>%
      mutate(across(-1, as.numeric)) %>%
      column_to_rownames(var = names(.)[1]) %>%
      upset(nsets = 3, 
            mainbar.y.label = "Unique IDs Count",
            sets.x.label = "Set Size",
            text.scale = 1.2)
  }
  
  cell_specific_CSPA %>%
    dplyr::select(1, 4 , 5) %>%
    mutate(across(-1, as.numeric)) %>%
    column_to_rownames(var = names(.)[1]) %>%
    upset(nsets = 3, 
          mainbar.y.label = "Unique IDs Count",
          sets.x.label = "Set Size",
          text.scale = 1.2)

}