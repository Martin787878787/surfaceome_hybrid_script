# belongs to shs_downstream
heat <- function(input_data) {
  
  # ===========================================================================================================================================================================
  # subset and adjust dataframes for downstream plotting
  input_medianlog2_w_matrix <- input_data %>% 
    select(entry, condition, log2_median, overlap)  %>%
    pivot_wider(id_cols = entry, names_from = condition, values_from = log2_median) %>%
    column_to_rownames("entry") %>%
    as.matrix()
  
  input_zscore_w_matrix <- input_data %>% 
    select(entry, condition, log2_z_score, overlap)   
    pivot_wider(id_cols = entry, names_from = condition, values_from = log2_z_score) %>%
      column_to_rownames("entry") %>%
      as.matrix()
  

  
    ################################################################################################################################################
    ## heatmaps & hierarchical clustering  #########################################################################################################
    dist_matrix_medlog2 <- dist(data_medianlog2_w_matrix, method = "euclidean")    # Use binary distance
    hc_medlog2          <- hclust(dist_matrix_medlog2   , method = "ward.D2")   # methods "ward.D2" (best) "average" (ok, looks like complete) or "complete" (ok, looks like average)
    # Cut the tree into a specified number of clusters
    n_clusters_zscore   <- length(unique(data_zscore_L$overlap    )) # Adjust this value based on your needs
    clusters_zscore     <- cutree(hc_zscore , k = n_clusters_zscore)
    # Create a named vector for cluster colors
    cluster_colors_zscore   <- rainbow(n_clusters_zscore)
    names(cluster_colors_zscore)   <- as.character(1:n_clusters_zscore)  # Ensure names match cluster levels
    # Create a row annotation for clusters
    row_annotation_medlog2 <- rowAnnotation(
      Cluster = factor(clusters_medlog2),
      col = list(Cluster = cluster_colors_medlog2)  # Use named vector for colors
    )

    # Create a single heatmap with row annotation
    png(paste0(directory_output, "/_heatmap_data_csc_medianlog2_.png"), width = 1500, height = 15000, res = 300)
    Heatmap(data_csc_medianlog2_w_matrix,
            name = "Expression",
            show_row_names = TRUE,
            show_column_names = TRUE,
            cluster_rows = TRUE,
            column_order = condition_levels_CSC,
            cluster_columns = FALSE,
            row_names_gp = gpar(fontsize = 8),
            col = colors_csc_log2median,
            column_title = "T Cell Types",
            row_title = "Proteins",
            right_annotation = row_annotation_CSC_medlog2 # Use right_annotation for row annotations
    )  
    dev.off() # Close the device
  

    
    # # ===========================================================================================================================================================================
    # # subset and adjust dataframes for downstream plotting
    # data_csc_medianlog2_L <- data_CSC_prot %>% 
    #   select(entry, condition, log2_median, overlap) # %>%    
    # data_csc_zscore_L <- data_CSC_prot %>% 
    #   select(entry, condition, log2_z_score, overlap) # %>%    
    # # long to wide format
    # data_csc_medianlog2_w <- data_csc_medianlog2_L %>%
    #   pivot_wider(id_cols = entry, names_from = condition, values_from = log2_median)
    # data_csc_zscore_w <- data_csc_zscore_L %>%
    #   pivot_wider(id_cols = entry, names_from = condition, values_from = log2_z_score)
    # # wide to matrix
    # data_csc_medianlog2_w_matrix <- data_csc_medianlog2_w %>%
    #   column_to_rownames("entry") %>%
    #   as.matrix()
    # data_csc_zscore_w_matrix <- data_csc_zscore_w %>%
    #   column_to_rownames("entry") %>%
    #   as.matrix()
    
    
    
    
    
  
  
  
  
  input_data
  
  hm <- "asdf"
  
  
  return(hm)
}