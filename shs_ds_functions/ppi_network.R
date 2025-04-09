# belongs to shs_downstream
# provide input as 1 PPI per row
ppi_network <- function(data = ppi_str %>% filter(comparison %in% c("nCD4", "nnCD4", "nCD8", "nnCD8")), 
                        graph_layout = ""
                        ) {
  set.seed(123)  # Any fixed integer - otherwise random (=irreducible) clustering algorithm initialization
  
  ## data preprocessing ____________________________________________________________________________________
  data <- data %>% dplyr::select(protein1_entry_name, protein2_entry_name, comparison)
  
  # initialize df with all unique proteins 
  condition_matrix <- data.frame(protein = unique(c(data$protein1_entry_name, data$protein2_entry_name)))
  # Loop through each condition to create a new column in condition_matrix
  for (conp in unique(data$comparison)) {
    # Filter rows matching the current condition
    condition_data <- subset(data, comparison == conp)
    
    # Identify proteins present in the current condition
    condition_proteins <- unique(c(condition_data$protein1_entry_name, condition_data$protein2_entry_name))
    
    # Create a binary column: 1 if protein matches, 0 otherwise
    condition_matrix[[conp]] <- ifelse(condition_matrix$protein %in% condition_proteins, 1, 0)
  }
  # store protein names to rownames
  condition_matrix <- condition_matrix %>%
    # Create concatenated column (excluding protein)
    mutate(presence_absence = apply(dplyr::select(., -protein), 1, paste, collapse = "-")) %>%
    # Sort by presence_absence in descending order
    arrange(desc(presence_absence)) %>%
    # Move protein to row names (must be last step)
    column_to_rownames("protein")
  
  # _________________________________________________________________________________________________________
  
  # kick out ppi duplications
  data_unique <- data %>% 
    mutate(protein_min = pmin(protein1_entry_name, protein2_entry_name),
           protein_max = pmax(protein1_entry_name, protein2_entry_name)
    ) %>%
    distinct(protein_min, protein_max, .keep_all = TRUE) %>%
    dplyr::select(-protein_min, -protein_max) 
  
  # install.packages("igraph")
  library(igraph)
  graph <- graph_from_data_frame(data_unique,
                             directed = FALSE)
  ## walktrap_result clustering ----------------------------------------------------------------------------------------------------------------
  # Test ideal parameters
  t_results <- sapply(3:10, function(t) {
    wt <- cluster_walktrap(graph, steps = t)
    c(steps = t, 
      modularity = max(wt$modularity),
      clusters = length(unique(membership(wt))))
  })
  # sparse network 3-4; dense network 5-7
  optimal_wt_steps <- t_results[1, which.max(t_results["modularity",])]
  paste0("Applying Walktrep clustering steps:    ", optimal_wt_steps)
  walktrap_result <- cluster_walktrap(graph, 
                         weights = NULL, # If NULL and no such attribute is present, then the edges will have equal weights.
                         steps   = optimal_wt_steps,
                         merges  = TRUE
                         )
  paste0("Walktrap clustering modulairty =    ", round(modularity(walktrap_result),2))
  # Validate cluster quality    >0.3 = meaningfull structure; 0.42 = moderate but non-random clustering
  if(modularity(walktrap_result) < 0.4) {  # modularity ranges -1 to 1. higher = stronger community structure. (calculation by comparison of observed edge-weights within clusters vs. expected weights for random edge distribution)
    warning("Consider alternative algorithms for low modularity networks")
  }
  
  ## Louvian clustering ----------------------------------------------------------------------------------------------------------------
  # Test ideal parameters
  resolutions <- seq(0.1, 1.5, by = 0.1)
  mod_scores <- sapply(resolutions, function(r) {
    max(cluster_louvain(graph, resolution = r)$modularity)
  })
  optimal_louvain_res <- resolutions[which.max(mod_scores)]  # Pick resolution with highest modularity
  paste0("Applying Louvain clustering resolution:    ", optimal_louvain_res)
  louvain_result <- cluster_louvain(graph,
                               resolution = optimal_louvain_res)
  paste0("Louvain clustering modulairty =    ", round(modularity(louvain_result),2))
    # Validate cluster quality    >0.3 = meaningfull structure; 0.42 = moderate but non-random clustering
  if(modularity(louvain_result) < 0.4) {  # modularity ranges -1 to 1. higher = stronger community structure. (calculation by comparison of observed edge-weights within clusters vs. expected weights for random edge distribution)
    warning("Consider alternative algorithms for low modularity networks")
  }
  
  # Leiden clustering ----------------------------------------------------------------------------------------------------------------
  # Test ideal parameters
  resolutions <- seq(0.1, 1.0, by=0.1)
  modularity_scores <- sapply(resolutions, function(r) {
    leiden <- cluster_leiden(graph, resolution_parameter = r)
    max(leiden$modularity)
  })
  optimal_leiden_res <- resolutions[which.max(modularity_scores)]
  paste0("Applying Leiden resolution:    ", optimal_leiden_res)
  leiden_result <- cluster_leiden(graph, resolution_parameter = optimal_leiden_res)
  paste0("Leiden clustering does not calculate modulairty ...")
  
  # cross_tab <- table(Walktrap=membership(walktrap_result), Leiden=membership(leiden_result))
  # pheatmap(log10(cross_tab+10), cluster_rows=F, cluster_cols=F)

  # compare colustering =======
  library(pheatmap)
  cross_tab_wt_lou <- table(Walktrap=membership(walktrap_result), Louvain=membership(louvain_result))
  pheatmap(log10(cross_tab_wt_lou+10), cluster_rows=F, cluster_cols=F)
  cross_tab_wt_lei <- table(Walktrap=membership(walktrap_result), Louvain=membership(leiden_result))
  pheatmap(log10(cross_tab_wt_lei+10), cluster_rows=F, cluster_cols=F)
  
  ## combined clustering methods??
  
  
  ## clustering end ______________________________________________________________________________________________________
  ################################################################################################################################
  
  
  # plotting  ###################################################################################################################
  
  # Define fixed quadrant positions (clockwise from top-left)
  quadrant_order <- c("nCD4", "nnCD4", "nCD8", "nnCD8") # Top-left, Top-right, Bottom-right, Bottom-left

  # Define colors for each condition in quadrant order
  condition_colors <- c("#92c5de", "#f4a582", "#ca0020", "#0571b0")
  
  # Create pie values (always 4 quadrants)
  vertex_pie_values <- lapply(1:nrow(condition_matrix), function(i) rep(1, 4))
  
  # Assign colors based on condition presence (white for 0)
  vertex_pie_colors <- lapply(1:nrow(condition_matrix), function(i) {
    ifelse(condition_matrix[i, quadrant_order] == 1, 
           condition_colors, 
           "white")
  })
  names(vertex_pie_colors) <- rownames(condition_matrix)
  
  # Assign to graph
  V(graph)$pie <- vertex_pie_values
  V(graph)$pie.color <- vertex_pie_colors[V(graph)$name]
  
  
  graph_layout <- "fr" # move this to function parameters
  # determine the layout 
  if (graph_layout == "tree") {
    # 1. Use deterministic tree layout for hierarchical structures
    fixed_layout <- layout_as_tree(graph, root=1)  # Specify root node [6]
  } else if (graph_layout == "fr") {
    # 2. For force-directed layouts, add constraints
    fixed_layout <- layout_with_fr(graph, 
                                   niter=1000,
                                   start.temp=0.1,
                                   grid="nogrid")  # Disable approximate grid [4]
  } 
  # else if (graph_layout == "nicely") {
  #   # 3. Store and reuse layout coordinates
  #   my_layout <- layout_nicely(graph)
  # }
  
  # full network plot
  plot(graph,
       vertex.shape = "pie",
       vertex.pie = V(graph)$pie,
       vertex.pie.color = V(graph)$pie.color,
       vertex.size = 4,
       vertex.label = V(graph)$name,
       edge.width = 1.5,
       layout = fixed_layout,
       vertex.pie.angle = pi/4) # Rotate 45° to align quadrants
  
  
  # install.packages("ggraph")
  # install.packages("tidygraph")
  # library(ggraph)
  # library(tidygraph)
  # graph <- as_tbl_graph(data_unique)  
  # 
  # ggraph(graph, layout = "fr") +  # Fruchterman-Reingold (force-directed)
  #   geom_edge_link(alpha = 0.5) +  
  #   geom_node_point() +
  #   theme_void()
  
  
  
  
  
  
  
  return(data_ordered) 
}