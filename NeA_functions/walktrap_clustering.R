walktrap_clustering <- function(walktrap_step_number){
 
   ##### load data
  load(file = paste0(result_directory, "intermediate/", "pageRank_network.rds"))
  load(file = paste0(result_directory, "intermediate/", "layout_pageRank.rds" ))
  
  
  ##### run walktrap clustering based on specific input value defining step number or based on automated parameter optimization ("auto")
  if (!is.numeric(walktrap_step_number) && walktrap_step_number != "auto") {
    stop("walktrap_step_number must be numeric or 'auto'")
  } 
  if (is.numeric(walktrap_step_number)) {
    walktrap_clusters <- cluster_walktrap(network_propagated, 
                                          steps      = walktrap_step_number,
                                          merges     = FALSE, 
                                          modularity = TRUE, 
                                          membership = TRUE)
  } else if (walktrap_step_number == "auto") {
    # Test ideal parameters
    t_results <- sapply(3:10, function(t) {
      wt <- cluster_walktrap(network_propagated, 
                             steps = t, 
                             merges = FALSE, 
                             modularity = TRUE, 
                             membership = TRUE)
      c(steps = t, 
        modularity = max(wt$modularity),
        clusters = length(unique(membership(wt))))
    })
    # sparse network 3-4; dense network 5-7
    optimal_wt_steps <- t_results[1, which.max(t_results["modularity",])]
    message("Applying Walktrep clustering steps:    ", optimal_wt_steps)
    walktrap_clusters <- cluster_walktrap(network_propagated, 
                                        steps      = optimal_wt_steps,
                                        merges     = FALSE, 
                                        modularity = TRUE, 
                                        membership = TRUE
                                        
    )
    paste0("Walktrap clustering modulairty =    ", round(modularity(walktrap_clusters),2))
    # Validate cluster quality    >0.3 = meaningfull structure; 0.42 = moderate but non-random clustering
    if(modularity(walktrap_clusters) < 0.4) {  # modularity ranges -1 to 1. higher = stronger community structure. (calculation by comparison of observed edge-weights within clusters vs. expected weights for random edge distribution)
      warning("Consider alternative algorithms for low modularity networks")
    }
    
  } 


  V(network_propagated)$Cluster_membership <- membership(walktrap_clusters)
  
  pdf(paste0(result_directory, "plots/", "Walktrap_cluster_size.pdf"), pointsize = 2)
    hist(table(membership(walktrap_clusters)), breaks = 50)
  dev.off()
  
  ##### export results
  cluster_export_table <- as.matrix(membership(walktrap_clusters))
  write.table(cluster_export_table, paste0(result_directory, "intermediate/", "walktrap_clusters.tsv"), quote = F , sep = "\t", row.names = T)
}