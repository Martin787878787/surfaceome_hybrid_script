walktrap_clustering <- function(walktrap_step_number){
 
   ##### load data
  load(file = paste0(result_directory, "intermediate/", "pageRank_network.rds"))
  load(file = paste0(result_directory, "intermediate/", "layout_pageRank.rds" ))
  
  ##### run walktrap clustering
  walktrap_clusters <- cluster_walktrap(network_propagated, steps = walktrap_step_number,
                                        merges = FALSE, modularity = TRUE, membership = TRUE)
  
  
  V(network_propagated)$Cluster_membership <- membership(walktrap_clusters)
  
  pdf(paste0(result_directory, "plots/", "Walktrap_cluster_size.pdf"), pointsize = 2)
    hist(table(membership(walktrap_clusters)), breaks = 50)
  dev.off()
  
  ##### export results
  cluster_export_table <- as.matrix(membership(walktrap_clusters))
  write.table(cluster_export_table, paste0(result_directory, "intermediate/", "walktrap_clusters.tsv"), quote = F , sep = "\t", row.names = T)
}