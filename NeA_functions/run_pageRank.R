run_pageRank <- function(dampening_factor, pageRank_percentile_cutoff, pageRank_retain_all_changes){
  ##### load data
  load( file = paste0(result_directory, "intermediate/", "full_network.rds"))
  load( file = paste0(result_directory, "intermediate/", "layout_full.rds" ))
  load( file = paste0(result_directory, "intermediate/", "changes.rds"     ))
  
  #### run page rank
  # Identifies special nodes (changes) and assigns them higher importance in PageRank. 
  seeds             <- rep(0, length(V(network_full)$name))
  index_seed        <- which(V(network_full)$name %in% changes)
  seeds[index_seed] <- 1
  
  pageRank_full     <- page_rank(network_full, algo = "prpack",  personalized = seeds, # Custom jump probabilities (biased toward 'changes' nodes)
            directed = FALSE,           # Treats all edges as bidirectional (common in social/biological networks).
            damping = dampening_factor) # Typically 0.85 (85% chance of following links, 15% chance of jumping).
  
  boxplot_data      <- data.frame(cbind(pageRank_full[[1]], # PageRank scores for ALL nodes
                                        rep("All Proteins", length(V(network_full)$name)))) # Labels
  # Extract scores for the subset (changes)
  pageRank_pbody    <- pageRank_full[[1]][names(pageRank_full[[1]]) %in% changes] 
  temp              <- data.frame( cbind(pageRank_pbody, rep("Assembly Changes Proteins", length(pageRank_pbody))) )
  colnames(temp)    <- colnames(boxplot_data)
  
  boxplot_data      <- rbind(boxplot_data, temp)
  boxplot_data[, 1] <- as.numeric(boxplot_data[, 1])

  p <- ggplot(boxplot_data, aes(x = X2, y = X1, fill = X2)) +            # Applying ggplot function
        geom_boxplot() +
        scale_fill_manual(values = c("grey75", "#0057b7")) +
        labs(y= "PageRank score") +
        ggtitle("Personalized pageRank score distribution") +
        theme_classic() +
        theme(axis.title.x=element_blank(), legend.title = element_blank())
  ggsave(paste0(result_directory, "plots/pageRank_scores.pdf"), p, width = 21, height = 29.7, units = "cm")
  
  # p <- ggplot(boxplot_data, aes(x = X1, fill = X2)) + 
  #   geom_histogram(position = "identity", alpha = 0.7, bins = 100, color = "black") +  # Increased bins for finer x-axis granularity
  #   scale_fill_manual(values = c("grey75", "#0057b7")) + 
  #   scale_y_log10() +  # Logarithmic scale for y-axis
  #   labs(x = "PageRank Score", y = "Number of Proteins (log scale)") + 
  #   ggtitle("Personalized PageRank Score Distribution with Log Scale") + 
  #   theme_classic() + 
  #   theme(legend.title = element_blank())
  # ggsave(paste0(result_directory, "plots/pageRank_scores_hist.pdf"), p, width = 21, height = 29.7, units = "cm")
  
  #### plot full network colored by pageRank score
  V(network_full)$pageRank <-  pageRank_full[[1]]  # Assign PageRank Scores to Nodes
  fine                     <- 500 # this will adjust the resolving power.
  pal                      <- colorRampPalette(c("white", "#0057b7","#001e64")) # color gradient
  
  # Plot Layout Boundaries
  upper   <- 1.1 * apply(layout_full, 2, max)
  lower   <- 1.1 * apply(layout_full, 2, min)
  
  graphCol  <- pal(fine)[as.numeric(cut(V(network_full)$pageRank, breaks = fine))] # Map PageRank Scores to Colors
  
  V(network_full)$size    <- 1200       # Uniform node size
  V(network_full)$color   <- graphCol   # Apply PageRank-based colors
  V(network_full)[which(V(network_full)$changes == TRUE)]$frame.color   <- "#ffc900"   # Highlight "changes" nodes
  # plot
  pdf(paste0(result_directory, "plots/full_network_colored_by_pageRank.pdf"), pointsize = 2)
    plot(network_full, layout = layout_full, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)   # need to adjust limits with new data
    title("Full input network - colored by pageRank score", cex.main = 2)
  dev.off()

  #### get propagated network
  # 1 Identify Significant Nodes via PageRank Cutoff
  pageRank_score_cutoff <- quantile(pageRank_full[[1]], pageRank_percentile_cutoff)  # Step 1: Calculate a cutoff score using a percentile (e.g., top 10% if pageRank_percentile_cutoff = 0.9).
  index_pagerank        <- which(pageRank_full[[1]] >= pageRank_score_cutoff)        # Step 2: Find nodes whose PageRank scores exceed this cutoff.
  
  # 2. Extract Node Names and Indices
  if (pageRank_retain_all_changes == FALSE) {
    propagted_proteins   <- names(pageRank_full[[1]])[index_pagerank]                # Gets the names of significant nodes.
  } else if (pageRank_retain_all_changes == TRUE) {
    propagted_proteins <- union(names(pageRank_full[[1]])[index_pagerank], changes)
  }
  index_pagerank       <- which(V(network_full)$name %in% propagted_proteins)         # Reindexing: Maps these names back to their positions in the original network (network_full).
 
  # 3. Create the Filtered (Propagated) Subnetwork
  network_propagated   <- induced_subgraph(network_full, index_pagerank)
  layout_propagated    <- layout_full[index_pagerank, ]
  
  #### plot propagated network
  upper    <- 1.1 * apply(layout_full, 2, max)
  lower    <- 1.1 * apply(layout_full, 2, min)
  
  pdf(paste0(result_directory, "plots/pageRank_network.pdf"), pointsize = 2)
    plot(network_propagated, layout = layout_propagated, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)   
    title("Network after pageRank", cex.main = 2)
  dev.off()
 
  removed_changes <- unique(setdiff(unique(changes), union(V(network_propagated)$name, missing_in_igraph)))
  print(paste0("..... PageRank report .........................................................................................."))
  print(paste0("PageRank filtering yields   ",  length(V(network_propagated)$name), " nodes & ", gsize(network_propagated), " edges (",   round(length(V(network_propagated)$name)/length(V(network_full)$name)*100,2), " % of orignially ", length(V(network_full)$name), " nodes)"))
  print(paste0("Out of ", length(unique(changes))-length(missing_in_igraph), " proteins from input list (changes-igraphexcluded), ", length(removed_changes), " were excluded from output network due to pageRank quantile filtering"))
  print(paste0("  (Note:   *igraphexcluded* = ",  length(missing_in_igraph), " proteins that were exluded due to 0 edgest (isolated nodes))"))  
  print(paste0("................................................................................................................"))

  write.table(data.frame(name = V(network_propagated)$name), file = paste0(result_directory, pageRank_percentile_cutoff , "_","network_propagated.csv"), col.names = FALSE, row.names = FALSE)
  
  #### export results
  save(network_propagated, file = paste0(result_directory, "intermediate/pageRank_network.rds"))
  save(layout_propagated,  file = paste0(result_directory, "intermediate/layout_pageRank.rds" ))
  save(pageRank_full,      file = paste0(result_directory, "intermediate/pageRank_score.rds"  ))
  write.table(propagted_proteins, paste0(result_directory, "intermediate/propagted_proteins.tsv"), quote = F , sep = "\t", row.names = F)
}
