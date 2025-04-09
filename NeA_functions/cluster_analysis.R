cluster_analysis <- function(cluster_method, cluster_min_size, cluster_max_size, cluster_min_changes_above_2_subunits, cluster_min_changes_with_2_subunits,
                             cluster_col_pal = "Dark2", cluster_q_value_cutoff, complex_portal){  # , string_input gene_names,
  ##### load data
  load(file = paste0(result_directory, "intermediate/", "pageRank_network.rds"))
  load(file = paste0(result_directory, "intermediate/", "layout_pageRank.rds" ))
  load(file = paste0(result_directory, "intermediate/", "layout_full.rds"     ))
  load(file = paste0(result_directory, "intermediate/", "changes.rds"         ))
  
  if(cluster_method == "walktrap"){
    cluster <- read.table(paste0(result_directory, "intermediate/","walktrap_clusters.tsv"  ), stringsAsFactors = FALSE)
    V(network_propagated)$Cluster_membership <- cluster$V1 
    
  } else if (cluster_method == "markov"){
    cluster <- read.table(paste0(result_directory, "intermediate/","markov_clusters.tsv"    ), stringsAsFactors = FALSE)
    V(network_propagated)$Cluster_membership <- cluster$V.network_propagated..Cluster_membership 
    
  } else if (cluster_method == "betweenness"){
    cluster <- read.table(paste0(result_directory, "intermediate/","betweenness_clusters.tsv"), stringsAsFactors = FALSE)
    V(network_propagated)$Cluster_membership <- cluster$V1
  }
  
  #### filter clusters based on size
  cluster_sizes    <- table(V(network_propagated)$Cluster_membership)
  index_cluster    <- which(cluster_sizes >= cluster_min_size & cluster_sizes <= cluster_max_size)
  number_of_clusters <- length(index_cluster)
  
  clusters_to_use    <- names(cluster_sizes[index_cluster])
  index_cluster      <- which(V(network_propagated)$Cluster_membership %in% clusters_to_use)
  
  network_clustered  <- induced_subgraph(network_propagated, index_cluster)
  layout_clustered   <- layout_propagated[index_cluster, ]
  
  # count the number of assembly changes in each cluster and condition
  V(network_clustered)$markers_in_cluster <- rep(0, length(V(network_clustered)$name))
  
  for(c in unique(V(network_clustered)$Cluster_membership)){
    nb_markes   <- length(which(V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$changes == TRUE))
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$markers_in_cluster <- nb_markes
  }
  
  # remove the clusters that do not have the required number of marker changes 
  clusters_to_use <- unique(V(network_clustered)$Cluster_membership[which(V(network_clustered)$markers_in_cluster >= cluster_min_changes_above_2_subunits)])
  
  cluster_sizes      <- table(V(network_clustered)$Cluster_membership)
  cluster_2_subunits <- names(cluster_sizes)[which(cluster_sizes == 2)]
  cluster_to_use_two_subunits   <- unique(V(network_clustered)$Cluster_membership[which(V(network_clustered)$markers_in_cluster >= cluster_min_changes_with_2_subunits)])
  cluster_to_use_two_subunits   <- intersect(cluster_2_subunits, cluster_to_use_two_subunits)
  
  clusters_to_use    <- unique(c(clusters_to_use, cluster_to_use_two_subunits))
  
  index_cluster      <- which(V(network_clustered)$Cluster_membership %in% clusters_to_use)
  network_clustered  <- induced_subgraph(network_clustered, index_cluster)
  layout_clustered   <- layout_clustered[index_cluster, ]
  number_of_clusters <- length(clusters_to_use)
  
  # add the name of the complex from the complex portal database
  V(network_clustered)$complex <- rep(NA, length(V(network_clustered)$name))
  df_cluster_name <- data.frame()
  
  for(c in clusters_to_use){
    # get the proteins in the cluster
    proteins_in_complex_propageted <- V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$name
    # get the corresponding complex names
    complexes_propagated <- as.vector(unlist(sapply(proteins_in_complex_propageted, function(x) return(complex_portal[grep(x, complex_portal$`Identifiers (and stoichiometry) of molecules in complex`), 2]))))
    
    # get the proteins with FiLiP marker changes in cluster
    proteins_in_complex_FiLiP <- V(network_clustered)[which(V(network_clustered)$Cluster_membership == c & V(network_clustered)$changes == TRUE)]$name
    # get the corresponding complex names
    complexes_FiliP <- as.vector(unlist(sapply(proteins_in_complex_FiLiP, function(x) return(complex_portal[grep(x, complex_portal$`Identifiers (and stoichiometry) of molecules in complex`), 2]))))
    
    # extract the complexes for which a FiLiP marker was detected
    unique_complexes <- unique(complexes_FiliP)
    if(TRUE %in% grepl("]",unique_complexes)){
      unique_complexes[grepl("]",unique_complexes)] <- strsplit(unique_complexes[grepl("]",unique_complexes)], "\\]")[[1]][2]
    }
    complex_subunits <- as.vector(unlist(sapply(unique_complexes, function(x) return(str_count(as.character(complex_portal[grep(x, complex_portal$`Recommended name`), 5]), "\\|") +1))))
    
    # select the most abundant complex name
    # complex <- names(sort(table(complexes), decreasing = TRUE)[1])
    
    df_complex_summary <- data.frame(Complex = unique_complexes, Subunits = complex_subunits)
    
    # get the number of subunits changing in each complex
    df_complex_summary$detected_subunits_FiLiP <- sapply(df_complex_summary$Complex, function(x) return(length(grep(x, complexes_FiliP))))
    df_complex_summary$detected_subunits_propagated <- sapply(df_complex_summary$Complex, function(x) return(length(grep(x, complexes_propagated))))
    df_complex_summary$fraction_detected <- apply(df_complex_summary, 1, function(x) return(as.numeric(x[3])/as.numeric(x[2])))
    
    df_complex_summary <- df_complex_summary[order(df_complex_summary$detected_subunits_propagated, decreasing = TRUE), ]
    df_complex_summary <- df_complex_summary[order(df_complex_summary$detected_subunits_FiLiP, decreasing = TRUE), ]
    
    #complex <- paste(df_complex_summary[which(df_complex_summary$fraction_detected == max(df_complex_summary$fraction_detected)), "Complex"], collapse = ";")
    #complex <- df_complex_summary[which(df_complex_summary$detected_subunits_FiLiP == max(df_complex_summary$detected_subunits_FiLiP)), "Complex"]
                     
    # if there are more than 2 complexes with the highest number of FiLiP marker changes, chose the two that have the highest number of propagated proteins 
    #if(length(complex)>2){
     #complex <- df_complex_summary[which(df_complex_summary$Complex %in% complex & df_complex_summary$detected_subunits_propagated >= sort(df_complex_summary$detected_subunits_propagated, decreasing = TRUE)[2]), "Complex"]
    #}
    
    # merge into one name
     
    complex <- paste(unlist(apply(df_complex_summary[, c("Complex", "detected_subunits_FiLiP", "detected_subunits_propagated")], 1, function(x) return(paste(x, collapse = "|")))), collapse = ";")
    complex <- strsplit(complex, "\\|")[[1]][1]
    
      
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)[1]]$complex <- complex
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$complex_name <- complex
    
    df_cluster_name <- rbind(df_cluster_name, c(c, as.character(complex)))
  }
  
  df_cluster_name           <- as.data.frame(df_cluster_name)
  colnames(df_cluster_name) <- c("Cluster", "ComplexName")
  write.csv(df_cluster_name, paste0(result_directory, "final/", "cluster_complex_names.csv"), row.names = FALSE, quote = FALSE)
  
  # create a summary df with all the information
  df_summary    <- data.frame(V(network_clustered)$name)
  colnames(df_summary)   <- c("Protein")
  #df_summary$Gene <- sapply(df_summary$Protein, function(x) return(gene_names[grep(x, gene_names$Entry), "Gene.names"]))
  df_summary$Cluster     <- V(network_clustered)$Cluster_membership
  df_summary$Complex     <- V(network_clustered)$complex_name
  
  df_summary$FLiP_Change       <-  V(network_clustered)$FLiP_Change
  df_summary$Fraction_shift    <-  V(network_clustered)$Fraction_shift
  df_summary$homomultimer_link <-  V(network_clustered)$homomultimer_link
 
  #df_summary[which(df_summary$Complex == "PRP19-associated complex"), "Complex"] <- "Splicesome LSM1-7 Pat1"
  write.table(df_summary, paste0(result_directory, "final/", "Complex_Summary.tsv"), quote = FALSE, row.names = FALSE, sep = '\t')
  
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  mycolors <- colorRampPalette(brewer.pal(8, cluster_col_pal))(number_of_clusters)
  
  membership <- unique(V(network_clustered)$Cluster_membership)
  lookup = setNames(mycolors, as.character(membership))
  
  df <- data.frame(as.character(V(network_clustered)$Cluster_membership))
  colnames(df) <- "membership"
  df <- transform(df, membership=lookup[membership], stringsAsFactors=FALSE)
  
  #V(network_clustered)$color <- df[,1]
  V(network_clustered)$color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  V(network_clustered)$size <- 1200
  #V(network_clustered)$frame.color <- df[,1]
  V(network_clustered)$frame.color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  V(network_clustered)[which(V(network_clustered)$FiLiP_changes == TRUE)]$frame.color <- "black"
  
  pdf(paste0(result_directory, "plots/", "clustered_network_size_filter_cluster_label.pdf"), pointsize = 2)
  plot(network_clustered, layout = layout_clustered, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors,vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)   
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  pdf(paste0(result_directory, "plots/", "clustered_network_size_filter.pdf"), pointsize = 2)
  plot(network_clustered, layout = layout_clustered, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  
  #### change coordinates
  layout_clustered_new <- layout_nicely(network_clustered,  )
  
  upper <- 1.1 * apply(layout_clustered_new, 2, max)
  lower <- 1.1 * apply(layout_clustered_new, 2, min)
  
  
  V(network_clustered)$size <- 35
  pdf(paste0(result_directory, "plots/", "clustered_network_size_filter_new_layout.pdf"), pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  pdf(paste0(result_directory, "plots/", "clustered_network_size_filter_new_layout_labeled.pdf"), pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors,vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)    
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  ###### plot LiP anf FLiP changes
  V(network_clustered)$color       <- "grey80"
  V(network_clustered)$frame.color <- "grey80"
  #V(network_clustered)$color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  #V(network_clustered)$frame <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  V(network_clustered)$size <- 40
  
  V(network_clustered)$FLiP_changes <- V(network_clustered)$name %in% changes
  
  V(network_clustered)[which(V(network_clustered)$FLiP_changes == TRUE)]$frame.color      <- "#F7941D" 
  V(network_clustered)[which(V(network_clustered)$Fraction_shift == TRUE)]$color          <- "#1C75BC"
  V(network_clustered)[which(V(network_clustered)$homomultimer_link == TRUE)]$frame.color <- "#EF4136" 
      
  
  #### plot network with clusters
  upper <- 1.1 * apply(layout_clustered_new, 2, max)
  lower <- 1.1 * apply(layout_clustered_new, 2, min)
  
  pdf(paste0(result_directory, "plots/", "clustered_network_significant_with_changes.pdf"), pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)
  title("Changes (yellow), fraction shift (blue), homomultimer link (red)", cex.main = 2)
  dev.off()
  
  pdf(paste0(result_directory, "plots/", "clustered_network_significant_with_changes_labeled.pdf"), pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors, vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5) 
  title("Changes (yellow), fraction shift (blue), homomultimer link (red)", cex.main = 2)
  dev.off()
  
  
  ##### export results
  # write.table(cluster_q, "final/cluster_q_values.tsv", quote = F , sep = "\t", row.names = F)
  save(network_clustered   , file = paste0(result_directory, "final/", "cluster_sig_network.rds"))
  save(layout_clustered_new, file = paste0(result_directory, "final/", "layout_cluster_sig.rds"))
  #save(layout_sig, file = "final/layout_sig_original.rds")
}



