create_igraph_objects <- function(add_non_connected = FALSE,
                                  network_edgelist_file, changes,
                                  remove_loop_edges = TRUE){

  ##### load data
  network_edgelist <- as.matrix(fread(network_edgelist_file))
  
  #### create igraph object
  network_full <- graph_from_edgelist(network_edgelist, directed = FALSE)
  
  #### simplify
  network_full <- simplify(network_full, remove.multiple = TRUE, remove.loops = remove_loop_edges)
  
  
  #### add not-connected vertices from p-body data
  # if (add_non_connected){
  #   index_pbody <- which(FiLiP_changes %in% V(network_full)$name)
  #   missing_v <- FiLiP_changes[-index_pbody]
  #   
  #   network_full <- add_vertices(network_full, length(missing_v), name = missing_v)
  #   
  # }
  # 
  
  #### plot the network as force directed network
  # layout_full <- layout_with_drl(network_full,  )
  layout_full <- layout_with_graphopt(network_full,  )
  
  
  #### define network plot properties
  # V(network_full)$size <- 150
  V(network_full)$size        <- 600       # dot size
  V(network_full)$color       <- "grey80"
  V(network_full)$frame.color <- "grey80"
  V(network_full)$label       <- ""
  E(network_full)$color       <- "grey90"
  
  
  #### get limits for plot
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  pdf(paste0(result_directory,"plots/full_input_network.pdf"), pointsize = 2)
    plot(network_full, layout = layout_full, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)   # need to adjust limits with new data
    title("Full input network", cex.main = 2)
  dev.off()
  
  
  ##### add info change

  ##### add info from FLiP marker changes
  index_changes                              <- which(V(network_full)$name %in% changes)
  V(network_full)$changes                    <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$changes[index_changes]     <- TRUE
  V(network_full)$size[index_changes]        <- 1200
  V(network_full)$frame.color[index_changes] <- "#2680FF"
  V(network_full)$size[index_changes]        <- 1200
  
  pdf(paste0(result_directory, "plots/full_input_network_FiLiP_changes.pdf"), pointsize = 2)
    plot(network_full, layout = layout_full, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)
    title("Full input network - FiLiP Changes highlighted", cex.main = 2)
  dev.off()
  
  # igraph_diagnostics ##############################################################################################
  missing_in_igraph = setdiff(changes, V(network_full)$name)
  present_in_igraph = intersect(changes, V(network_full)$name)
  unconnected_changes_nodes = present_in_igraph[degree(network_full)[present_in_igraph] == 0]
  # Create diagnostic output
  igraph_diagnostics <- list(
    missing_in_igraph = missing_in_igraph,
    unconnected_changes_nodes = unconnected_changes_nodes,
    present_in_igraph = present_in_igraph
  )
  print(paste0("=== Protein Loss Report ==="))
  # Section 1: Proteins completely missing from network_full
  if (length(missing_in_igraph) > 0) {
    print(paste0("[Critical] Proteins missing from network_full (", length(missing_in_igraph), "):"))
    print(missing_in_igraph)
    cat("These proteins were not found in the edgelist and weren't added to the network_full.")
  } else {
    print(paste0("[OK] All changes proteins exist in the network_full.\n"))
  }
  
  # Section 2: Proteins with no connections
  if (length(unconnected_changes_nodes) > 0) {
    print(paste0("[Warning] Isolated changes proteins (", length(unconnected_changes_nodes), "):"))
    print(unconnected_changes_nodes)
    print(paste0("These exist in the network_full but have no interactions (degree = 0)."))
  } else {
    print(paste0("[OK] No isolated changes proteins found."))
  }
  # Section 3: Summary statistics
  print(paste0("=== Summary ==="))
  print(paste0("Total changes proteins in input:    " , length(changes)                  ))
  print(paste0("Present in network_full:            " , length(present_in_igraph)        ))
  print(paste0("Completely missing:                 " , length(missing_in_igraph)        ))
  print(paste0("Present but unconnected (0 edges):  " , length(unconnected_changes_nodes)))
  ################################################################################################
  
  #### save network and layout
  save(network_full, file = paste0(result_directory, "intermediate/full_network.rds"))
  save(layout_full,  file = paste0(result_directory, "intermediate/layout_full.rds" ))
  save(changes,      file = paste0(result_directory, "/intermediate/changes.rds"    ))
  save(igraph_diagnostics,  file = paste0(result_directory, "intermediate/igraph_diagnostics.rds"))
  
  
}

