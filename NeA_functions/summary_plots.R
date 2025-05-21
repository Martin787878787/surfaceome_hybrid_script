summary_plots <- function(){
  ##### load data
  load(file =                 paste0(result_directory, "intermediate/", "Changes.rds"            ))
  propagated_proteins <- fread(paste0(result_directory, "intermediate/", "propagated_proteins.tsv" ))
  load(file =                 paste0(result_directory, "final/"       , "cluster_sig_network.rds"))
  load(file =                 paste0(result_directory, "intermediate/", "full_network.rds"       ))
  
  
  
  #### plot 1: number of proteins at different steps
  number_of_input_hits              <- length(intersect(changes, names(V(network_full))))
  number_of_propagated_proteins      <- length(unique(propagated_proteins$x))
  number_of_proteins_in_sig_custers <- length(unique(V(network_clustered)$name))
  
  plot_1_data <- data.frame("category" = c("PBI Marker Change",
                                           "Propagated Network",
                                           "Final Clusters"),
                            "number_of_proteins" = c(number_of_input_hits,
                                                     number_of_propagated_proteins,
                                                     number_of_proteins_in_sig_custers))
  
  plot_1_data$category <- factor(plot_1_data$category, levels = plot_1_data$category)
  
  p <- ggplot(plot_1_data, aes(x = category, y = number_of_proteins)) +
    geom_bar(stat="identity", fill = c("#0057b7", "#0057b7", "#ffc900")) +
    geom_text(aes(label=number_of_proteins), vjust=1.6, color="white", size=8)+
    theme_classic() +
    xlab( "Number of proteins in ...") +
    ggtitle("Number of proteins at key steps") +
    theme( axis.title.y=element_blank(), legend.title = element_blank(),
           axis.title.x = element_text(size =16), text = element_text(size = 20))
  ggsave(paste0(result_directory, "plots/", "number_of_proteins_at_key_steps.pdf"), p, width = 21, height = 21, units = "cm")
  
  
  ### plot overview_ cluster size
  plot_2_data <- as.data.frame(table(V(network_clustered)$Cluster_membership))
  
    p <- ggplot(plot_2_data, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = plot_2_data$Var1) +
    geom_text(aes(label=Freq), vjust=1.6, color="white", size=8)+
    theme_classic() +
    xlab( "Cluster") +
    ylab("Number of proteins per cluster") +
    theme( axis.title.y=element_blank(), legend.position = "none",
           axis.title.x = element_text(size =16), text = element_text(size = 20), axis.text = element_text(size = 24))
  ggsave(paste0(result_directory, "plots/", "number_of_proteins_per_cluster.pdf"), p, width = 60, height = 21, units = "cm")
  
  
}
