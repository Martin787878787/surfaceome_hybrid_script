load_protein_data_diff <- function(directory_input) { ##########################
  FC_cutoff     = 1    #########################################################
  pvalue_cutoff = 0.05 #########################################################
  ##############################################################################
  
  data_prot_diff <- read.csv(paste0(directory_input, "/_data_prot_diff_abundance.csv"), header = TRUE) %>%
    dplyr::select(entry, entry_name, comparison, imputed_comparison, log2FC, adj_pvalue) %>%
    filter(log2FC >= FC_cutoff | log2FC <= FC_cutoff) %>%  # fold change filtering 1 = 2x   0.5849625 = 1.5x
    filter(adj_pvalue <= pvalue_cutoff) %>%         # p value filtering
    distinct()
  
  # sort dataframe by comparison, then determine overlap 
  data_prot_diff <- data_prot_diff %>%
    mutate(comparison = factor(comparison, levels = unique(data_prot_diff$comparison))) %>%
    arrange(comparison) %>%
    group_by(entry) %>%
    mutate(overlap = paste(sort(unique(comparison[!is.na(log2FC) & log2FC != 0])), collapse = ".")) %>%
    ungroup() 
  
  # QC plots
  print(
    ggplot(data_prot_diff, aes(x = comparison, y = log2FC)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "Log2 Fold Change by Comparison",
           x = "Comparison",
           y = "Log2 Fold Change") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  print(
    ggplot(data_prot_diff, aes(x = comparison, y = adj_pvalue)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "Log2 Fold Change by Comparison",
           x = "Comparison",
           y = "Log2 Fold Change") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  print(
    ggplot(data_prot_diff %>%
             group_by(comparison, imputed_comparison) %>%
             summarise(count = n(), .groups = 'drop'), 
           aes(x = comparison, y = count, fill = imputed_comparison)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("yes" = "red", "no" = "black")) +
      theme_minimal() +
      labs(title = "Imputed vs. Not Imputed Counts by Comparison",
           x = "Comparison",
           y = "Count",
           fill = "Imputation Status") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  
  if (length(unique(data_prot_diff$comparison)) > 1) { 
    # Upset plot prep
    comparisons <- unique(data_prot_diff$comparison)
    sets <- lapply(comparisons, function(comp) {
      data_prot_diff$entry[data_prot_diff$comparison == comp & abs(data_prot_diff$log2FC)  > FC_cutoff & data_prot_diff$adj_pvalue < pvalue_cutoff] # abs(data_prot_diff$log2FC) to include all regulated proteins
    })
    names(sets) <- comparisons
    # UpSet plot
    print(
      upset(fromList(sets), 
            nsets = length(comparisons),
            order.by = "freq",
            mainbar.y.label = "Intersection Size",
            sets.x.label = "Set Size",
            text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
            point.size = 3,
            line.size = 1)  
    )
    grid.text("SigRegulated (Up+Down) Overlap", x = 0.65, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    #
    comparisons <- unique(data_prot_diff$comparison)
    sets <- lapply(comparisons, function(comp) {
      data_prot_diff$entry[data_prot_diff$comparison == comp & data_prot_diff$log2FC > FC_cutoff & data_prot_diff$adj_pvalue < pvalue_cutoff] # abs(data_prot_diff$log2FC) to include all regulated proteins
    })
    names(sets) <- comparisons
    # UpSet plot
    print(
      upset(fromList(sets), 
            nsets = length(comparisons),
            order.by = "freq",
            mainbar.y.label = "Intersection Size",
            sets.x.label = "Set Size",
            text.scale = c(3, 3, 3, 3, 3, 3),
            point.size = 3,
            line.size = 1)
    )
    grid.text("SigUp Overlap", x = 0.65, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
  }
  
  return(data_prot_diff)
}


