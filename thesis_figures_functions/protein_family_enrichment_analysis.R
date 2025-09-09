# belongs to thesis figures functions
protein_family_enrichment_analysis <- function(plot_heading_name, protein_fam_analysis_entry_input_df) {

  # load resource
  protein_family_upsp_resource <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/protein_families/human_upsp_proteinfamilies_2025_07_23.csv", header = TRUE) %>%
    dplyr::rename("entry" = "Entry","protein_families"= "Protein.families") %>%
    dplyr::select(entry, protein_families) %>%
    mutate(protein_families = as.character(protein_families)) %>%
    filter(!is.na(protein_families) & protein_families != "")
  # map protein families to data
  protein_fam_data_annotated <- protein_fam_analysis_entry_input_df %>%
    left_join(protein_family_upsp_resource, by = "entry") %>%
    as.data.frame() %>%
    filter(!is.na(protein_families) & protein_families != "")
  
  # extract protein family info from rows - each name entry pair by expansion to multiple rows were needed to not miss any terms ---------
  protein_family_upsp_resource_unpacked <- protein_family_upsp_resource %>%
    separate_rows(protein_families, sep = ";") %>%
    filter(!is.na(protein_families) & protein_families != "") %>%
    as.data.frame() %>%
    distinct()
  # write.csv(protein_family_upsp_resource_unpacked %>% dplyr::select(protein_families) %>% distinct(), "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/thesis_figures_functions/unique_protein_families.csv", row.names=FALSE)
  #
  protein_fam_data_annotated_unpacked <- protein_fam_data_annotated %>%
    separate_rows(protein_families, sep = ";") %>%
    filter(!is.na(protein_families) & protein_families != "") %>%
    distinct()
  #
  cat(nrow(protein_fam_data_annotated_unpacked) - nrow(protein_fam_data_annotated), "mutli-term annotations expanded indataset\n")
  
  # analysis ---------
  protein_family_analysis <- as.data.frame(table(protein_fam_data_annotated_unpacked$protein_families)) %>%
    setNames(c("protein_families", "data_group_size")) %>% 
    left_join(
      as.data.frame(table(protein_family_upsp_resource_unpacked$protein_families)) %>%
        setNames(c("protein_families", "bg_upsp_group_size")),
      by = "protein_families"
    ) %>%
    filter(protein_families != "", !is.na(data_group_size), !is.na(bg_upsp_group_size)) %>%
    mutate(
      total_detected   = sum(data_group_size, na.rm = TRUE),
      total_background = sum(bg_upsp_group_size, na.rm = TRUE),
      recall           = data_group_size / bg_upsp_group_size
    ) %>%
    mutate(
      detected_freq   = data_group_size    / total_detected,
      background_freq = bg_upsp_group_size / total_background,
      enrichment      = if_else(bg_upsp_group_size == 0, NA_real_, detected_freq / background_freq)
    ) %>%
    # Calculate p-value using Fisher's exact test
    rowwise() %>%
    mutate(
      p_value = fisher.test(
        matrix( c(data_group_size, 
                  bg_upsp_group_size, 
                  total_detected - data_group_size, 
                  total_background - bg_upsp_group_size),
                nrow = 2,
                byrow = FALSE),
        alternative = "greater")$p.value, # single sided test = "greater"
      Significance = if_else(p_value < 0.05, "yes", "no") ) %>%
    ungroup() %>%
    dplyr::select(-total_detected, -total_background) 
    
  
  protein_family_analysis <- protein_family_analysis %>%
    mutate(
      log2FC    = log2(enrichment),
      negLog10P = -log10(p_value),
      plot_label = case_when(
        log2FC >= filter_log2fc_cutoff & negLog10P >= -log10(filter_sig_cutoff) ~ as.character(protein_families),
        TRUE ~ ""),
      plot_label = as.character(plot_label)
    )

  ## plotting ..................................................................................................
  foldchangelimit_max <- max(protein_family_analysis$log2FC, na.rm = TRUE)
  foldchangelimit_min <- min(protein_family_analysis$log2FC, na.rm = TRUE)
  max_pvalue          <- max(-log10(protein_family_analysis$p_value), na.rm = TRUE)
  # Helper frame for significance lines
  segmentation <- data.frame(x = c(foldchangelimit_min, filter_log2fc_cutoff, -filter_log2fc_cutoff, filter_log2fc_cutoff),
                             y = c(-log10(filter_sig_cutoff), -log10(filter_sig_cutoff), -log10(filter_sig_cutoff), -log10(filter_sig_cutoff)),
                             xend = c(-filter_log2fc_cutoff, foldchangelimit_max, -filter_log2fc_cutoff, filter_log2fc_cutoff),
                             yend = c(-log10(filter_sig_cutoff), -log10(filter_sig_cutoff), max_pvalue, max_pvalue),
                             col = rep("black", times=4),
                             linetype = rep("dashed", times=4))
  
  plot <- ggplot(protein_family_analysis, aes(x = log2FC, y = negLog10P, label = plot_label)) +
    geom_point(aes(fill = Significance), color = "white", shape = 21, size = 3, stroke = 0.6) +
    scale_fill_manual(values = c("yes" = "black", "no" = "grey")) +
    geom_segment(x=segmentation$x[1], y=segmentation$y[1], xend=segmentation$xend[1], yend=segmentation$yend[1], linetype="dashed", col="grey", linewidth = 0.4) +
    geom_segment(x=segmentation$x[2], y=segmentation$y[2], xend=segmentation$xend[2], yend=segmentation$yend[2], linetype="dashed", col="grey", linewidth = 0.4) +
    geom_segment(x=segmentation$x[3], y=segmentation$y[3], xend=segmentation$xend[3], yend=segmentation$yend[3], linetype="dashed", col="grey", linewidth = 0.4) +
    geom_segment(x=segmentation$x[4], y=segmentation$y[4], xend=segmentation$xend[4], yend=segmentation$yend[4], linetype="dashed", col="grey", linewidth = 0.4) +
    geom_text(aes(label = plot_label), hjust = 1) +
    xlab("log2(Fold Change)") +
    ylab("-log10(p-value)") +
    ggtitle(paste0(plot_heading_name)) +    # <-- Add this line for the title
    plot_theme()   
  
  sig_enrichment_table <- protein_family_analysis %>% filter(log2FC >= 1, p_value <= 0.05) %>% arrange(desc(log2FC)) %>%
    dplyr::select(-c(detected_freq, background_freq, enrichment, Significance, plot_label)) %>%
    mutate(plot_heading = plot_heading_name)
  
  sig_enrichment_proteins <- protein_fam_data_annotated_unpacked %>% 
    dplyr::filter(protein_families %in% sig_enrichment_table$protein_families) %>%
    left_join(sig_enrichment_table %>% dplyr::select(protein_families, log2FC), by = "protein_families") %>%
    arrange(desc(log2FC))  %>%
    dplyr::select(-log2FC) %>%
    mutate(plot_heading = plot_heading_name)
    
  result <- list(plot, sig_enrichment_table, sig_enrichment_proteins)
  names(result) <- c(paste0(plot_heading_name, "__plot"     ),
                     paste0(plot_heading_name, "__table"    ),
                     paste0(plot_heading_name, "__proteins"))
  
  
  return(result)
  
}
