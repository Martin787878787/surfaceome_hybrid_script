pF_pD_enrichment_analysis <- function(data_subset_name, protein_entry_vector, mode = c("protein_family", "domain_ft"),
    filter_log2fc_cutoff = 1,  filter_sig_cutoff = 0.05) {
  
  # -------- test input ----------------------------------------------------------------------------
  # mode = "domain_ft"
  # data_subset_name <- "panT"
  # protein_entry_vector <- v31_LUX_data_prot_diff_abundance_sigup_ssEXTENDED %>%
  #   dplyr::filter(plot_heading == "panT") %>% # comment this line out to get full dataset view
  #   pull(entry) %>%
  #   unique()
  # ------------------------------------------------------------------------------------------------
  
  mode <- match.arg(mode)
  resource_file = "/Users/mgesell/Desktop/currentR/git/shs_resources/human_upsp_202508.csv"
    
  # Mode-specific configuration
  term_col <- switch(mode,
                     protein_family = "protein_families",
                     domain_ft      = "domain_ft")
  
  # Load resource with relevant annotation column-------------------------------
  resource <- read_protti(resource_file, header = TRUE) %>%
    dplyr::select(entry, entry_name, !!sym(term_col)) %>%
    mutate(term = as.character(!!sym(term_col))) %>%
    filter(!is.na(term) & term != "")
  
  # Unpack multi-term annotations -----------------------------------------------
  resource_unpacked <- resource %>%
    separate_rows(term, sep = ";") %>%
    filter(!is.na(term) & term != "") %>%
    as.data.frame() %>% 
    mutate(term = gsub("^ "     , "", term)) %>%   # kick space 
    distinct()
  #
  if (mode == "domain_ft") { # domain info specific trimming of info
    resource_unpacked <- resource_unpacked %>%
      filter(str_detect(term, "/note=")) %>%        # kick other elements of annotation (DOMAIN, note, EVIDENCE,)
      mutate(term = gsub("/note=", "", term),       # trim
             term = gsub("\""    , "", term)) %>%   # trim
      distinct() %>% # require distinct before removing domain occurance count (next row) --> # of domains / protein is NOT reduced to 1 (unique occurrence)
      mutate(term = sub(" \\d+$" , "", term)) # trim to remove domain count e.g. SUSHI 12 and SUSHI 23 are all SUSHI
  }
  
  # subset the resource to the proteins detected in current query data
  data_annotated <- resource_unpacked %>%
      dplyr::filter(entry %in% protein_entry_vector)
  
  # Enrichment analysis ============================================================
  enrichment_analysis <- as.data.frame(table(data_annotated$term)) %>%
    setNames(c("term", "data_group_size")) %>% 
    left_join(
      as.data.frame(table(resource_unpacked$term)) %>%
        setNames(c("term", "bg_resource_group_size")),
      by = "term"
    ) %>%
    filter(term != "", !is.na(data_group_size), !is.na(bg_resource_group_size)) %>%
    mutate(
      total_detected   = sum(data_group_size, na.rm = TRUE),
      total_background = sum(bg_resource_group_size, na.rm = TRUE),
      recall           = data_group_size / bg_resource_group_size
    ) %>%
    mutate(
      detected_freq   = data_group_size    / total_detected,
      background_freq = bg_resource_group_size / total_background,
      enrichment      = if_else(bg_resource_group_size == 0, NA_real_, detected_freq / background_freq)
    ) %>% 
    rowwise() %>%
    mutate(
      p_value = fisher.test(
        matrix(c(data_group_size, bg_resource_group_size, total_detected - data_group_size, total_background - bg_resource_group_size),
               nrow = 2, byrow = FALSE), alternative = "greater")$p.value,
      Significance = if_else(p_value < filter_sig_cutoff, "yes", "no")
    ) %>%
    ungroup() %>%
    dplyr::select(-total_detected, -total_background)  
  
  enrichment_analysis <- enrichment_analysis %>%
    mutate(
      log2FC    = log2(enrichment),
      negLog10P = -log10(p_value),
      plot_label = case_when(
        log2FC >= filter_log2fc_cutoff & negLog10P >= -log10(filter_sig_cutoff) ~ as.character(term),
        TRUE ~ ""),
      plot_label = as.character(plot_label)
    )
  
  ## plotting
  foldchangelimit_max <- max(enrichment_analysis$log2FC, na.rm = TRUE)
  foldchangelimit_min <- min(enrichment_analysis$log2FC, na.rm = TRUE)
  max_pvalue          <- max(-log10(enrichment_analysis$p_value), na.rm = TRUE)
  segmentation <- data.frame(
    x = c(foldchangelimit_min, filter_log2fc_cutoff, -filter_log2fc_cutoff, filter_log2fc_cutoff),
    y = rep(-log10(filter_sig_cutoff), 4),
    xend = c(-filter_log2fc_cutoff, foldchangelimit_max, -filter_log2fc_cutoff, filter_log2fc_cutoff),
    yend = c(-log10(filter_sig_cutoff), -log10(filter_sig_cutoff), max_pvalue, max_pvalue),
    col = rep("black", 4),
    linetype = rep("dashed", 4)
  )
  
  plot <- ggplot(enrichment_analysis, aes(x = log2FC, y = negLog10P)) +
    geom_point(aes(fill = Significance), color = "white", shape = 21, size = 3, stroke = 0.6) +
    scale_fill_manual(values = c("yes" = "black", "no" = "grey")) +
    geom_segment(data = segmentation, 
                 aes(x = x, y = y, xend = xend, yend = yend), linetype = "dashed", col = "grey", linewidth = 0.4) +
    geom_text(aes(label = plot_label), hjust = 1) +
    xlab("log2(Fold Change)") +
    ylab("-log10(p-value)") +
    ggtitle(data_subset_name) +
    plot_theme()
  
  sig_enrichment_table <- enrichment_analysis %>% filter(log2FC >= 1, p_value <= 0.05) %>% arrange(desc(log2FC)) %>%
    dplyr::select(-c(detected_freq, background_freq, enrichment, Significance, plot_label)) %>%
    mutate(plot_heading = data_subset_name,
           term_summary = paste0(round(recall*100, 0), "% (", data_group_size, " of ", bg_resource_group_size,")")  ) %>%
    dplyr::select(plot_heading, term, term_summary, data_group_size, bg_resource_group_size, recall, p_value, log2FC, negLog10P)
  
  sig_enrichment_proteins <- data_annotated %>% 
    distinct() %>% # this removes double occurances for domains thate were before domain_x_occurance1, domain_x_occcurance2
    filter(term %in% unique(sig_enrichment_table$term)) %>%
    left_join(sig_enrichment_table %>% dplyr::select(term, log2FC), by = "term") %>%
    arrange(desc(log2FC)) %>%
    dplyr::select(-log2FC) %>%
    mutate(plot_heading = data_subset_name)
  
  result <- list(plot, sig_enrichment_table, sig_enrichment_proteins)
  names(result) <- c(paste0(data_subset_name, "__plot"),
                     paste0(data_subset_name, "__table"),
                     paste0(data_subset_name, "__proteins"))
  return(result)
}
