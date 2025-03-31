# belongs to shs_downstream
load_protein_data <- function(directory_input, condition_levels) {
  data_prot <- read.csv(paste0(directory_input, "/_data_prot_level.csv"), header = TRUE) %>%
    dplyr::select(entry, condition, con_rep, imputed_prot_intensity_log2) %>%
    distinct() %>%
    group_by(entry, condition) %>%
    mutate(log2_median = median(imputed_prot_intensity_log2, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(entry) %>%
    mutate(
      log2_mean_protein = mean(imputed_prot_intensity_log2, na.rm = TRUE),
      log2_sd_protein   = sd(imputed_prot_intensity_log2, na.rm = TRUE),
      log2_z_score      = ifelse(
        !is.na(log2_median) & log2_sd_protein != 0,
        (log2_median - log2_mean_protein) / log2_sd_protein,
        NA
      )
    ) %>%
    ungroup() %>%
    dplyr::select(-con_rep, -imputed_prot_intensity_log2, -log2_mean_protein, -log2_sd_protein) %>%
    distinct()
  
  # Complete data frame for all entry-condition combinations
  all_combinations <- expand.grid(
    entry = unique(data_prot$entry),
    condition = unique(data_prot$condition)
  )
  data_prot <- all_combinations %>%
    left_join(data_prot, by = c("entry", "condition")) %>%
    mutate(
      log2_median  = coalesce(log2_median , NA),
      log2_z_score = coalesce(log2_z_score, NA)
    ) %>%
    mutate(imputed = if_else(is.na(log2_median), TRUE, FALSE))
  
  # sort dataframe by conditions, then determine overlap 
  data_prot <- data_prot %>%
    mutate(condition = factor(condition, levels = condition_levels)) %>%
    arrange(condition) %>%
    group_by(entry) %>%
    mutate(overlap = paste(sort(unique(condition[!is.na(log2_median) & log2_median != 0])), collapse = "_")) %>%
    ungroup() 
  
  # impute NA (after completing df)
  NA_impute_zscore <- min(data_prot$log2_z_score, na.rm = TRUE) - # substract 10% data offset from min for NA imputation (white in plot)
    (max(data_prot$log2_z_score, na.rm = TRUE) - min(data_prot$log2_z_score, na.rm = TRUE))*0.1
  data_prot <- data_prot %>%
    mutate(log2_z_score    = replace_na(log2_z_score   , NA_impute_zscore),   # zscore impute by value 10% lower from min than overall data range
           log2_median = replace_na(log2_median, 0                   ))   # median log2 impute by 0# LUX ---
  
  # z-score qc output
  boxplot(data_prot$log2_z_score)
  boxplot(data_prot$log2_median)
  print(summary(data_prot$log2_z_score))
  
  # ===========================================================================================================================================================================
  # QC plot IDs/condition 
  print(
    ggplot(data_prot %>% 
             filter(!imputed) %>%
             group_by(condition) %>%
             summarize(unique_entries = n_distinct(entry)) %>%
             mutate(condition = factor(condition, levels = condition_levels)), # Convert condition to a factor with the specified order
           aes(x = condition, y = unique_entries)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(x = "Condition", y = "Number of Unique Entries", 
           title = "Unique Entries per Condition") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  # median intensity plots
  print(
    ggplot(data_prot %>% filter(log2_median > 0), aes(x = condition, y = log2_median)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "log2_median Fold Change by Comparison",
           x = "Comparison",
           y = "Log2 Fold Change") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  # Upset
  conditions <- unique(data_prot$condition)
  sets <- lapply(conditions, function(condi) {
    data_prot %>% filter(log2_median != 0 & data_prot$condition == condi) %>% pull(entry)
  })
  names(sets) <- conditions
  # UpSet plot
  print(
    upset(fromList(sets), 
          nsets = length(conditions),
          order.by = "freq",
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size",
          text.scale = c(3, 3, 3, 3, 3, 3),
          point.size = 3,
          line.size = 1)
  )
  grid.text("Protein Level ID Overlap", x = 0.65, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
  
  
  return(data_prot)
  
}
