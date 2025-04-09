calculate_overlaps <- function(data, id_col, condition_col) {
  # Get unique conditions
  conditions <- unique(data[[condition_col]])
  n_conditions <- length(conditions)
  
  if (n_conditions <= 1) {
    stop("Number of unique conditions must be greater than 1")
  }
  
  # Generate all possible combinations of conditions, sorted by number of conditions
  combinations <- as.data.frame(expand.grid(replicate(n_conditions, c(1, 0), simplify = FALSE)))
  combinations <- combinations[order(rowSums(combinations), decreasing = TRUE), ]

  
  # Initialize a list to store assigned entries
  assigned_entries <- list()
  
  # Calculate overlaps for each combination
  overlaps <- lapply(1:nrow(combinations), function(i) {
    selected_conditions <- conditions[as.logical(combinations[i,])]
    if (length(selected_conditions) > 0) {
      entries <- data[data[[condition_col]] %in% selected_conditions, id_col]
      unique_entries <- setdiff(unique(entries), unlist(assigned_entries))
      if (length(unique_entries) > 0) {
        overlap <- Reduce(intersect, split(data[[id_col]], data[[condition_col]])[selected_conditions])
        overlap <- setdiff(overlap, unlist(assigned_entries))
        assigned_entries <<- c(assigned_entries, list(overlap))
        list(
          conditions = paste(selected_conditions, collapse = " & "),
          entries = paste(unique_entries, collapse = ", "),
          overlap = paste(overlap, collapse = ", "),
          count = length(unique_entries),
          overlap_count = length(overlap)
        )
      }
    }
  })
  
  # Remove NULL entries
  overlaps <- overlaps[!sapply(overlaps, is.null)]
  
  # Create result dataframe
  result_df <- do.call(rbind, lapply(seq_along(overlaps), function(i) {
    x <- overlaps[[i]]
    data.frame(
      row = i,
      conditions = x$conditions,
      entries = x$entries,
      overlap = x$overlap,
      count = x$count,
      overlap_count = x$overlap_count,
      stringsAsFactors = FALSE
    )
  }))
  
  # Reset row names
  rownames(result_df) <- NULL
  
  return(result_df)
}
