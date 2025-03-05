
plot_dual_distance_bubble <- function(data, group_filter, min_recall, min_p_value, 
                                      grouping, term_column, distance_column, distance_method,  # distance parameters
                                      x_var, y_var, size_var, fill_var, title_var) {                             # plot parameters
 
  ## prepare data oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # subset data
  if (is.null(group_filter)) {
    grouping_levels <- "clustering"  # no data subset used --> plotting order will rely on clustering of grouping
  } else {  
    data <- data %>% 
      filter(!!sym(grouping) %in% group_filter)
    grouping_levels <- group_filter               # recycle group_filter as data levels if given
  }

  # only terms that meat filter criteria in at least one condition shall be displayed in plot
  if (is.null(min_p_value)) {
    terms_oi <- data %>%
      filter(!!sym(distance_column) >= min_recall) %>%
      pull(!!sym(term_column)) %>%
      unique()
  } else {
  terms_oi <- data %>%
    filter(!!sym(distance_column)  >= min_recall & p_value <= min_p_value) %>%
    pull(!!sym(term_column)) %>%
    unique()
  }
  # Filter and arrange data
  data <- data %>%
    filter(!!sym(term_column) %in% terms_oi)
  # in case one condition does not have any terms matching the filter conditions its lost 
  # --> below recreate artifical data for that group so it is displayed in plot (recall = 0 sufficient to not see in plot but still be displayed on axis)
  if(!is.null(setdiff(group_filter, unique(data %>% pull(!!sym(grouping))))) && 
              length(setdiff(group_filter, unique(data[[grouping]]))) > 0   ) { 
    data_dummy <- data %>%
      mutate(!!sym(grouping)        := setdiff(group_filter, unique(data %>% pull(!!sym(grouping)))),
             !!sym(distance_column) := 0) %>%
      distinct()
    data <- rbind(data, data_dummy) # append dummy rows 
  }
  ## distance calculation oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # brind df to wide format for distance calucation
  data_w <- data %>%
    dplyr::select(!!sym(grouping), !!sym(term_column), !!sym(distance_column)) %>%
    distinct() %>%
    pivot_wider(names_from = !!sym(grouping), values_from = !!sym(distance_column), values_fill = 0)
  
  # Create a binary presence/absence matrix
  data_w_binary <- data_w
  data_w_binary[, -1] <- ifelse(data_w_binary[, -1] > 0, 1, 0)                                # distance layer 1 = presence absence pattern
  
  # Calculate distance based on binary data
  dist_matrix <- dist(data_w_binary[, -1], method = distance_method)
  hclust_result <- hclust(dist_matrix)
  term_order <- data_w %>% pull(!!sym(term_column)) %>% .[hclust_result$order]
  
  # Calculate combined recall for secondary sorting
  data_w$combined_recall <- rowSums(data_w[, -1]) 
  
  # Calculate presence/absence distances
  presence_absence_dist <- as.matrix(dist(data_w_binary[, -1], method = distance_method))     # distance layer 2 = optional = combined recall
  
  # Custom ordering function - here the mandatory (presence absence = layer 1 ) and optional (combined recall) layer get combined
  # optional distance affects order only where the binary distance is same for term_column values
  custom_order <- function(term_order, presence_absence_dist, combined_recall) {
    n               <- length(term_order)  # number of terms
    ordered_indices <- order(term_order)   # order
    
    for (i in 2:n) { # loop trought n terms  & start from second term
      j <- i         # loop counter
      # inner loop
      # loop moves backwards through the ordered terms, until ...   We haven't reached the beginning (j > 1)    or   presence/absence distance between the current term and the previous term is zero.
      while (j > 1 && presence_absence_dist[ordered_indices[j-1], ordered_indices[j]] == 0) {  
         # If the combined recall of the current term is higher than the previous term, their positions are swapped.
         if (combined_recall[ordered_indices[j-1]] < combined_recall[ordered_indices[j]]) {
          temp <- ordered_indices[j-1]
          ordered_indices[j-1] <- ordered_indices[j]
          ordered_indices[j] <- temp
         }
        # move to the next pair of terms (moving backwards)
        j <- j - 1
      }
    }
    
    return(ordered_indices)
  }
  
  # Apply custom ordering
  term_order_indices <- custom_order(
    match(data_w[[term_column]], term_order),
    presence_absence_dist,
    data_w$combined_recall
  )
  
  term_order_final <- data_w[[term_column]][term_order_indices]
  
  # Distance of comparisons
  if (ncol(data_w) > 3) {  # if multiple comparisons present (term_name combined_recall comparison_1 comparison_2 --> n = 4)
    comparison_dist <- dist(t(data_w_binary[, -1]), method = distance_method)
    comparison_hclust <- hclust(comparison_dist)
    comparison_order <- colnames(data_w)[-1][comparison_hclust$order]
  } else {
    comparison_order <- colnames(data_w)[-1]
  }
  
  # order groups based on grouping or clustering according to input param
  if (length(grouping_levels) == 1 && grouping_levels == "clustering") {
    data_ordered <- data %>%
      mutate(
        !!sym(term_column) := factor(!!sym(term_column), levels = term_order_final),
        !!sym(grouping)    := factor(!!sym(grouping)   , levels = comparison_order)  # clustering
      )
  } else { 
  data_ordered <- data %>%
    mutate(
      !!sym(term_column) := factor(!!sym(term_column), levels = term_order_final),
      !!sym(grouping)    := factor(!!sym(grouping)   , levels = grouping_levels)
    )
  }
  
  # plotting ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  if (is.null(fill_var)) {  # if there is no p_value or other continuous column to define color gradient 
    plot <- ggplot(data_ordered, aes(x = !!sym(x_var), y = !!sym(y_var))) +
      geom_point(aes(size = !!sym(size_var)), shape = 21, color = "darkgrey", fill = "black") +
      scale_size(range = c(0, 10), limits = c(0.33, 1)) +
      labs(x = "GO Term", y = "Condition", size = "Recall", title = title_var) +  
      theme(panel.grid.major = element_line(color = "darkgrey", linetype = "dotted"), panel.grid.minor = element_line(color = "grey", linetype = "dotted")    )
    # plot <-  ggplot(data_ordered, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    #   geom_point(aes(size = !!sym(size_var)), shape = 21, color = "darkgrey", fill = "black") +
    #   scale_size(range = c(0, 10), limits = c(0.33, 1)) + # ensure size scales comparable accross plots (same recall is same size) independent of min(recall) of data subset
    #   theme_minimal() +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #   labs(x = "GO Term", y = "Condition", size = "Recall", 
    #        title = title_var) +
    #   theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
  } else {                  # p_value or other continuous column defines color gradient 
    plot <-  ggplot(data_ordered, 
                    aes(x = !!sym(x_var), y = !!sym(y_var))) +
      geom_point(aes(size = !!sym(size_var), fill = !!sym(fill_var)), shape = 21, color = "darkgrey") +
      scale_size(range = c(0, 10), limits = c(0.33, 1)) +
      scale_fill_gradient(low = "firebrick", high = "#fee0d2") +
      labs(x = "GO Term", y = "Condition", size = "Recall", fill = "P-value", 
           title = title_var) +  
      theme(panel.grid.major = element_line(color = "darkgrey", linetype = "dotted"), panel.grid.minor = element_line(color = "grey", linetype = "dotted")    )
    
    # plot <-  ggplot(data_ordered, 
    #                 aes(x = !!sym(x_var) , y = !!sym(y_var))) +
    #   geom_point(aes(size = !!sym(size_var), fill = !!sym(fill_var)), shape = 21, color = "darkgrey") +
    #   scale_size(range = c(0, 10), limits = c(0.33, 1)) + # ensure size scales comparable accross plots (same recall is same size) independent of min(recall) of data subset
    #   scale_fill_gradient(low = "firebrick", high = "#fee0d2") + ##fee0d2
    #   theme_minimal() +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #   labs(x = "GO Term", y = "Condition", size = "Recall", fill = "P-value", 
    #        title = title_var) +
    #   theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  }
  
  return(plot) 
  
}