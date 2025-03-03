# belongs to shs_downstream
# purpuse = calculate distances for x and y axis data for nice plots
calculate_dual_distance_pattern_magnitude <- function(data, terms, term_column, distance_column, distance_method) {

  # Filter and arrange data
    data <- data %>%
      filter(!!sym(term_column) %in% terms)
    
    # Calculate distance matrix (for similarity clustering)
    data_w <- data %>%
      dplyr::select(comparison, !!sym(term_column), !!sym(distance_column)) %>%
      distinct() %>%
      pivot_wider(names_from = comparison, values_from = !!sym(distance_column), values_fill = 0)
    
    # Create a binary presence/absence matrix
    data_w_binary <- data_w
    data_w_binary[, -1] <- ifelse(data_w_binary[, -1] > 0, 1, 0)
    
    # Calculate distance based on binary data
    dist_matrix   <- dist(data_w_binary[, -1], method = distance_method )# distance_method = "euclidean"
    hclust_result <- hclust(dist_matrix)
    term_order <- data_w %>% pull(!!sym(term_column)) %>% .[hclust_result$order]
    
    # Calculate combined recall for secondary sorting
    data_w$combined_recall <- rowSums(data_w[, -1])
    
    # Order terms based on clustering and combined recall
    term_order_final <- data_w %>%
      mutate(term = !!sym(term_column)) %>%
      arrange(match(term, term_order), desc(combined_recall)) %>%
      pull(term)
    
    # Distance of comparisons
    comparison_dist <- dist(t(data_w_binary[, -1]), method = distance_method)
    comparison_hclust <- hclust(comparison_dist)
    comparison_order <- colnames(data_w)[-1][comparison_hclust$order]
    
    # Order df based on clustering and combined recall
    data_ordered <- data %>%
      mutate(
        term_name = factor(!!sym(term_column), levels = term_order_final),
        comparison = factor(comparison, levels = comparison_order)
      )
    
    return(data_ordered) 
  }
  
  
  
#   
#   
#     
#   # Filter and arrange data
#   data <- data %>%
#     filter(!!sym(term_column) %in% terms)
#   
#   # Calculate distance matrix (for similarity clustering)
#   data_w <- data %>%
#     dplyr::select(comparison, !!sym(term_column), !!sym(distance_column)) %>%
#     distinct() %>%
#     pivot_wider(names_from = comparison, values_from = !!sym(distance_column), values_fill = 0)
#   
#   # Custom distance function
#   custom_dist <- function(x, y) {
#     x <- as.numeric(x)
#     y <- as.numeric(y)
#     x <- x[!is.na(x)]
#     y <- y[!is.na(y)]
#     if (length(x) == 0 || length(y) == 0) return(1)  # Maximum distance if all NA
#     pattern_dist <- 1 - cor(x, y, method = "spearman", use = "pairwise.complete.obs")
#     magnitude_dist <- abs(sum(x, na.rm = TRUE) - sum(y, na.rm = TRUE)) / (sum(x, na.rm = TRUE) + sum(y, na.rm = TRUE))
#     return(pattern_dist + magnitude_dist)
#   }
#   
#   # Calculate distance matrix using custom distance function
#   dist_matrix <- as.dist(outer(1:nrow(data_w), 1:nrow(data_w), 
#                                Vectorize(function(i, j) custom_dist(as.numeric(data_w[i, -1]), as.numeric(data_w[j, -1])))))
#   
#   # Perform hierarchical clustering
#   hclust_result <- hclust(dist_matrix, method = "complete")
#   term_order    <- data_w %>% pull(!!sym(term_column)) %>% .[hclust_result$order]
#   
#   # Distance of comparisons
#   comparison_dist <- as.dist(outer(2:ncol(data_w), 2:ncol(data_w), 
#                                    Vectorize(function(i, j) custom_dist(data_w[[i]], data_w[[j]]))))
#   comparison_hclust <- hclust(comparison_dist, method = "complete")
#   comparison_order <- colnames(data_w)[-1][hclust(comparison_dist, method = "complete")$order]
#   
#   # Order df based on clustering
#   data_ordered <- data %>%
#     mutate(
#       term_name  = factor(!!sym(term_column), levels = term_order),
#       comparison = factor(comparison, levels = comparison_order)
#     )
#   
#   return(data_ordered) 
# }





  
#   
#   
#   
#   function(data, terms, term_column, distance_column, distance_method) {
#   
#   # filter, arrange, factor for plotting
#   data <- data %>%
#     filter(!!sym(term_column) %in% terms)
#   
#   # calculate distance matrix (for similarity clustering)
#   data_w <- data %>%
#     dplyr::select(comparison, !!sym(term_column), !!sym(distance_column)) %>%
#     # replace any numeric value with 1 such that presence absence based clustering??
#     distinct() %>%
#     pivot_wider(names_from = comparison, values_from =  !!sym(distance_column), values_fill = 0)  # !!sym(distance_column)
#   dist_matrix        <- dist(data_w[, -1], method = distance_method)    # manhattan # euclidean
#   hclust_result      <- hclust(dist_matrix)       
#   term_order         <- data_w %>% pull(!!sym(term_column)) %>% .[hclust_result$order]
#   
#   # distance of cmparisons
#   comparison_dist   <- dist(t(data_w[, -1]), method = distance_method)  # manhattan # euclidean
#   comparison_hclust <- hclust(comparison_dist)
#   comparison_order  <- colnames(data_w)[-1][comparison_hclust$order]
#   
#   # order df based on clustering
#   data_ordered <- data %>%
#     mutate(term_name  = factor(data %>% pull(!!sym(term_column)), levels = term_order      ),
#            comparison = factor(data %>% pull(comparison)        , levels = comparison_order)  )
#   
#   return(data_ordered) 
# }