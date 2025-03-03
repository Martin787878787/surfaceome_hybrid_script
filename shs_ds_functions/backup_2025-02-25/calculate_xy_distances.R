# belongs to shs_downstream
# purpuse = calculate distances for x and y axis data for nice plots
calculate_xy_distances <- function(data, terms, term_column, distance_column, distance_method) {
  
  # filter, arrange, factor for plotting
  data <- data %>%
    filter(!!sym(term_column) %in% terms)
  
  # calculate distance matrix (for similarity clustering)
  data_w <- data %>%
    dplyr::select(comparison, !!sym(term_column), !!sym(distance_column)) %>%
    distinct() %>%
    pivot_wider(names_from = comparison, values_from =  !!sym(distance_column), values_fill = 0)  # !!sym(distance_column)
  dist_matrix        <- dist(data_w[, -1], method = distance_method)    # manhattan # euclidean
  hclust_result      <- hclust(dist_matrix)       
  term_order         <- data_w %>% pull(!!sym(term_column)) %>% .[hclust_result$order]
  
  # distance of cmparisons
  comparison_dist   <- dist(t(data_w[, -1]), method = distance_method)  # manhattan # euclidean
  comparison_hclust <- hclust(comparison_dist)
  comparison_order  <- colnames(data_w)[-1][comparison_hclust$order]
  
  # order df based on clustering
  data_ordered <- data %>%
    mutate(term_name  = factor(data %>% pull(!!sym(term_column)), levels = term_order      ),
           comparison = factor(data %>% pull(comparison)        , levels = comparison_order)  )
  
  return(data_ordered) 
}