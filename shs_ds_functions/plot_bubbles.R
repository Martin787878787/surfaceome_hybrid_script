# belongs to shs_downstream
# 
plot_bubbles <- function(data, x_var, y_var, size_var, fill_var, title_var) {
  
  if (is.null(fill_var)) {
    plot <-  ggplot(data, aes(x = !!sym(x_var), y = !!sym(y_var))) +
      geom_point(aes(size = !!sym(size_var)), shape = 21, color = "darkgrey", fill = "black") +
      scale_size_continuous(range = c(1, 10)) + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "GO Term", y = "Condition", size = "Recall", 
           title = title_var) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
  } else {
    plot <-  ggplot(data,
                    aes(x = !!sym(x_var) , y = !!sym(y_var))) +
      geom_point(aes(size = !!sym(size_var), fill = !!sym(fill_var)), shape = 21, color = "darkgrey") +
      scale_size_continuous(range = c(1, 10)) + 
      scale_fill_gradient(low = "firebrick", high = "#fee0d2") + ##fee0d2
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "GO Term", y = "Condition", size = "Recall", fill = "P-value", 
           title = title_var) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  }
  
  return(plot) 
}