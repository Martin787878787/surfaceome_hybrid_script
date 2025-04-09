initialize_folders <- function(delete_existing_files = FALSE){
  
  dir.create(paste0(result_directory, "final")       , showWarnings = TRUE, mode = '777')
  dir.create(paste0(result_directory, "intermediate"), showWarnings = TRUE)
  dir.create(paste0(result_directory, "plots")       , showWarnings = TRUE)
  
  if(delete_existing_files){
    unlink(paste0(result_directory, "final/*"))
    unlink(paste0(result_directory, "intermediate/*"))
    unlink(paste0(result_directory, "plots/*"))
  }
  
  
  
}