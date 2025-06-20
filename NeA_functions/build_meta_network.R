# essentially dump together output of other filter_..._network functions and filter for distinct ppi relations
build_meta_network <- function(network_output_filename) {
  
  ## create networks
  # create string network -----------
  string_network        <- "physical"   # "physical" of "full"
  string_score_cutoff   <- 700 # cathy used 700 as default     string combined score ranges 0-1000; 400 threshold for medium confidence, 700 for high confidence;    low throughput high accuracy experiment data results score usually +/- 600; high throughput low accuracy data scores <= 250
  # string_linkage_categs <- c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining")   # select info levels to be used these are available:   c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining") 
  network_output_filename_string <- "network_interactions_string"   # 
  
  filter_string_network(string_network, string_score_cutoff, 
                        changes, 
                        network_output_filename = network_output_filename_string)
 
  # create biogrid network -----------
  complexPortal_input_file         <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/complex_portal_9606.tsv" #" /Users/mgesell/Desktop/currentR/Network_plotting_Cathy/input/complex_portal_9606.tsv"
  complexPortal_mapping_table_file <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/uniprotkb_AND_reviewed_true_AND_model_o_2024_11_08.tsv" #"/Users/mgesell/Desktop/currentR/Network_plotting_Cathy/input/uniprotkb_AND_reviewed_true_AND_model_o_2024_11_08.tsv"
  network_output_filename_biogrid   <- "network_interactions_complexPortal"   # 
  
  filter_complexPortal_network(complexPortal_input_file, complexPortal_mapping_table_file,
                               changes,
                               network_output_filename = network_output_filename_biogrid)
  
  # create complex portal network -----------
  biogrid_input_file              <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/BIOGRID-ALL-4.4.242_20250218.mitab.txt" #"/Users/mgesell/PhD/local_resources/PPIs_and_complexes/BIOGRID-ALL-4.4.242_20250218.mitab.txt"
  network_output_filename_cp  <- "network_interactions_biogrid_physical"   # 
  
  filter_biogrid_network(biogrid_input_file, 
                         changes, 
                         network_output_filename = network_output_filename_cp)
  
  ## load network results networks
  network_string         <- read.table(paste(result_directory, "intermediate/", network_output_filename_string , ".tsv", sep = "" ), header = TRUE)
  network_biogrid        <- read.table(paste(result_directory, "intermediate/", network_output_filename_biogrid, ".tsv", sep = "" ), header = TRUE)
  network_complex_portal <- read.table(paste(result_directory, "intermediate/", network_output_filename_cp     , ".tsv", sep = "" ), header = TRUE)
  
  network_meta <- rbind(network_string, 
                        network_biogrid,
                        network_complex_portal
                        ) %>%
    filter(!is.na(protein1_entry), !is.na(protein2_entry)) %>%  # filter out NAs in mapping ------
    filter(protein1_entry != protein2_entry) %>%                # filter out self-loops  A = A --
    mutate(protein_min = pmin(protein1_entry, protein2_entry),
           protein_max = pmax(protein1_entry, protein2_entry) ) %>%
    distinct(protein_min, protein_max, .keep_all = TRUE) %>%   # kick AB BA duplicates ----------
    dplyr::select(-protein_min, -protein_max) %>%
    distinct() # Remove duplicates

  ## export -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  write.table(network_meta
              , paste(result_directory, "intermediate/", network_output_filename, ".tsv", sep = "" ), quote = F , sep = "\t", row.names = F)
  
  ## step 4: reporting  --------------------------------------------------------------------------------------------------------------------------------------------------------------
  message(paste0("_____ using meta network of ", length(  unique( c(network_meta$protein1_entry, network_meta$protein2_entry) )  )  , " unique proteins and ", nrow(network_meta), " ppis ________ "))
  
  }