rm(list = ls())
#
library('data.table')
library('ggplot2')
library("RColorBrewer")
library("plyr")
library("dplyr")
library("stringr")
library("purrr")
library("igraph")
library("MCL")
library("protti")
set.seed(123)
#######################################################################################################################################################################################
#######################################################################################################################################################################################
script_version   = "NeA1.4" # Network Analysis human 1.0
ppi_network      = "meta_pS700_CP_pBG"  # "string", "complexPortal", "biogrid_physical", "meta_pS700_CP_pBG" (=physical string 700 cutoff, complex portal, physical biogrid)
cluster_method   = "walktrap"          # "walktrap", "markov", "betweenness"
# write results to ...
output_dir       = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements"
experiment_name  = "v31_NeA"
test_parameter   = "" # e.g. physical / __CD4network / __CD8network
today_str        = format(Sys.Date(), "%Y-%m-%d")
base_directory   =  paste0(output_dir, "/", experiment_name, "/", today_str, "_", ppi_network, "_", test_parameter)
dir.create(paste0(base_directory), showWarnings = FALSE)
## PageRank
pageRank_dampening_values = c(0.85)  # default 0.85
pageRank_cutoff_values    = c(1)     # 0.9, 0.95, 0.975, 0.99, 1)  # the higher to more stringent the protein (shell of grey proteins). e.g. set to 1 and setpageRank_retain_all_changes = TRUE to see plots for only input list (no interactors)
pageRank_retain_all_changes = TRUE   # chose based on your goal: TRUE~"keep all input nodes"; FALSE~"keep only nodes with highest connectivity 
## defines max cluster size (all larger than this excluded for downstream analysis)
param = 50 # start high (e.g. 500) than check "number_of_proteins_per_cluster.pdf" for orientation & titrate down (e.g. 100, 30) until "clustered_network_significant_with_changes_labeled.pdf" looks nice.  

###### "changes" input data ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# option 1 direct input (default)
changes   <- read.delim("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/LUX_metachanges_shs2.22.csv", header = FALSE, stringsAsFactors = FALSE, sep = ',') %>% pull(1)

# option 2 & 3    2: prior knowlege extraction for T cell linages.    3: prior knowledge + condition (T cell subset) specific information retrieval
  # result for 2: .../intermediate/network_interactions_meta_pS700_CP_pBG.tsv
  # result for 3: .../intermediate/network_interactions_meta_pS700_CP_pBG_changes_sub-network.tsv
# keep both below empty unless specific subset network / subset prior knowledge network extraction desired
subset     = "" # make sure to match subset if desired
LUX_subset = ""  # "nCD4_TCR_vs_nCD4_Iso"   "nnCD4_TCR_vs_nnCD4_Iso"   "nCD8_TCR_vs_nCD8_Iso"    "nnCD8_TCR_vs_nnCD8_Iso"

## extract indirect pior knowledge TCR ppi network ..........................................................................................................................................................
if (subset == "CD4") {
  changes   <- read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_NeA/2025-05_CD4_CD8_network_proxis.csv", header = TRUE) %>%
    filter(!is.na(entry_CD4)) %>% select(entry_CD4) %>% pull(1)
} else if (subset == "CD8") {
  changes   <- read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_NeA/2025-05_CD4_CD8_network_proxis.csv", header = TRUE) %>%
    filter(!is.na(entry_CD8)) %>% select(entry_CD8) %>% pull(1)
}

if (LUX_subset != "") {
  # ## extract LUX sig-up abTCR subset network
  data_prot_diff_v31_LUX_sigup <- read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_uniqueImpute_shs2.22/_data_prot_diff_abundance.csv", header = TRUE) %>%
    filter(log2FC >= 1 & adj_pvalue <= 0.05 &   # extract sig up and ... 
             comparison == LUX_subset) %>%      # condition specific data  
    select(entry) %>% pull(1) %>% unique()
  # merge subset specific and LUX subset specific IDs (to later on see interactions between LUX and prior knowledge sets)
  changes <- unique(c(changes, data_prot_diff_v31_LUX_sigup))
}# ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

#######################################################################################################################################################################################
#######################################################################################################################################################################################
for (loop_var1 in pageRank_dampening_values) {
  for (loop_var2 in pageRank_cutoff_values) {
    print(paste0("pageRank_dampening: ", loop_var1, "    |    pageRank_cutoff: ",  loop_var2))
    setwd(base_directory)
    result_directory = paste0(base_directory, "/",  script_version, "_", ppi_network, "-", loop_var1, "-", loop_var2, "-", param, "/")
    dir.create(result_directory, showWarnings = FALSE)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # source all functions
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    files.sources = list.files("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_functions", full.names=TRUE, pattern="*.R$")   # previously /Users/mgesell/Desktop/currentR/Network_plotting_Cathy/functions
    sapply(files.sources, source)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 0) house keeping: create empty folders for output
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # delete_existing_files   <- TRUE
    # initialize_folders(delete_existing_files)
    dir.create(paste0(result_directory, "final")       , showWarnings = FALSE)
    dir.create(paste0(result_directory, "intermediate"), showWarnings = FALSE)
    dir.create(paste0(result_directory, "plots")       , showWarnings = FALSE)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 1) load and filter input ppi_network
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" _____ section 1 _____________________________________________________________________________________________________________________________________________________")
    if (ppi_network == "string" ){             # ppi_network = string
      string_network        <- "physical"   # "physical" of "full"
      string_score_cutoff   <- 700 # cathy used 700 as default     string combined score ranges 0-1000; 400 threshold for medium confidence, 700 for high confidence;    low throughput high accuracy experiment data results score usually +/- 600; high throughput low accuracy data scores <= 250
      # string_linkage_categs <- c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining")   # select info levels to be used these are available:   c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining") 
      network_output_filename <- "network_interactions_string"   # 
      
      filter_string_network(string_network, string_score_cutoff, 
                            changes, 
                            network_output_filename)
      
    } else if (ppi_network =="complexPortal"){ # ppi_network = complexPortal
      complexPortal_input_file         <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/complex_portal_9606.tsv" #" /Users/mgesell/Desktop/currentR/Network_plotting_Cathy/input/complex_portal_9606.tsv"
      complexPortal_mapping_table_file <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/uniprotkb_AND_reviewed_true_AND_model_o_2024_11_08.tsv" #"/Users/mgesell/Desktop/currentR/Network_plotting_Cathy/input/uniprotkb_AND_reviewed_true_AND_model_o_2024_11_08.tsv"
      network_output_filename_biogrid   <- "network_interactions_complexPortal"   # 
      
      filter_complexPortal_network(complexPortal_input_file, complexPortal_mapping_table_file,
                                   changes,
                                   network_output_filename)
    } else if (ppi_network == "biogrid_physical") {
      biogrid_input_file              <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/BIOGRID-ALL-4.4.242_20250218.mitab.txt" #"/Users/mgesell/PhD/local_resources/PPIs_and_complexes/BIOGRID-ALL-4.4.242_20250218.mitab.txt"
      network_output_filename_cp  <- "network_interactions_biogrid_physical"   # 
     
       filter_biogrid_network(biogrid_input_file, 
                              changes, 
                              network_output_filename)
    } else if (ppi_network == "meta_pS700_CP_pBG") {
      network_output_filename      <- "network_interactions_meta_pS700_CP_pBG"   # 
      build_meta_network(network_output_filename)
    }
    
    extract_changes_internal_network(changes, network_output_filename)
    
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 2) create i-graph objects
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" _____ section 2 _____________________________________________________________________________________________________________________________________________________")
    #   settings
    p_cutoff          <- 0.05
    filter_column     <- "q"   # ideally set to to q once final data is here
    add_non_connected <- TRUE
    
    network_edgelist_file <- paste0(result_directory, "intermediate/", network_output_filename, ".tsv")
    # run the function
    create_igraph_objects(network_edgelist_file=network_edgelist_file, changes = changes)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 3) propagate the ppi_network
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" _____ section 3 _____________________________________________________________________________________________________________________________________________________")
    dampening_factor           <- loop_var1   # Keeping the damping factor around 0.85, which is commonly used and provides a balance between computation efficiency and ranking effectiveness
    pageRank_percentile_cutoff <- loop_var2   # string is messy --> be more selective
  
    # run the function
    run_pageRank(dampening_factor, pageRank_percentile_cutoff, pageRank_retain_all_changes)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 4) Clustering
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" _____ section 4 _____________________________________________________________________________________________________________________________________________________")
    #    settings walktrap
    walktrap_step_number         <- "auto"  # number or "auto" # cathy default = 6   4-5 generally recommended for ppi networks
    #    settings markov
    markov_expansion_parameter   <- 2
    markov_inflation_coefficient <- 2
    #    settings betweenness
    
    # run the function
    if (cluster_method == "walktrap"){
      walktrap_clustering(walktrap_step_number)
    } else if (cluster_method == "markov"     ){
      markov_clustering(markov_expansion_parameter, markov_inflation_coefficient)
    } else if (cluster_method == "betweenness"){   # slow
      betweenness_clustering()
    }
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 5) Check cluster for enrichment of original hits
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" _____ section 5 _____________________________________________________________________________________________________________________________________________________")
    #    settings 
    if (ppi_network == "string" ){                 #  string
      cluster_min_size = 2
      cluster_max_size = param 
      cluster_min_changes_above_2_subunits <- 1 
      cluster_min_changes_with_2_subunits  <- 1
    } else if (ppi_network =="complexPortal"){     # complexPortal
      cluster_min_size = 2
      cluster_max_size = param 
      cluster_min_changes_above_2_subunits <- 1 
      cluster_min_changes_with_2_subunits  <- 1
    } else if (ppi_network =="biogrid_physical"){  # biogrid_physical
      cluster_min_size = 2
      cluster_max_size = param 
      cluster_min_changes_above_2_subunits <- 1 
      cluster_min_changes_with_2_subunits  <- 1
    } else if (ppi_network =="meta_pS700_CP_pBG"){  # meta_pS700_CP_pBG
      cluster_min_size = 2
      cluster_max_size = param 
      cluster_min_changes_above_2_subunits <- 1 
      cluster_min_changes_with_2_subunits  <- 1
    }
    
    ###
    cluster_col_pal        <- "Pastel1"
    cluster_q_value_cutoff <- 0.05
    complex_portal2        <- fread("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/complex_portal_9606.tsv")   # "/Users/mgesell/Desktop/currentR/Network_plotting_Cathy/input/complex_portal_9606.tsv"    read protein complex names from the complex portal database
    complex_portal         <- fread("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_resources/unpacked_tweeked_complex_portal_mg/ComplexPortal_human9606_20250218.tsv") # "/Users/mgesell/Desktop/currentR/Network_plotting_Cathy/input/unpacked_tweeked_complex_portal_mg/ComplexPortal_human9606_20250218.tsv"  unboxed complexes by mg
    # run the function
    cluster_analysis(cluster_method, cluster_min_size, cluster_max_size, cluster_min_changes_above_2_subunits, cluster_min_changes_with_2_subunits,
                     cluster_col_pal, cluster_q_value_cutoff, complex_portal = complex_portal)# , 
                     # gene_names, string_input)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 6) Filter reference network for pageRank list --> node list (with interactors in second column) for import into cytoscape
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # run the function
    print(" _____ section 6 _____________________________________________________________________________________________________________________________________________________")
    extract_resource_subnetwork(ppi_network, pageRank_retain_all_changes)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 7) Global statistics
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # run the function
    print(" _____ section 7 _____________________________________________________________________________________________________________________________________________________")
    summary_plots()
    save.image(file= paste0(result_directory, "final/", "settings_used.RData"))
  
  }
}
