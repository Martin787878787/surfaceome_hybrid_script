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
#
rm(list = ls())           # Purge workspace
set.seed(123)
#######################################################################################################################################################################################
#######################################################################################################################################################################################
script_version   = "NeA_1.3" # Network Analysis human 1.0
ppi_network      = "string"      # c("string", "complexPortal")
cluster_method   = "walktrap"    # c("walktrap", "markov", "betweenness")
## pageRank
pageRank_dampening_values = c(0.85)  # default 0.85
pageRank_cutoff_values    = c(0.5, 0.6, 0.7, 0.75, 0.8)  # the more "changes" proteins the higher number is recommended. (parameters defines which edges least connected to input nodes - higher --> more stringent cutting)
pageRank_retain_all_changes = FALSE   # chose based on your goal: TRUE~"keep all input nodes"; FALSE~"keep only nodes with highest connectivity 
## clustering
walktrap_step_number_param = "auto"    # "auto" (auto_optimize) or *number*
# defines max cluster size (all larger than this excluded for downstream analysis)
param = 50 # start high (e.g. 500) than check "number_of_proteins_per_cluster.pdf" for orientation & titrate down (e.g. 100, 30) until "clustered_network_significant_with_changes_labeled.pdf" looks nice.  

# directory
base_directory  =  paste0("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_NeA")
dir.create(paste0(base_directory), showWarnings = FALSE)
# input data
# changes         <- read.delim("/Users/mgesell/Desktop/currentR/Network_plotting_Cathy/v31_1-2_imputed_sigup-meta_nounique.txt", header = FALSE, stringsAsFactors = FALSE, sep = '\t') %>% pull(1)
# the file is/can be generated with shs_downstream_2.2 see cytoscape section
changes         <- read.delim("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/LUX_metachanges_up_shs2.22.csv", header = FALSE, stringsAsFactors = FALSE, sep = ',') %>% pull(1)

#######################################################################################################################################################################################
#######################################################################################################################################################################################

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# source all functions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
files.sources = list.files("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_functions", full.names=TRUE, pattern="*.R$")
sapply(files.sources, source)

for (loop_var1 in pageRank_dampening_values) {
  for (loop_var2 in pageRank_cutoff_values) {  # }}
    print(" ============= starting loop =============")
    print(paste0("pageRank_dampening: ", loop_var1, "    |    pageRank_cutoff: ",  loop_var2))
    setwd(base_directory)
    result_directory = paste0(base_directory, "/",  script_version, "_", ppi_network, "-", loop_var1, "-", loop_var2, "-", param, "-walktrap", walktrap_step_number_param, "/")
    dir.create(result_directory, showWarnings = FALSE)

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
    print(" ------------ ppi network ------------")
    if (ppi_network == "string" ){             # ppi_network = string
      string_network        <- "uniprot_swissprot"   # uniprot_unreviewed or uniprot_swissprot
      string_score_to_use   <- "combined_score"
      string_score_cutoff   <- 400 # cathy used 700 as default     string combined score ranges 0-1000; 400 threshold for medium confidence, 700 for high confidence;    low throughput high accuracy experiment data results score usually +/- 600; high throughput low accuracy data scores <= 250
      # string_linkage_categs <- c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining")   # select info levels to be used these are available:   c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining") 
      network_output_filename <- "network_interactions_string"   # 
      
      filter_string_network(string_network, string_score_to_use, string_score_cutoff, 
                            changes, 
                            network_output_filename)
      
    } else if (ppi_network =="complexPortal"){ # ppi_network = complexPortal
      complexPortal_input_file         <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_input/complex_portal_9606.tsv"
      complexPortal_mapping_table_file <- "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_input/uniprotkb_AND_reviewed_true_AND_model_o_2024_11_08.tsv"
      network_output_filename          <- "network_interactions_complexPortal"   # 
      ø
      filter_complexPortal_network(complexPortal_input_file, complexPortal_mapping_table_file,
                                   network_output_filename)
    }
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 2) create i-graph objects
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" ------------ i-graph ------------")
    #   settings
    p_cutoff          <- 0.05
    filter_column     <- "q"   # ideally set to to q once final data is here
    add_non_connected <- TRUE
    
    network_edgelist_file <- paste0(result_directory, "intermediate/", network_output_filename, ".tsv")
    # run the function
    create_igraph_objects(network_edgelist_file=network_edgelist_file, changes = changes, 
                          add_non_connected = FALSE) #
    
    load( file = paste0(result_directory, "intermediate/igraph_diagnostics.rds"))
    missing_in_igraph <- igraph_diagnostics$missing_in_igraph

    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 3) propagate the ppi_network
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" ------------ pageRank network propagation ------------")
    dampening_factor           <- loop_var1   # Keeping the damping factor around 0.85, which is commonly used and provides a balance between computation efficiency and ranking effectiveness
    pageRank_percentile_cutoff <- loop_var2   # string is messy --> be more selective
  
    # run the function
    run_pageRank(dampening_factor, pageRank_percentile_cutoff, pageRank_retain_all_changes)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 4) Clustering
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" ------------ network clustering ------------")
    #    settings walktrap
    walktrap_step_number         <- walktrap_step_number_param  # 6 was in cathys default function
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
    print(" ------------ relating network clustering to original hits ------------")
    #    settings 
    if (ppi_network == "string" ){             # ppi_network = string
      cluster_min_size = 2
      cluster_max_size = param 
      cluster_min_changes_above_2_subunits <- 1 
      cluster_min_changes_with_2_subunits  <- 1
    } else if (ppi_network =="complexPortal"){ # ppi_network = complexPortal
      cluster_min_size = 2
      cluster_max_size = param 
      cluster_min_changes_above_2_subunits <- 1 
      cluster_min_changes_with_2_subunits  <- 1
    }
    
    ###
    cluster_col_pal        <- "Pastel1"
    cluster_q_value_cutoff <- 0.05
    complex_portal2        <- fread("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_input/complex_portal_9606.tsv")   # read protein complex names from the complex portal database
    complex_portal         <- fread("/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/NeA_input/unpacked_tweeked_complex_portal_mg/ComplexPortal_human9606_20250218.tsv") # unboxed complexes by mg
    # run the function
    cluster_analysis(cluster_method, cluster_min_size, cluster_max_size, cluster_min_changes_above_2_subunits, cluster_min_changes_with_2_subunits,
                     cluster_col_pal, cluster_q_value_cutoff, complex_portal = complex_portal)# , 
                     # gene_names, string_input)
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # 6) Global statistics
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(" ------------ global stats ------------")
    # run the function
    summary_plots()
    save.image(file= paste0(result_directory, "final/", "settings_used.RData"))
  
  }
}
