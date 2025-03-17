###  Hdownstream processing of shs script
###  2025-02-20
###  Martin Gesell
####################################################################################################################################
# 1.4  v31 csc lux hybrid script
# 1.5  v16.2 Jurkat act timecourse heatmap and 
# 2.0  include GO, interaction network data (biogrid, string) and complex data (complex_portal, compact)

rm(list = ls())                # clear workspace
script_version = "_shs-ds2.0"  # version

# Load libraries  ----------------------------------------------------------------------------------------------------------
library(dplyr)      # Data wrangling
library(tidyr)      # General functions
library(tidyverse)  # same      
library(ggplot2)      # Plotting
library(viridis)
library(scales)       # for alpha function; otherwise alpha applies for fill and outline in suvo
library(ggrepel)      # Labelling of points in plots without overlaps
library(ggpubr)       # Multi panel plotting
library(ggthemes)     # Plotting theme
library(plotly)       # interactive plots
library(htmlwidgets)  # html plots
library(gridExtra)    # plotting
library(ggplotify)     # for as.grob() transformation of qc17 qc_sample_correlation
library(cowplot)       # required for protti heatmap to be savable as composite
library(UpSetR)     # for overlap assessment of upregulated protein lists accross conditions (suvo_candi --> sigupuni)
library(grid)       # UpSetR title
library(pheatmap)
library(ComplexHeatmap)  # heatmap
library(circlize)        # heatmap subsets
library(dendextend)      # heatmap subsets
library(gprofiler2)      # GO enrichment 
#devtools::install_github("jpquast/protti", dependencies = TRUE, ref = "developer") # download developer version. specifically for calculate_protein_abundance peptide option
library(protti)     # https://cran.r-project.org/web/packages/protti/index.html
library(qvalue)       # qvalue calculation - note it is not the same as adj_pvalue. pvalue is 
library(magrittr)     # piping
library(purrr)        # helps handling lists
library(ggnewscale)       # dual heatmap plot
library(biomaRt)
set.seed(1234)
#
## load R functions
sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/shs_ds_functions", pattern = "\\.R$", full.names = TRUE), source)
#

## script parameters
test_parameter    = ""
poi_list          = "LUX_targets.CPIs" # select one colnames(poi):  
string_network    = "uniprot_swissprot"
#
# directories

directory_input_CSC_16 = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v16.2_CSC_deam_7aa__Jurkat-Tact_timecourse/_output_1-2_ludo_adjp_0.7string_shs2.22")
# directory_input_CSC   = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v16.2_CSC_deam_7aa__Jurkat-Tact_timecourse/_output2025-01-30_1-2_ludo_adjp_0.7string_shs2.21")
directory_input_CSC_31           = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_CSC_FP20_deam_semi_7aa_4ss/_output_1-2_ludo_adjp_0.7string_exclude_nCD4-2_nCD8-2_shs2.22")
directory_input_CSC_31_meta      = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_CSC_FP20_deam_semi_7aa_4ss/_output_1-2_ludo_adjp_0.7string_exclude_nCD4-2_nCD8-2_meta_shs2.22")
# directory_input_LUX   = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.21")
directory_input_LUX_31           = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_uniqueImpute_shs2.22")
directory_input_LUX_31_meta      = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_meta_uniqueImpute_shs2.22")
directory_input_LUX_31_full_meta = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_uniqueImpute_full_meta_shs2.22")

directory_input_LUX_24     = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v24_FP20_HoxHoxox_semi_6aa__panT_LUX_TCR-CD4-CD8/_output_1-2_ludo_adjp_0.7string_shs2.22")

directory_root        = c("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_shs_ds_2.1")
directory_output      = paste0(directory_root, "/", script_version, test_parameter)
if (!dir.exists( directory_output)) {  dir.create(directory_output)  }
setwd(paste(directory_root))   
#
# condition_levels_CSC      = c("NAct", "Act_0", "Act_0_5h", "Act_1h", "Act_2h", "Act_4h", "Act_24h", "Act_72h") # experiment specific condition & order desired in plots (& dataframes)
condition_levels_CSC_16 = c("NAct", "Act_0", "Act_0_5h" , "Act_1h",  "Act_2h",  "Act_4h",  "Act_24h",  "Act_72h")
condition_levels_LUX_24         = c("TCR", "CD4", "CD8") # experiment specific condition & order desired in plots (& dataframes)
condition_levels_LUX_24_all     = c("TCR", "CD4", "CD8", "Iso") # experiment specific condition & order desired in plots (& dataframes)
#
condition_levels_CSC_31          = c("nCD4", "nnCD4", "nCD8", "nnCD8") # experiment specific condition & order desired in plots (& dataframes)
condition_levels_LUX_31          = c("nCD4_TCR", "nnCD4_TCR", "nCD8_TCR", "nnCD8_TCR") # experiment specific condition & order desired in plots (& dataframes)
condition_levels_LUX_31_all      = c("nCD4_TCR", "nnCD4_TCR", "nCD8_TCR", "nnCD8_TCR","nCD4_Iso", "nnCD4_Iso", "nCD8_Iso", "nnCD8_Iso") # experiment specific condition & order desired in plots (& dataframes)
#
condition_levels_CSC_31_meta      = c("nCD4", "nnCD4", "nCD8", "nnCD8", "CD4meta", "CD8meta", "nMeta", "nnMeta") # experiment specific condition & order desired in plots (& dataframes)
condition_levels_LUX_31_meta      = c("nCD4_TCR", "nnCD4_TCR", "nCD8_TCR", "nnCD8_TCR", "CD4meta_TCR", "CD8meta_TCR", "nMeta_TCR", "nnMeta_TCR") # experiment specific condition & order desired in plots (& dataframes)
condition_levels_LUX_31_meta_all  = c("nCD4_TCR", "nnCD4_TCR", "nCD8_TCR", "nnCD8_TCR","nCD4_Iso", "nnCD4_Iso", "nCD8_Iso", "nnCD8_Iso", "CD4meta_TCR", "CD8meta_TCR", "nMeta_TCR", "nnMeta_TCR", "CD4meta_Iso", "CD8meta_Iso", "nMeta_Iso", "nnMeta_Iso") # experiment specific condition & order desired in plots (& dataframes)
#
condition_levels_LUX_31_meta_full  = c("metaTCR", "metaIso") # experiment specific condition & order desired in plots (& dataframes)


## human string - parameters ===================================================================================================================================================================================================================================================================
# usually don't touch
string_linkage_categs = c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining")   # select info levels to be used these are available:   c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining") 
string_min_comb_score = 700 # string combined score ranges 0-1000; 400 threshold for medium confidence, 700 for high confidence;    low throughput high accuracy experiment data results score usually +/- 600; high throughput low accuracy data scores <= 250
string_network        = "uniprot_unreviewed"   # select uniprot_unreviewed (default*) or uniprot_swissprot           *data searched for upsp >> will only look at upsp protein1 anyhow. any non-upsp network info (protein2)  might be helpfull for further analysis so keep itin
if (string_network == "uniprot_unreviewed") {  # --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  string <- read.csv("~/PhD/local_resources/string/string_net_unrev.csv", header = TRUE)  # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_unrev.csv
} else if (string_network == "uniprot_swissprot") {
  string <- read.csv("~/PhD/local_resources/string/string_net_upsp.csv",  header = TRUE)  # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_upsp.csv
}
# apply string info filtering
string <- string %>% 
  dplyr::select(all_of(c(string_linkage_categs, c("protein1_entry", "protein2_entry", "protein1_entry_name", "protein2_entry_name", "reviewed_protein1", "meta_surfaceome_protein1", "combined_score")))) %>%  # filter out unwanted categories  ...
  filter(!if_all(any_of(string_linkage_categs), ~ . == 0 | is.na(.))) %>%      # ... then remove all rows that lost their evidence due to columsn no longer present
  filter(combined_score >= string_min_comb_score) %>% 
  mutate(CD4_LUX = case_when(protein1_entry_name %in% c("CD3D_HUMAN", "CD3E_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN", "CD4_HUMAN")                              ~ 1, TRUE ~ 0),
         CD8_LUX = case_when(protein1_entry_name %in% c("CD3D_HUMAN", "CD3E_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN", "CD8A_HUMAN", "CD8B_HUMAN", "CD8B2_HUMAN")               ~ 1, TRUE ~ 0),
         TCR_LUX = case_when(protein1_entry_name %in% c("CD3D_HUMAN", "CD3E_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN",  "CD4_HUMAN", "CD8A_HUMAN", "CD8B_HUMAN", "CD8B2_HUMAN") ~ 1, TRUE ~ 0)) %>%
  filter(CD4_LUX > 0 | CD8_LUX > 0 | TCR_LUX > 0) # reduce dataframe to interactors (otherwise huge)
sum(string$CD4_LUX) # fyi 
sum(string$CD8_LUX) # fyi 
sum(string$TCR_LUX) # fyi 
#
# read resources
poi_reference        = read.csv("/Users/mgesell/PhD/local_resources/POI_lists/POI_lists.csv"   , header = TRUE, sep = ",")  # colnames(poi)
poi                  = as.vector(poi_reference[[poi_list]][poi_reference[[poi_list]] != ""])   
protein_families     = read_protti("/Users/mgesell/PhD/local_resources/protein_families.csv"   , header = TRUE, sep = ",")
surface_annotations  = read_protti("/Users/mgesell/PhD/local_resources/surface.annotations.csv", header = TRUE, sep = ",")
# #ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl" ) 
# #useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "https://may2024.archive.ensembl.org")
# #saveRDS(ensembl, file = "/Users/mgesell/PhD/local_resources/ensembl_feb_2025_mart.rds")
ensembl <- readRDS("/Users/mgesell/PhD/local_resources/ensembl_feb_2025_mart.rds")
unloadNamespace("biomaRt")
#
## proteome file loaded and pcois annotated
proteome <- load_uniprot_and_annotate(proteome = read_protti("/Users/mgesell/PhD/local_resources/human_upsp_202501.csv"  , header = TRUE, sep = ","))

##############################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################
## load data 
# protein level data
data_CSC_prot_v16         <-   load_protein_data(directory_input = directory_input_CSC_16                 , condition_levels = condition_levels_CSC_16)

data_CSC_prot_v31         <-   load_protein_data(directory_input = directory_input_CSC_31                 , condition_levels = condition_levels_CSC_31)
data_CSC_prot_v31_meta    <-   rbind(data_CSC_prot_v31,
                                     load_protein_data(directory_input = directory_input_CSC_31_meta      , condition_levels = condition_levels_CSC_31_meta) ) 
data_LUX_prot_v31          <-  load_protein_data(directory_input = directory_input_LUX_31                 , condition_levels = condition_levels_LUX_31)
data_LUX_prot_v31_meta     <-  rbind(data_LUX_prot_v31,
                                     load_protein_data(directory_input = directory_input_LUX_31_meta      , condition_levels = condition_levels_LUX_31_meta)  )
data_LUX_prot_v31_meta_full <- rbind(data_LUX_prot_v31_meta,
                                     load_protein_data(directory_input = directory_input_LUX_31_full_meta , condition_levels = condition_levels_LUX_31_meta_full))

data_LUX_prot_v24   <- load_protein_data(directory_input = directory_input_LUX_24, condition_levels = condition_levels_LUX_24_all)

# diff abundance data
data_CSC_prot_diff_v31           <- load_protein_data_diff(directory_input       = directory_input_CSC_31           )
data_CSC_prot_diff_v31_meta      <- rbind(data_CSC_prot_diff_v31,
                                          load_protein_data_diff(directory_input = directory_input_CSC_31_meta      )  )
data_LUX_prot_diff_v31           <- load_protein_data_diff(directory_input       = directory_input_LUX_31           )
data_LUX_prot_diff_v31_meta      <- rbind(data_LUX_prot_diff_v31,
                                          load_protein_data_diff(directory_input = directory_input_LUX_31_meta      )  )
data_LUX_prot_diff_v31_meta_full <- rbind(data_LUX_prot_diff_v31_meta,
                                          load_protein_data_diff(directory_input = directory_input_LUX_31_full_meta )  )

data_LUX_prot_diff_v24            <- load_protein_data_diff(directory_input       = directory_input_LUX_24           )
data_LUX_prot_diff_v24_v31_total  <- rbind(data_LUX_prot_diff_v24,
                                           data_LUX_prot_diff_v31_meta_full)


################################################################################################################################################
# set plot theme for script (from Amanda K)
# my plot theme to fit most journal figure requirements:
# black axis ticks and borders;   
# 5-7pt or 6-8pt text for a 2-panel figure    
# Line width 0.5-1.5 pt
# Colour blind friendly colour scheme (would recommend defining colours in source functions and use throughout) e.g. no red and green used together
# Gray fills between 10-80%
# Figure resolution at least 300 dpi ( I like to export in pdf not svg or tiff if using ggsave so figure is small and vectorized)
plot_theme <- function() {
  theme_classic() %+replace%
    theme(
      axis.text    = element_text(size = 16, face = "plain"),
      axis.text.x  = element_text(),
      axis.title   = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(),
      title        = element_text(size = 16, face = "plain"),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16),
      legend.text  = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin  = unit(c(t = 1, r = 1, b = 0.1, l = 1), "cm")
    )
}
theme_set(plot_theme()) # all plots generated from this script thosuld have same theme now 

sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/shs_ds_functions", pattern = "\\.R$", full.names = TRUE), source)
###################################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################################
################################################## End of Input Section ###########################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################################################


# v31 specific - check abundance vs LUX sig_up and LUX TCR core ###########################################################################################################################################################################
data_CSC_prot_LUX_filtered <- data_CSC_prot_v31 %>% 
  filter(entry %in% unique(data_LUX_prot_v31$entry) & log2_median > 0) 
data_CSC_prot_LUX_filtered <- data_CSC_prot_LUX_filtered  %>%
  mutate(# all enriched hits
    TCR_LUX_sigup   = ifelse(entry %in% unique(data_LUX_prot_diff_v31 %>%                                                        filter(log2FC >1) %>% pull(entry) ), "sig_up", "non_sig"),
    # positive lUX FC and mentioned 4x (core enriched)
    TCR_LUX_core    = ifelse(entry %in% unique(data_LUX_prot_diff_v31 %>% group_by(entry) %>% filter(n() == 4) %>% ungroup() %>% filter(log2FC >1) %>% pull(entry) ), "sig_up_core", ""),
    LUX_signal_plot = if_else(!is.na(TCR_LUX_core) & TCR_LUX_core != "", TCR_LUX_core, TCR_LUX_sigup)
  )
# Create a boxplot of log2_median grouped by LUX_signal_plot
ggplot(data_CSC_prot_LUX_filtered, aes(x = LUX_signal_plot, y = log2_median, fill = LUX_signal_plot)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("sig_up" = "#0570b0", "non_sig" = "grey", "sig_up_core" = "#31a354")) +
  labs(
    title = "v31 CSC abundances grouped by v31 LUX enrichment categories",
    x = "LUX Signal Plot",
    y = "log2 Median"
  ) +
  theme(legend.position = "none")
##########################################################################################################################################################################################################################################



### GO analysis #####################################################################################################################################################################################
## LUX -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ensembl id translation
library(biomaRt)
gost_overall_LUX <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full                                                        %>% pull(entry)), set = "overall" , max_term_size = 100)
gost_nCD4        <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD4_TCR_vs_nCD4_Iso")       %>% pull(entry)), set = "nCD4"    , max_term_size = 100)
gost_nnCD4       <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD4_TCR_vs_nnCD4_Iso")     %>% pull(entry)), set = "nnCD4"   , max_term_size = 100)
gost_nCD8        <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD8_TCR_vs_nCD8_Iso")       %>% pull(entry)), set = "nCD8"    , max_term_size = 100)
gost_nnCD8       <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD8_TCR_vs_nnCD8_Iso")     %>% pull(entry)), set = "nnCD8"   , max_term_size = 100)
gost_CD4meta       <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD4meta_TCR_vs_CD4meta_Iso") %>% pull(entry)), set = "CD4meta" , max_term_size = 100)
gost_CD8meta       <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD8meta_TCR_vs_CD8meta_Iso") %>% pull(entry)), set = "CD8meta" , max_term_size = 100)
gost_nMeta         <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nMeta_TCR_vs_nMeta_Iso")     %>% pull(entry)), set = "nMeta"   , max_term_size = 100)
gost_nnMeta        <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnMeta_TCR_vs_nnMeta_Iso")   %>% pull(entry)), set = "nnMeta"  , max_term_size = 100)
gost_TCRmeta       <- go_gost(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "metaTCR_vs_metaIso")         %>% pull(entry)), set = "metaTCR" , max_term_size = 100)
gost_panT_TCRLUX      <- go_gost(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "TCR_vs_Iso")         %>% pull(entry)), set = "panT_TCRLUX" , max_term_size = 100)
gost_panT_CD4LUX      <- go_gost(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD4_vs_Iso")         %>% pull(entry)), set = "panT_CD4LUX" , max_term_size = 100)
gost_panT_CD8LUX      <- go_gost(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD8_vs_Iso")         %>% pull(entry)), set = "panT_CD8LUX" , max_term_size = 100)

# combine GO sets to one long df
gost_LUX <- rbind(gost_overall_LUX, gost_nCD4, gost_nnCD4, gost_nCD8, gost_nnCD8, 
                  gost_CD4meta, gost_CD8meta, gost_nMeta, gost_nnMeta, gost_TCRmeta, 
                  gost_panT_TCRLUX, gost_panT_CD4LUX, gost_panT_CD8LUX) %>%
  group_by(term_name) %>%
  mutate(overlap = paste(sort(unique(comparison)), collapse = "_")) %>%
  ungroup()
# plotting ------
common_params <- list( # Define common plot parameters
  data = gost_LUX,  min_recall = 0.5,    min_p_value = 0.05,    grouping = "comparison",    term_column = "term_name",    distance_column = "recall",    distance_method = "euclidean", # distance parameters
  x_var = "comparison",    y_var = "term_name",    size_var = "recall",    fill_var = "p_value",    title_var = "GO Term Enrichment" # plot parameters
)
group_filters <- list( # Define different group filters
  NULL,
  c("overall", "metaTCR"),
  c("nCD4"   , "nnCD4"      , "nCD8"   , "nnCD8"),
  c("CD4meta", "CD8meta"    , "nMeta"  , "nnMeta"     , "metaTCR"),
  c("metaTCR", "panT_TCRLUX", "CD4meta", "panT_CD4LUX", "CD8meta", "panT_CD8LUX")
)
# Loop through group filters and create plots
plots_go_LUX <- lapply(group_filters, function(filter) {
  plot <- do.call(plot_dual_distance_bubble, c(common_params, list(group_filter = filter)))
  if (is.null(filter)) {     plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))   }
  print(plot)  # display plot
  return(plot) # store plot 
})

## CSC  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gost_overall_CSC        <- go_gost(query_list = unique(data_CSC_prot_diff_v31_meta                                                 %>% pull(entry)), set = "overall"            , max_term_size = 100)
gost_nCD4_vs_nCD8       <- go_gost(query_list = unique(data_CSC_prot_diff_v31_meta %>% filter(comparison == "nCD4_vs_nCD8"  )      %>% pull(entry)), set = "nCD4_vs_nCD8"       , max_term_size = 100)
gost_nCD4_vs_nnCD4      <- go_gost(query_list = unique(data_CSC_prot_diff_v31_meta %>% filter(comparison == "nCD4_vs_nnCD4" )      %>% pull(entry)), set = "nCD4_vs_nnCD4"      , max_term_size = 100)
gost_nnCD4_vs_nnCD8     <- go_gost(query_list = unique(data_CSC_prot_diff_v31_meta %>% filter(comparison == "nnCD4_vs_nnCD8")      %>% pull(entry)), set = "nnCD4_vs_nnCD8"     , max_term_size = 100)
gost_nCD8_vs_nnCD8      <- go_gost(query_list = unique(data_CSC_prot_diff_v31_meta %>% filter(comparison == "nCD8_vs_nnCD8" )      %>% pull(entry)), set = "nCD8_vs_nnCD8"      , max_term_size = 100)
gost_nMeta_vs_nnMeta      <- go_gost(query_list = unique(data_CSC_prot_diff_v31_meta %>% filter(comparison == "nMeta_vs_nnMeta" )    %>% pull(entry)), set = "nMeta_vs_nnMeta"    , max_term_size = 100)
gost_CD4meta_vs_CD8meta   <- go_gost(query_list = unique(data_CSC_prot_diff_v31_meta %>% filter(comparison == "CD4meta_vs_CD8meta" ) %>% pull(entry)), set = "CD4meta_vs_CD8meta" , max_term_size = 100)
# combine GO sets to one long df
gost_CSC <- rbind(gost_overall_CSC, gost_nCD4_vs_nCD8, gost_nCD4_vs_nnCD4, gost_nnCD4_vs_nnCD8, gost_nCD8_vs_nnCD8, gost_nMeta_vs_nnMeta, gost_CD4meta_vs_CD8meta) %>%
  group_by(term_name) %>%
  mutate(overlap = paste(sort(unique(comparison)), collapse = "_")) %>%
  ungroup()
# plotting ------
common_params <- list( # Define common plot parameters
  data = gost_CSC,
  min_recall = 0.5,  min_p_value = 0.05,  grouping = "comparison",  term_column = "term_name", distance_column = "recall",  distance_method = "euclidean", # distance parameters
  x_var = "comparison",  y_var = "term_name",  size_var = "recall",  fill_var = "p_value",  title_var = "GO Term Enrichment"                               # distance parameters
)
group_filters <- list( # Define different group filters
  NULL,
  c("overall"),
  c("CD4meta_vs_CD8meta", "nMeta_vs_nnMeta"),
  c("nCD4_vs_nnCD4"     , "nCD8_vs_nnCD8", "nCD4_vs_nCD8", "nnCD4_vs_nnCD8")
)
# Loop through group filters and create plots
plots_go_CSC <- lapply(group_filters, function(filter) {
  plot <- do.call(plot_dual_distance_bubble, c(common_params, list(group_filter = filter)))
  if (is.null(filter)) { plot <- plot + theme(axis.text.x = element_text(angle = 25, hjust = 1))  }
  print(plot)
  return(plot)
})
# END gost ________________________________________________________________________________________________________________________________________________________________________

## PPIs #########################################################################################################################################################################
tcr_chains <- poi_reference %>% filter(!tcr_chains_manual_entry == "") %>% pull(tcr_chains_manual_entry)
# biogrid PPIs -----------------------------------------------------------------------------------
ppi_bg_overall_LUX <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% pull(entry))                                                     , mode = "physical", set = "overall")
ppi_bg_nCD4        <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD4_TCR_vs_nCD4_Iso"    ) %>% pull(entry)), mode = "physical", set = "nCD4"   )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_nnCD4       <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD4_TCR_vs_nnCD4_Iso"  ) %>% pull(entry)), mode = "physical", set = "nnCD4"  )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_nnCD8       <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD8_TCR_vs_nnCD8_Iso"  ) %>% pull(entry)), mode = "physical", set = "nnCD8"  )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_nCD8        <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD8_TCR_vs_nCD8_Iso"    ) %>% pull(entry)), mode = "physical", set = "nCD8"   )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_nnMeta     <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnMeta_TCR_vs_nnMeta_Iso"  ) %>% pull(entry)), mode = "physical", set = "nnMeta" )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_nMeta      <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nMeta_TCR_vs_nMeta_Iso"    ) %>% pull(entry)), mode = "physical", set = "nMeta"  )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_CD4meta    <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD4meta_TCR_vs_CD4meta_Iso") %>% pull(entry)), mode = "physical", set = "CD4meta")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_CD8meta    <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD8meta_TCR_vs_CD8meta_Iso") %>% pull(entry)), mode = "physical", set = "CD8meta")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_TCRmeta    <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "metaTCR_vs_metaIso")         %>% pull(entry)), mode = "physical", set = "metaTCR")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_panT_TCRLUX  <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "TCR_vs_Iso")         %>% pull(entry)), mode = "physical", set = "panT_TCRLUX")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_panT_CD4LUX  <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD4_vs_Iso")         %>% pull(entry)), mode = "physical", set = "panT_CD4LUX")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_bg_panT_CD8LUX  <- ppi_biogrid(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD8_vs_Iso")         %>% pull(entry)), mode = "physical", set = "panT_CD8LUX")  # mode is dummy parameter for now - considering to update later to allow non-physical search too

# combine bg ppi sets to one long df
ppi_bg <- rbind(ppi_bg_overall_LUX, ppi_bg_nCD4, ppi_bg_nnCD4, ppi_bg_nnCD8, ppi_bg_nCD8, 
                ppi_bg_nnMeta, ppi_bg_nMeta, ppi_bg_CD4meta, ppi_bg_CD8meta, ppi_bg_TCRmeta,
                ppi_bg_panT_TCRLUX, ppi_bg_panT_CD4LUX, ppi_bg_panT_CD8LUX) %>%
  group_by(node_bg) %>%
  mutate(overlap = paste(sort(unique(comparison)), collapse = "_")) %>%
  ungroup()
#
ppi_bg %>%
  ggplot(aes(x = comparison, y = node_degree_biogrid_global)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "Comparison", y = "Node Degree (BioGrid)", title = "BioGrid: Global Degree of Query Proteins") +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
               hjust = 0.5, vjust = -0.5, color = "black", size = 5) 
ppi_bg %>%
  ggplot(aes(x = comparison, y = node_degree_biogrid_query_subset)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "Comparison", y = "Node Degree (BioGrid)", title = "BioGrid: Query Internal Degree") +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
               hjust = 0, vjust = -0.5, color = "black", size = 5) 
# not sure how to continue analysis from there
# ???

# string PPIs -----------------------------------------------------------------------------------
ppi_str_overall_LUX <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>%                                                    pull(entry)), mode = "physical", set = "overall" , string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_nCD4        <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD4_TCR_vs_nCD4_Iso"  ) %>% pull(entry)), mode = "physical", set = "nCD4"    , string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_nnCD4       <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD4_TCR_vs_nnCD4_Iso") %>% pull(entry)), mode = "physical", set = "nnCD4"   , string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_nnCD8       <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD8_TCR_vs_nnCD8_Iso") %>% pull(entry)), mode = "physical", set = "nnCD8"   , string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_nCD8        <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD8_TCR_vs_nCD8_Iso"  ) %>% pull(entry)), mode = "physical", set = "nCD8"    , string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_nnMeta    <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnMeta_TCR_vs_nnMeta_Iso"  ) %>% pull(entry)), mode = "physical", set = "nnMeta" , string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_nMeta     <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nMeta_TCR_vs_nMeta_Iso"    ) %>% pull(entry)), mode = "physical", set = "nMeta"  , string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_CD4meta   <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD4meta_TCR_vs_CD4meta_Iso") %>% pull(entry)), mode = "physical", set = "CD4meta", string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_CD8meta   <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD8meta_TCR_vs_CD8meta_Iso") %>% pull(entry)), mode = "physical", set = "CD8meta", string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_metaTCR   <- ppi_string(query_list = unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "metaTCR_vs_metaIso")         %>% pull(entry)), mode = "physical", set = "metaTCR", string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_panT_TCRLUX  <- ppi_string(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "TCR_vs_Iso")         %>% pull(entry)), mode = "physical", set = "panT_TCRLUX", string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_panT_CD4LUX <- ppi_string(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD4_vs_Iso")         %>% pull(entry)), mode = "physical", set = "panT_CD4LUX", string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
ppi_str_panT_CD8LUX <- ppi_string(query_list = unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD8_vs_Iso")         %>% pull(entry)), mode = "physical", set = "panT_CD8LUX", string_phys_confidence = "medium")  # mode is dummy parameter for now - considering to update later to allow non-physical search too

# combine bg ppi sets to one long df
ppi_str <- rbind(ppi_str_overall_LUX, ppi_str_nCD4, ppi_str_nnCD4, ppi_str_nnCD8, ppi_str_nCD8, ppi_str_nnMeta, 
                 ppi_str_nMeta, ppi_str_CD4meta, ppi_str_CD8meta, ppi_str_metaTCR,
                 ppi_str_panT_TCRLUX, ppi_str_panT_CD4LUX, ppi_str_panT_CD8LUX) %>%
  group_by(node_str) %>%
  mutate(overlap = paste(sort(unique(comparison)), collapse = "_")) %>%
  ungroup()

ppi_str %>%
  ggplot(aes(x = comparison, y = node_degree_string_phys_global)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "Comparison", y = "Node Degree (STRING)", title = "STRING: Global Physical Degree of Query Proteins") +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
               hjust = 0.5, vjust = -0.5, color = "black", size = 5) 
ppi_str %>%
  ggplot(aes(x = comparison, y = node_degree_string_phys_query_subset)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "Comparison", y = "Node Degree (STRING)", title = "STRING: Query Internal Physical Degree") +# not sure how to continue analysis from there
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
               hjust = 0, vjust = -0.5, color = "black", size = 5) 
# ???
## PPIs end _____________________________________________________________________________________________________________________________________________________________________

################################################################################################################################################################################################################################
################################################################################################################################################################################################################################
## Complexes ###################################################################################################################################################################################################################
sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/shs_ds_functions", pattern = "\\.R$", full.names = TRUE), source)
# complex_protal complexes ------------------------------------------------------------------------------------------------------------------------------------------------------------
comp_cp_overall_LUX <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full                                                      %>% pull(entry)) ), set = "overall")  
comp_cp_nCD4        <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD4_TCR_vs_nCD4_Iso"  )   %>% pull(entry)) ), set = "nCD4" )  
comp_cp_nnCD4       <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD4_TCR_vs_nnCD4_Iso")   %>% pull(entry)) ), set = "nnCD4")  
comp_cp_nCD8        <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD8_TCR_vs_nCD8_Iso"  )   %>% pull(entry)) ), set = "nCD8" )  
comp_cp_nnCD8       <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD8_TCR_vs_nnCD8_Iso")   %>% pull(entry)) ), set = "nnCD8")  
comp_cp_nMeta     <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nMeta_TCR_vs_nMeta_Iso"    ) %>% pull(entry)) ), set = "nMeta"  )  
comp_cp_nnMeta    <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnMeta_TCR_vs_nnMeta_Iso"  ) %>% pull(entry)) ), set = "nnMeta" )  
comp_cp_CD4meta   <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD4meta_TCR_vs_CD4meta_Iso") %>% pull(entry)) ), set = "CD4meta")  
comp_cp_CD8meta   <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD8meta_TCR_vs_CD8meta_Iso") %>% pull(entry)) ), set = "CD8meta")  
comp_cp_CD8meta   <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD8meta_TCR_vs_CD8meta_Iso") %>% pull(entry)) ), set = "CD8meta")  
comp_cp_metaTCR   <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "metaTCR_vs_metaIso")         %>% pull(entry)) ), set = "metaTCR")  
comp_cp_panT_TCRLUX    <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "TCR_vs_Iso")         %>% pull(entry)) ), set = "panT_TCRLUX")  
comp_cp_panT_CD4LUX   <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD4_vs_Iso")         %>% pull(entry)) ), set = "panT_CD4LUX")  
comp_cp_panT_CD8LUX   <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD8_vs_Iso")         %>% pull(entry)) ), set = "panT_CD8LUX")  
# merge df
comp_cp_LUX <- rbind(comp_cp_overall_LUX, comp_cp_nCD4, comp_cp_nnCD4, comp_cp_nnCD8, comp_cp_nCD8, 
                     comp_cp_nnMeta, comp_cp_nMeta, comp_cp_CD4meta, comp_cp_CD8meta, comp_cp_metaTCR,
                     comp_cp_panT_TCRLUX, comp_cp_panT_CD4LUX, comp_cp_panT_CD8LUX) %>%
  group_by(complex_name_cp) %>%
  mutate(overlap = paste(sort(unique(comparison)), collapse = "_")) %>%
  ungroup() %>%
  arrange(desc(complex_recall_cp))
# plotting Complex Portal ------
common_params <- list( # Define common plot parameters
  data = comp_cp_LUX,  min_recall = 0.5,    min_p_value = NULL,    grouping = "comparison",    term_column = "recommended_name",    distance_column = "complex_recall_cp",    distance_method = "euclidean", # distance parameters
  x_var = "comparison",    y_var = "recommended_name",    size_var = "complex_recall_cp",  fill_var = NULL,    title_var = "Complex Portal" # plot parameters
)
group_filters <- list( # Define different group filters
  NULL,
  c("overall", "metaTCR"),
  c("nCD4"   , "nnCD4"      , "nCD8"   , "nnCD8"),
  c("CD4meta", "CD8meta"    , "nMeta"  , "nnMeta"     , "metaTCR"),
  c("metaTCR", "panT_TCRLUX", "CD4meta", "panT_CD4LUX", "CD8meta", "panT_CD8LUX")
)
plots_cp <- lapply(group_filters, function(filter) { # Loop through group filters and create plots
  plot <- do.call(plot_dual_distance_bubble, c(common_params, list(group_filter = filter)))
  if (is.null(filter)) {     plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))   }   # Add custom theme for NULL filter
  print(plot) # Display the plot in RStudio's Plots window
  return(plot) # Store the plot in a list for further use
})

# corum complexes --------------------------------------------------------------------------------------------------------------------------------------------------------------------
comp_corum_overall_LUX <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full                                                      %>% pull(entry)) ), set = "overall")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_nCD4        <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD4_TCR_vs_nCD4_Iso"  )   %>% pull(entry)) ), set = "nCD4" )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_nnCD4       <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD4_TCR_vs_nnCD4_Iso")   %>% pull(entry)) ), set = "nnCD4")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_nCD8        <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nCD8_TCR_vs_nCD8_Iso"  )   %>% pull(entry)) ), set = "nCD8" )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_nnCD8       <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnCD8_TCR_vs_nnCD8_Iso")   %>% pull(entry)) ), set = "nnCD8")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_nMeta     <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nMeta_TCR_vs_nMeta_Iso"    ) %>% pull(entry)) ), set = "nMeta"  )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_nnMeta    <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "nnMeta_TCR_vs_nnMeta_Iso"  ) %>% pull(entry)) ), set = "nnMeta" )  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_CD4meta   <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD4meta_TCR_vs_CD4meta_Iso") %>% pull(entry)) ), set = "CD4meta")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_CD8meta   <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "CD8meta_TCR_vs_CD8meta_Iso") %>% pull(entry)) ), set = "CD8meta")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_metaTCR   <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full %>% filter(comparison == "metaTCR_vs_metaIso")         %>% pull(entry)) ), set = "metaTCR")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_panT_TCRLUX  <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "TCR_vs_Iso")         %>% pull(entry)) ), set = "panT_TCRLUX")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_panT_CD4LUX  <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD4_vs_Iso")         %>% pull(entry)) ), set = "panT_CD4LUX")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
comp_corum_panT_CD8LUX  <- complex_corum(query_list = sort(unique(data_LUX_prot_diff_v24_v31_total %>% filter(comparison == "CD8_vs_Iso")         %>% pull(entry)) ), set = "panT_CD8LUX")  # mode is dummy parameter for now - considering to update later to allow non-physical search too
# merge df
comp_cor <- rbind(comp_corum_overall_LUX, comp_corum_nCD4, comp_corum_nnCD4, comp_corum_nCD8, comp_corum_nnCD8, 
                  comp_corum_nMeta, comp_corum_nnMeta, comp_corum_CD4meta, comp_corum_CD8meta, comp_corum_metaTCR,
                  comp_corum_panT_TCRLUX, comp_corum_panT_CD4LUX, comp_corum_panT_CD8LUX) %>%
  group_by(complex_name_cor) %>%
  mutate(overlap = paste(sort(unique(comparison)), collapse = "_")) %>%
  ungroup() %>%
  arrange(desc(complex_recall_cor))

# plotting Corum Complexes ------
common_params <- list( # Define common plot parameters
  data = comp_cor,  min_recall = 0.5,    min_p_value = NULL,    grouping = "comparison",    term_column = "complex_name_cor",    distance_column = "complex_recall_cor",    distance_method = "euclidean", # distance parameters
  x_var = "comparison",    y_var = "complex_name_cor",    size_var = "complex_recall_cor",    fill_var = NULL,    title_var = "Corum Complexes" # plot parameters
)
group_filters <- list( # Define different group filters
  NULL,
  c("overall", "metaTCR"),
  c("nCD4"   , "nnCD4"      , "nCD8"   , "nnCD8"),
  c("CD4meta", "CD8meta"    , "nMeta"  , "nnMeta"     , "metaTCR"),
  c("metaTCR", "panT_TCRLUX", "CD4meta", "panT_CD4LUX", "CD8meta", "panT_CD8LUX")
)
# Loop through group filters and create plots
plots_cor <- lapply(group_filters, function(filter) {
  plot <- do.call(plot_dual_distance_bubble, c(common_params, list(group_filter = filter)))
  # Add custom theme for NULL filter
  if (is.null(filter)) {     plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))   }
  print(plot) # Display the plot in RStudio's Plots window
  return(plot) # Store the plot in a list for further use
})



























































# ===========================================================================================================================================================================
# color palettes (define after imputation!)
# colors_csc_zscore <- colorRamp2(c(NA_impute_CSC_zscore, max(data_CSC_prot_v31$log2_z_score, na.rm = TRUE)),
# c("white","firebrick"))  
colors_csc_log2median <- colorRamp2(c(min(data_CSC_prot_v31$log2_median, na.rm = TRUE), max(data_CSC_prot_v31$log2_median, na.rm = TRUE)),
                                    c("white","firebrick"))  




##################################################################################################################################################
## pcoi heatmaps #################################################################################################################################

cluster_based_on = "CSC_based_clustering" # "CSC_based_clustering" or "LUX_binary"
# loop for plotting PCOI heatmaps
for (pcoi in pcoi_list) {
  ## subset dataframep pbased on current CPOI  
  df <- data_merge_w_pcoi_tile_plot %>% 
    filter(!!sym(pcoi) == 1) %>% 
    dplyr::select(1:8) %>% # filter for pcoi of interest and subset for data columns (unselect pcoi columns)
    {rownames(.) <- gsub("_HUMAN", "", rownames(.)); .}  # for simpler plots 
  # Skip if matrix is empty
  if (nrow(df) <= 0) {
    print(paste("Skipping", pcoi, "- no data"))
    next
  }
  # Check for unique entries
  if (nrow(df) == 1) {
    # If there's only one unique entry, create a heatmap without clustering
    df_long <- df %>%
      rownames_to_column("entry_name") %>%
      pivot_longer(cols = all_of(column_order),
                   names_to = "condition",
                   values_to = "value") %>%
      mutate(dataset   = ifelse(grepl("LUX$", condition), "LUX", "CSC"),
             condition = factor(condition, levels = column_order),
             color     = case_when(
               dataset == "LUX" & value == 1 ~ "black",
               dataset == "LUX" & value == 0 ~ "white",
               dataset == "CSC" ~ colors_csc_zscore(value)
             ))
    # Create and save heatmap with clustering
    png(filename = paste0(directory_output, "/_heatmap_merge/_", pcoi, ".png"), width = 800, height = 200 + length(unique(df_long$entry_name))*50)
    print(
      ggplot(df_long, aes(x = condition, y = entry_name)) +
        geom_tile(aes(fill = color), color = "grey") +
        scale_fill_identity(guide = "none") +
        new_scale_fill() +
        geom_point(data = df_long[df_long$dataset == "LUX",], 
                   aes(fill = factor(value)), shape = 22, size = 3, stroke = 0) +
        scale_fill_manual(values = lux_colors, name = "LUX") +
        new_scale_color() +
        geom_point(data = df_long[df_long$dataset == "CSC",], 
                   aes(color = value), shape = 21, size = 3, stroke = 0) +
        scale_color_gradient(low = "white", high = "firebrick", name = "CSC", limits = range(data_CSC_prot_v31$log2_z_score, na.rm = TRUE)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              panel.grid = element_blank()) +
        labs(x = "Condition", y = "Entry_name", title = paste(pcoi))
    )
    dev.off()
    print(paste("Heatmap merge saved for:", pcoi))
  } else {
    ### Clustering ================================================================================
    if(cluster_based_on == "CSC_continous") {             # for continuous data
      # Extract LUX data and convert to matrix for clustering
      cluster_data <- as.matrix(df[, grep("LUX$", names(df))])
      ## continuous value clustering -------------------------
      # Perform hierarchical clustering
      dist_matrix <- dist(cluster_data, method = "euclidean")    # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
      hc          <- hclust(dist_matrix)   # For z-scores, Euclidean distance is commonly used. However, if your data has a lot of zeros or is sparse, consider using Manhattan distance or correlation-based distances (like Pearson or Spearman) to capture relationships more effectively.
      # Get the order of rows based on clustering
      row_order <- hc$order
    } else if (cluster_based_on == "LUX_binary") {     # for binary values
      # Extract LUX data and convert to matrix for clustering
      cluster_data <- as.matrix(df[, grep("LUX$", names(df))])
      ## binary clustering -----------------------------------
      # Perform hierarchical clustering using Jaccard distance (suitable for binary data)
      dist_matrix <- dist(cluster_data, method = "binary")    # Use binary distance
      hc          <- hclust(dist_matrix, method = "ward.D2")   # methods "ward.D2" (best) "average" (ok, looks like complete) or "complete" (ok, looks like average)
      # Get the order of rows based on clustering
      row_order <- hc$order
      # ================================================================================================
    }
    
    ## enforce the cluster order into dataframe  
    df_ordered <- df[row_order, ]
    # Prepare the data
    df_long <- df_ordered %>%
      rownames_to_column("entry_name") %>%
      pivot_longer(cols = all_of(column_order),
                   names_to = "condition",
                   values_to = "value") %>%
      mutate(dataset   = ifelse(grepl("LUX$", condition), "LUX", "CSC"),
             condition = factor(condition, levels = column_order),
             entry_name   = factor(entry_name, levels = unique(rownames(df_ordered))),  # Set entry_name levels based on clustering
             color     = case_when(
               dataset == "LUX" & value == 1 ~ "black",
               dataset == "LUX" & value == 0 ~ "white",
               dataset == "CSC" ~ colors_csc_zscore(value)
             ))
    # Create and save heatmap with clustering
    png(filename = paste0(directory_output, "/_heatmap_merge/_", pcoi, ".png"), width = 800, height = 200 + length(unique(df_long$entry_name))*10)
    print(
      ggplot(df_long, aes(x = condition, y = entry_name)) +
        geom_tile(aes(fill = color), color = "grey") +
        scale_fill_identity(guide = "none") +
        new_scale_fill() +
        geom_point(data = df_long[df_long$dataset == "LUX",], 
                   aes(fill = factor(value)), shape = 22, size = 3, stroke = 0) +
        scale_fill_manual(values = lux_colors, name = "LUX") +
        new_scale_color() +
        geom_point(data = df_long[df_long$dataset == "CSC",], 
                   aes(color = value), shape = 21, size = 3, stroke = 0) +
        scale_color_gradient(low = "white", high = "firebrick", name = "CSC", limits = range(data_CSC_prot_v31$log2_z_score, na.rm = TRUE)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              panel.grid = element_blank()) +
        labs(x = "Condition", y = "Entry_name", title = paste(pcoi))
    )
    dev.off()
    print(paste("Heatmap merge saved for:", pcoi))
  }
}

########################################################################################################
########################################################################################################
# HEATMAP LUX only 
if (!dir.exists( paste0(directory_output,"/_heatmap_lux") )) {  dir.create( paste0(directory_output,"/_heatmap_lux") )  }
# exlcude csc data
data_LUX_pcoi <- data_merge_w_pcoi[,-c(5:8)] 

# loop for plotting PCOI heatmaps
for (pcoi in pcoi_list) {
  # Subset and reshape data
  data_LUX_pcoi_ss <- data_LUX_pcoi %>% 
    filter(!!sym(pcoi) == "1") %>%
    dplyr::select(nCD4_LUX, nnCD4_LUX, nnCD8_LUX, nCD8_LUX) %>%
    mutate(row_sum = nCD4_LUX + nnCD4_LUX + nnCD8_LUX + nCD8_LUX) %>%  # filter row
    filter(row_sum > 0) %>%                                            # filter out NA values (introduced due to presence in CSC dataset)
    dplyr::select(-row_sum) %>%
    as.matrix()
  
  # Skip if matrix is empty
  if (nrow(data_LUX_pcoi_ss) == 0) {  # skip over zero line subsets
    print(paste("Skipping", pcoi, "- no data"))
    next
  }
  
  # Check for unique entries
  if (nrow(data_LUX_pcoi_ss) == 1) {
    # If there's only one unique entry, create a heatmap without clustering
    my_colors <- colorRampPalette(c("white", "black"))(2)
    
    png(filename = paste0(directory_output, "/_heatmap_lux/_", pcoi, ".png"), width = 800, height = 200 + nrow(data_LUX_pcoi_ss)*50)
    print(
      pheatmap(data_LUX_pcoi_ss,
               color = my_colors,
               cluster_rows = FALSE,  # No clustering
               cluster_cols = FALSE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               main = paste(pcoi),
               fontsize_row = 22,
               fontsize_col = 22,
               fontsize = 25,
               legend = FALSE)
    )
    dev.off()
    
    print(paste(pcoi, "without clustering"))
    
  } else {
    # Create color palette
    my_colors <- colorRampPalette(c("white", "black"))(2)
    
    # Create and save heatmap with clustering
    png(filename = paste0(directory_output, "/_heatmap_lux/_", pcoi, ".png"), width = 800, height = 200 + nrow(data_LUX_pcoi_ss)*50)
    print(
      pheatmap(data_LUX_pcoi_ss,
               color = my_colors,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               main = paste(pcoi),
               fontsize_row = 22,
               fontsize_col = 22,
               fontsize = 25,
               legend = FALSE)
    )
    dev.off()
    
    print(paste("LUX heatmap saved for:", pcoi))
    
  }
}
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# ================================================================================================================================================================
# ################################################################################################################################################################


# # heatmap plot
# ggplot(LUX_binary_heatmap_ss, aes(x = condition, y = entry_name, fill = factor(LUX_proximity))) +
#   geom_tile() +
#   scale_fill_manual(values = c("0" = "white", "1" = "black")) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         axis.text.y = element_text(size = 6)) +
#   labs(x = "Condition", y = "Entry Name", fill = "TCR Chain") +
#   coord_fixed(ratio = 1)



# ################################################################################################################################################################
# CSC line plots - median signal ###################################################################################################################################################
# Plotting global abundance (CSC data) of LUX enriched proteins across subsets ========================================================================================================================
data_csc_medianlog2_L <- data_csc_medianlog2_L %>%
  pivot_wider(id_cols = entry, names_from = condition, values_from = log2_median) %>%
  pivot_longer(cols = c(nCD4,  nCD8, nnCD8, nnCD4), names_to = "condition", values_to = "log2_median") %>%
  left_join(proteome %>% dplyr::select(-c("reviewed","protein_names", "gene_names", "organism", "length", "gene_ontology_go", "gene_ontology_biological_process", "gene_ontology_cellular_component", 
                                          "gene_ontology_molecular_function", "gene_ontology_i_ds", "subcellular_location_cc", "transmembrane", "function_cc", "pathway", "gene_names_primary", "disulfide_bond", 
                                          "glycosylation", "lipidation", "intramembrane", "topological_domain", "ensembl", "sequence", "interacts_with")),
            by = "entry") %>%
  mutate(entry_name = gsub("_HUMAN", "", entry_name)) %>%
  dplyr::select(sort(names(.))) %>%
  mutate(condition = factor(condition, levels = c("nCD4", "nnCD4", "nCD8", "nnCD8")))
# directory
if (!dir.exists( paste0(directory_output,"/_lineplot_csc" ) )) {  dir.create( paste0(directory_output,"/_lineplot_csc" ) )  }
# line plot for pcois lists --------------------------------------------------------
for (pcoi in pcoi_list) {
  # subest data for pcoi
  data_CSC_median_ss <- data_csc_medianlog2_L %>%
    filter(!!sym(pcoi) == 1) %>%
    dplyr::select(entry, entry_name, condition, log2_median) 
  # Skip if matrix is empty
  if (nrow(data_CSC_median_ss) == 0) {  # skip over zero line subsets
    print(paste("Skipping", pcoi, "- no data"))
    next
  } else {
    png(filename = paste0(directory_output, "/_lineplot_csc/_", pcoi, ".png"), width = 800, height = 200 + length(unique(data_CSC_median_ss$entry_name))*50)
    print(
      ggplot(data_CSC_median_ss, aes(x = condition, y = log2_median, color = entry_name, group = entry_name)) +
        geom_line(size = 1.1) +
        geom_point(size = 3, shape =18) +
        theme_minimal() +
        labs(x = "Condition", y = "Abundance", title = paste0("Median log2 Protein Abundances: ", pcoi) +
               theme(
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
                 axis.text.y = element_text(size = 15),
                 axis.title = element_text(size = 15),
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 12),
                 legend.key.size = unit(1.5, "lines")
               ) +
               scale_color_viridis(discrete = TRUE, option = "turbo") +
               ylim(min(data_csc_medianlog2_L$log2_median, na.rm = TRUE) - 0.1, 
                    max(data_csc_medianlog2_L$log2_median, na.rm = TRUE) + 0.1)
        )
    )
    dev.off()
    print(paste0("Plotted: ", pcoi))
    
  }
}
# line plot for custom protein list --------------------------------------------------------
protein_list <- poi_reference %>% filter(Thermo_Marker_Tact != "") %>% pull(Thermo_Marker_Tact)   #   protein_list <-  c("CD69_HUMAN", "IL2RA_HUMAN") #
data_CSC_median_ss <- data_CSC_prot_v16 %>% # data_CSC_prot_v16  data_CSC_prot_v31
  left_join(proteome %>% 
              dplyr::select(entry, entry_name), by = "entry") %>%
  filter(entry_name %in% protein_list) %>%
  dplyr::select(entry, entry_name, condition, log2_median)

# Plotting
ggplot(data_CSC_median_ss, aes(x = condition, y = log2_median, color = entry_name, group = entry_name)) +
  geom_line(size = 1.1) +
  geom_point(size = 3, shape = 18) +
  theme_minimal() +
  labs(x = "Condition", y = "Abundance", title = "Protein Abundance Across Conditions") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  ) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  ylim(15, 
       max(data_CSC_median_ss$log2_median, na.rm = TRUE) + 0.1)




## CSC UPSET PLOT  ==========================================================================================================================================================
# data formatting (matrix and binary)
data_CSC_prot_w_bin <- data_CSC_prot_w %>% 
  column_to_rownames("entry") %>%
  mutate(across(everything(), ~as.numeric(!is.na(.))))
# plot
png(paste0(directory_output, "/_upset_plot_CSC.png"), width = 800, height = 600)
upset(data_CSC_prot_w_bin , 
      sets = colnames(data_CSC_prot_w_bin), 
      order.by = "freq", 
      nsets = length(colnames(data_CSC_prot_w_bin)),
      nintersects = NA,
      text.scale = c(2, 2, 2, 2, 2, 1.5))  
grid.text("Sig-Up regulated Proteins", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))
dev.off() # Close the device













