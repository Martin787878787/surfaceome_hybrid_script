###                     __                                          
###                   / _|                                         
###    ___ _   _ _ __| |_ __ _  ___ ___  ___  _ __ ___   ___       
###   / __| | | | '__|  _/ _` |/ __/ _ \/ _ \| '_ ` _ \ / _ \      
###   \__ \ |_| | |  | || (_| | (_|  __/ (_) | | | | | |  __/      
###   |___/\__,_|_|_ |_| \__,__\___\___|\___/|_| |_| |__\___|  _   
###   | |         | |        (_)   | |                (_)     | |  
###   | |__  _   _| |__  _ __ _  __| |   ___  ___ _ __ _ _ __ | |_ 
###   | '_ \| | | | '_ \| '__| |/ _` |  / __|/ __| '__| | '_ \| __|
###   | | | | |_| | |_) | |  | | (_| |  \__ \ (__| |  | | |_) | |_ 
###   |_| |_|\__, |_.__/|_|  |_|\__,_|  |___/\___|_|  |_| .__/ \__|
###           __/ |                                     | |        
###          |___/                                      |_|        
###
###  Hybrid script for DDA/DIA, proteome/LUX/CSC analysis 
###  2024-09-25
###  Collaborative effort Martin & Marco 
####################################################################################################################################
rm(list = ls())           # Purge workspace
set.seed(123)
script_version = "_shs2.23" # version stamp on output directory
# 2.7   2024-07-29    add lux signature filter + corrected signature filtering (before not effectively due to usage of wrong df ... - not bad because no glyco anlysed recently)
# 2.8   2024-07-29    supervolcano - removal of missing degree of freedem middle plot (only required when non-limma t-test is used. limma now default >> kick out)
# 2.9   2024-07-30    Surface and Intracellular replaced with yes and no; supervolcano - CV based unique hit scatters
# 2.10  2024-08-01    GO term bar-plot for target enrichment
# 2.11  2024-08-05    STRING highlighting
# 2.12  2024-08-12    STRING highlighting v2
# 2.13  2024-08-      GO term bar-plot extended:       lectin (glyco binding); lipid binding; leucine rich repeat; IgSF
#                     Upset Plot + Subset-Heatmap for overlaps / unique IDs
# 2.14                ??
# 2.15  2024-09-25    individual comparison definintion; string data included in scatter.
# 2.16  2024-11-11    small updates; gLUX; scatter FC lines; % surface rank plot; fixed Ramos % CSPA plot
# 2.17  ?             ?
# 2.18  2024-12-04    BUG FIX - SIGNATURE COUNT - signature counting psm vs pg level (prior multiple signature-peptides/protein were summed up >> inflated protein and signature count)
# 2.19  2024-12-06    BUG FIX - Imputation annotation, median log2 intensity unique hit line, ...
# 2.20  2024-12-11    poi loop for multiple static volcanos. GO correction & refinement
# 2.21  2025-01-29    fix double labeling of datapoints. fix inf douple occurence of proteins (same proteins 2x with different CVs ...)
# 2.22  2025-02-24    unique hit imputation update
# 2.23  2025-03-20    CSPA_2.0 reintroduced and adjusted Marcos code as well as QCs (impact on CSPA reference for QC plots and in case of peptideCSC on continous updating of CSPA_2.0_*cell_type*.csv)

#Load libraries required ----------------------------------------------------------------------------------------------------------
library(dplyr)      # Data wrangling
library(tidyr)      # General functions
library(tidyverse)  # same      
library(stringr)    # Functions for handling string modifications
library(reshape2)   # Dataframe handling
library(ggplot2)      # Plotting
library(scales)       # for alpha function; otherwise alpha applies for fill and outline in suvo
library(ggrepel)      # Labelling of points in plots without overlaps
library(ggpubr)       # Multi panel plotting
library(ggthemes)     # Plotting theme
library(viridis)      # colors
library(plotly)       # interactive plots
library(htmlwidgets)  # html plots
library(gridExtra)    # plotting
library(ggplotify)     # for as.grob() transformation of qc17 qc_sample_correlation
library(cowplot)       # required for protti heatmap to be savable as composite
library(UpSetR)     # for overlap assessment of upregulated protein lists accross conditions (suvo_candi --> sigupuni)
#devtools::install_github("jpquast/protti", dependencies = TRUE, ref = "developer") # download developer version. specifically for calculate_protein_abundance peptide option
library(protti)     # https://cran.r-project.org/web/packages/protti/index.html
library(limma)      # https://bioconductor.org/packages/release/bioc/html/limma.html
library(MSstats)    # version 4.8.2        # BiocManager::install("MSstats")
library(foreach)      # parallelized for loop
library(doParallel)   # paraellized cluster setup
library(ranger)       # helper for particularization
library(BiocManager)    #
library(qvalue)       # qvalue calculation - note it is not the same as adj_pvalue. pvalue is 
library(magrittr)     # piping
library(purrr)        # helps handling lists

# ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
## ALLWAYS UPDATE BELOW  =============================================================================================================================================================================================================================
input_directory       = c("/Users/mgesell/Desktop/currentR/test/v24_LUX_TCR-CD4-CD8_FP20_HoxHoxox_6aa")
test_parameter        = "_1-2_ludo_adjp_0.7string" # leave empty unless you want to create a separate output directory for certain parameter change (e.g. "_iq_no-impute" provides output saved to "_output_iq_no-impute")
experiment_name       = "_"
data_type             = "DIA"      # DDA or DIA
replicates            = 3 # determines i) replicate ID filtering e.g. 0_vs_1 out for >1 replicate ii) t-testing. moderated t-test works only for >=2 remplicates --> in case of 1 replicate overwrite t-test protti variable with "t-test"
##
experiment_type       = "Lysate_ProtCSC_LUX_LUXCSC"     #  choose between "gLUX" "Lysate_ProtCSC_LUX_LUXCSC", "peptideCSC", and "hybrid_CSC_LUX" (last one can be used to require both deamidation and LUX ox on same peptides)
# CSPA_2.0 - only effective when:   experiment_type = "peptideCSC"
cell_type                  = "Jurkat"   # so far only "Jurkat" & "panT" implemented
dry_run                    = T # Put to false for a real run, this is important for the Jurkat CSPA, if you are only doing a test run with whatever dataset this needs to be True, so the results don't get added to the Jurkat_CSPA.csv file
CSPAv2_evidence_treshold   = 2 # In how many replicates should a PG be found at least to be considered for Jurkat CSPA
#
experiment_conditions = c("TCR", "CD4", "CD8") # order needs to match experiment_references order !!!    variable determines loop:   for (i in 1:length(experiment_conditions)) ... condition = experiment_conditions[i]
experiment_references = c("Iso", "Iso", "Iso") # order needs to match experiment_conditions order !!!    variable determines loop:   for (i in 1:length(experiment_references)) ... condition = experiment_references[i]
#


csc_filter_stringency = "consider_ONLY_nglyco_peptides"  # applied only when CSC selected; "consider_ONLY_nglyco_peptides" or "require_just_one_ngylco_peptide_in_any_sample_per_protein"
lux_filter_stringency = ""    # "lux_signature_peptides_only", "gLUX_SOG_483_signature_psm"  to .... ; usually keep "" or "non"
# stringent filter = highest data quality but missing values (loss of any non-glycosilated peptides for quantification)
qc_reference          = "meta_surfaceome" #choose one reference: cspa_2015  surfy_2018   tcsa_2021  uniprot_2023  cd_antigen   meta_surfaceome (= cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high)
poi_plotting_lists    = c("CPIs", "BindingBlockListPharmacoscopy_shilts2022_fig4d", "tcr_chains_cd3_mgmanual",
                          "Tact_markers", "Thermo_Marker_Tsubset", "tcr_signaling_GO.0050852", "LUX_targets.CPIs", "tetraspanins_galectins")  # select one colnames(poi):  
# hot lists:              "LUX_targets.CPIs"	"POI_Ben"	"FACS_candidates"	"Tcell_Act_GO.0042110"	"Bcell_Act_GO.0042113"
# Markers:                "Thermo_Marker_Tact"(1000)	"Thermo_Marker_Tsubset"	"Tact_markers" "BindingBlockListPharmacoscopy_shilts2022_fig4d"
# Surface Annotations:    "cspa_2015"	"surfy_2018"	"tcsa_2021"	"uniprot_2023"	"cd_antigen"	"VeneerProteome_High"	"CSPA2015_Surfy2018_TCSA_2021_Cdantigen_VeneerProteome_High"
#
## plotting parameters ========================================================================================================================================================================================================================================================================
# select any result subset column:    poi   meta_surfaceome      aggreg_left (sig, unique, nodf)  aggreg_right (...) left_right  (...)       cspa_2015  cd_antigen  surfy_2018   
plot_no_df_occupacy         =  0.7      # make lower confidence ids less bright
plot_label_volcano          =  "poi"               # default poi   
plot_label_inf              =  plot_label_volcano  # default plot_label_volcano
plot_fill_volcano           =  "meta_surfaceome"   # default meta_surfaceome  
plot_fill_label_inf         =  plot_fill_volcano   # default plot_fill_volcano                 
plot_fill_condition_scatter =  "meta_surfaceome" # x-y median condition signal plot (color surface or poi recommended)

## string parameters ===================================================================================================================================================================================================================================================================
# adjust colnames. they need to match what is defined in experiment_conditions; also length must be same for all (fill with "", "")
string_targets        <- data.frame(TCR = c("CD3D_HUMAN", "CD3E_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN", "CD4_HUMAN", "CD8A_HUMAN" , "CD8B_HUMAN"),  # define target & or fixed complex partners as well as condition
                                    CD4 = c("CD3D_HUMAN", "CD3E_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN", "CD4_HUMAN" , ""          , ""          ),
                                    CD8 = c("CD3D_HUMAN", "CD3E_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN", "CD8A_HUMAN", "CD8B_HUMAN", ""          ),# define target & or fixed complex partners as well as condition
                                    check.names = FALSE)  # this allows to set also numbers as colnames when they are surrounded by "" or '' 
# usually don't touch
string_linkage_categs = c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining")   # select info levels to be used these are available:   c("experimental", "database", "fusion", "neighborhood", 	"cooccurence",	"coexpression",	"textmining") 
string_min_comb_score = 700 # string combined score ranges 0-1000; 400 threshold for medium confidence, 700 for high confidence;    low throughput high accuracy experiment data results score usually +/- 600; high throughput low accuracy data scores <= 250
string_network        = "uniprot_swissprot"   # select uniprot_unreviewed (default*) or uniprot_swissprot           *data searched for upsp >> will only look at upsp protein1 anyhow. any non-upsp network info (protein2)  might be helpfull for further analysis so keep itin
if (string_network == "uniprot_unreviewed") {  # --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  string <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/string_net_unrev.csv", header = TRUE)  # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_unrev.csv
} else if (string_network == "uniprot_swissprot") {
  string <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_ppi/string/string_net_upsp.csv",  header = TRUE)  # or /Volumes/mgesell/03_DataProcessing/_resources/string/string_net_upsp.csv
}
# string interaction filter
string_filter <- string_targets %>% # reduce string to string interactor rows (otherwise huge dataframe occupying global environment)
  unlist() %>% # dataframe to vector
  unique() %>% # unique entry_names across conditions
  .[. != ""]   # omit ""
# apply string info filtering
string <- string %>%
  dplyr::select(all_of(c(string_linkage_categs, c("protein1_entry", "protein2_entry", "protein1_entry_name", "protein2_entry_name", "reviewed_protein1", "meta_surfaceome_protein1", "combined_score")))) %>%  # filter out unwanted categories  ...
  filter(!if_all(any_of(string_linkage_categs), ~ . == 0 | is.na(.))) %>%      # ... then remove all rows that lost their evidence due to columsn no longer present
  filter(combined_score >= string_min_comb_score) %>%
  filter(protein1_entry_name %in% string_filter
         |
           protein2_entry_name %in% string_filter
  )
rm(string_filter)

## data processing parameters - usually don't change ==================================================================================================================================================================================================================================
pvalue                = "adj_pvalue"          # select "pvalue" or "adj_pvalue" depending on how stringend you want to be.   recommendation for CSC is adj_pvalue    
protti_t_testing      = "moderated_t-test"# JP recommendation for protti: "moderated_t-test"; alternatives see documentation
data_processing       = "protti"          # select "protti" or "MSstats"; defines method of protein inference & testing 
unique_impute         = c("3_vs_0", "2_vs_0", "0_vs_2", "0_vs_3")             # workaround to complement protti missigness assignment:   unique hit ~ 2/3 or 4/6 replicates --> these get imputed
protti_imputation     = "ludovic"   # "no_imputation" / "ludovic" / "noise"
protti_normalization  = "median"          # chose between "median" (protti) and "none" (martins workaround see code); median = the normalised intensity is the original intensity minus the run median plus the global median. This is also the way it is implemented in the Spectronaut search engine.
protti_prot_quant     = "iq"              # choose between "iq" (maxLFQ like) or "sum";               use maxLFQ for whole proteome like data (assumes majority of features remain unchanged >> unclear if suited for enrichment data) 
protti_MAR            = 0.7  # default 0.7   e.g. 0.7x3 -->   > 2.1     # a numeric value that specifies the minimal degree of data completeness to be considered as MAR. Value has to be between 0 and 1, default is 0.7. It is multiplied with the number of replicates and then adjusted downward. The resulting number is the minimal number of observations for each condition to be considered as MAR. This number is always at least 1.
protti_MNAR           = 0.2  # default 0.2   e.g. 0.2x3 -->   < 0.6     # a numeric value that specifies the maximal degree of data completeness to be considered as MNAR. Value has to be between 0 and 1, default is 0.20. It is multiplied with the number of replicates and then adjusted downward. The resulting number is the maximal number of observations for one condition to be considered as MNAR when the other condition is complete.
protti_min_n_peptides = 1    # greater than condition.   controls developer version of protti::calculate_protein_abundance function:   filter(.data$n_peptides >= min_n_peptides)

# Specify exclcusion/inclusion criteria 
filter_exclude_files            = c("mgesell_Z_2411_118")       # e.g. c("mgesell_A_2402_186", "mgesell_A_2402_187")
filter_exclude_protein          = c("contam_", "iRT", "_YEAST", "P22629", "P00761") # use text to filter-out features by protein name;  P22629 = Streptavidin; P00761 = Trypsin
filter_include_protein          = c("")                        # specifically retain protein protein name;      
#filter_exclude_quantification  = c("TRUE")                   # !!! DISABLED !!!  (here and below)
filter_include_peptideSequence  = c()                          # use regular expression to filter-in features by peptide sequence, e.g: c("N[^P][ST]")
filter_exclude_peptideSequence  = c()                          # use to filter-out proteins by name, e.g: c("Oxidation")
filter_exclude_nonproteotypic   = c("FALSE")
## filter criteria
filter_log2fc_cutoff            = log(2, 2)   # Cutoff for Volcano plotting
filter_sig_cutoff               = 0.05        # Cutoff for p value in plotting (put p-value, -log10 is done as part of plotting)
filter_uni_hit_min_completeness = 0.5                              # >= condition; 0.5 ~ 2/3 or 2/4; set 0.6 to require for 2/3 and 3/4       fraction of replicates / condition required  (Inf assignment, filtering >> unique hit plot

# filtering IDs across replicates
if(replicates == 1) {
  filter_prot_comparison_pair_completeness    =  c("0_vs_0")  # specify comparisons that are to be dumped due to insufficient replicate coverage; format is:    replicate_ids_left "_vs_" replicate_ids_right
  protti_t_testing      = "t-test"
  } else if (replicates >= 2) {
  filter_prot_comparison_pair_completeness  =  c("0_vs_0", "1_vs_0", "0_vs_1", "1_vs_1")  # specify comparisons that are to be dumped due to insufficient replicate coverage; format is:    replicate_ids_left "_vs_" replicate_ids_right
} 

if (paste(experiment_type) == "Lysate_ProtCSC_LUX_LUXCSC") {   # for LUX use 2-0-2-0 because i) samples much dirty (crowded) ii) more peptides available
  workflow_psm_signature = "lux_signature_psm"
  workflow_signature     = "lux_signature"             # set workflow signature to highlight (csc_signature = deamidation N0.9840 within sequence motive
  ## Data comparison - Also defined in Pre-processing step
  filter_min_n_features_overall           = 1          # minimum number of features (peptideSequence+chargeState+fragmentIon combination) required per protein overall - Also defined in pre-processing
  filter_min_n_features_per_condition     = 0          # -""- required per protein per condition; usually don't set anything here - kills unique hits
  filter_min_n_features_per_sample        = 0          # -""- required per protein per  sample
  # Similar to evidence treshold but more globally
  filter_min_n_measurements_overall       = 2          #  minimum number of feature measurements (non-zero intensities) required per protein overall
  filter_min_n_measurements_condition     = 0          # -""- required per protein per condition; usually don't set anything here - kills unique hits
  filter_min_n_measurements_sample        = 0          # -""- required per protein per  sample
} else if (paste(experiment_type) == "peptideCSC"){ # default is 1-0-0_1-0-0.      2-0-0_2-0-0 causes loss of ~50% of glyco-PGs at feature filtering step. see data_prec_norm_feature_filtered
  workflow_psm_signature = "csc_signature_psm"
  workflow_signature     = "csc_signature"             # set workflow signature to highlight (csc_signature = deamidation N0.9840 within sequence motive
  ## Data comparison - Also defined in Pre-processing step
  filter_min_n_features_overall           = 1          # minimum number of features (peptideSequence+chargeState+fragmentIon combination) required per protein overall - Also defined in pre-processing
  filter_min_n_features_per_condition     = 0          # -""- required per protein per condition; usually don't set anything here - kills unique hits
  filter_min_n_features_per_sample        = 0          # -""- required per protein per  sample
  # Similar to evidence treshold but more globally
  filter_min_n_measurements_overall       = 1          # minimum number of feature measurements (non-zero intensities) required per protein overall
  filter_min_n_measurements_condition     = 0          # -""- required per protein per condition;  usually don't set anything here - kills unique hits
  filter_min_n_measurements_sample        = 0          # -""- required per protein per  sample
} else if (paste(experiment_type) == "gLUX"){ # default is 1-0-0_1-0-0.      2-0-0_2-0-0 causes loss of ~50% of glyco-PGs at feature filtering step. see data_prec_norm_feature_filtered
  workflow_psm_signature = "lux_signature_psm"
  workflow_signature    = "lux_signature"    # SOG signature (gLUX_SOG_483_signature_psm) is very rare --> use LUX signature instead
  ## Data comparison - Also defined in Pre-processing step
  filter_min_n_features_overall           = 2          # minimum number of features (peptideSequence+chargeState+fragmentIon combination) required per protein overall - Also defined in pre-processing
  filter_min_n_features_per_condition     = 0          # -""- required per protein per condition; usually don't set anything here - kills unique hits
  filter_min_n_features_per_sample        = 0          # -""- required per protein per  sample
  # Similar to evidence treshold but more globally
  filter_min_n_measurements_overall       = 2          # minimum number of feature measurements (non-zero intensities) required per protein overall
  filter_min_n_measurements_condition     = 0          # -""- required per protein per condition;  usually don't set anything here - kills unique hits
  filter_min_n_measurements_sample        = 0          # -""- required per protein per  sample
}

# CSPA2.0 - add whatever celline you like by just copying panT file (whole human upsp proteins) inserting a dummy column with only 0. call file *cell_type*_CSPA.csv. set cell_type = *your cell type* on top
glyco_standard_directory  = "/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_2.0"    # or c("/Volumes/mgesell/03_DataProcessing/3_PBMCs/primaryT8_CSPA")
cell_specific_CSPA        = read_protti(paste0(glyco_standard_directory, "/", cell_type, "_CSPA.csv")              , header = TRUE, sep = ",")  # Load "own" CSPA from standard directory
#
poi_table                 = read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv" , header = TRUE, sep = ",") 
# Load annotation files 
surface_annotations       = read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/surface.annotations.csv", header = TRUE, sep = ",")
proteome                  = read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/human_upsp_202501.csv"  , header = TRUE, sep = ",")
CSPA                      =    read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_per_cell_type.csv" , header = TRUE, sep = ",")
# reminder how large sets are
cat(paste("# Uniprot Surface: ", length(surface_annotations[!is.na(surface_annotations$uniprot_2023),]$uniprot_2023), " = " ,  round(length(surface_annotations[!is.na(surface_annotations$uniprot_2023),]$uniprot_2023)/length(proteome$entry)*100, 1), " % of total proteome",  sep = " "),
    paste("# CD Proteins:     ", length(surface_annotations[!is.na(surface_annotations$cd_antigen),]$cd_antigen)    , " = " ,  round(length(surface_annotations[!is.na(surface_annotations$cd_antigen),]$cd_antigen)/length(proteome$entry)*100, 1)    , "  % of total proteome",  sep = " "),
    paste("# CSPA Proteins:   ", length(surface_annotations[!is.na(surface_annotations$cspa_2015),]$cspa_2015)      , " = " ,  round(length(surface_annotations[!is.na(surface_annotations$cspa_2015),]$cspa_2015)/length(proteome$entry)*100, 1)      , "  % of total proteome",  sep = " "),
    paste("# Surfy Proteins:  ", length(surface_annotations[!is.na(surface_annotations$surfy_2018),]$surfy_2018)    , " = " ,  round(length(surface_annotations[!is.na(surface_annotations$surfy_2018),]$surfy_2018)/length(proteome$entry)*100, 1)    , " % of total proteome",  sep = " "),
    paste("# TCSA Proteins:   ", length(surface_annotations[!is.na(surface_annotations$tcsa_2021),]$tcsa_2021)      , " = " ,  round(length(surface_annotations[!is.na(surface_annotations$tcsa_2021),]$tcsa_2021)/length(proteome$entry)*100, 1)      , " % of total proteome",  sep = " "),
    paste("# Meta Surfaceome: ", length(surface_annotations[!is.na(surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high),]$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high),       " = ",  round(length(surface_annotations[!is.na(surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high),]$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high)/length(proteome$entry)*100, 1)      , " % of total proteome",  sep = " "),
    sep = "\n")
#
setwd(paste(input_directory))   
if (!dir.exists(    paste0(input_directory, "/_output", test_parameter, script_version))) {  dir.create(paste0(input_directory, "/_output", test_parameter, script_version))  }
output_directory  = paste0(input_directory, "/_output", test_parameter, script_version)
#
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste(input_directory))   
if (data_type == "DDA") { 
  # read, subset & trim data
  paste("expext error - need to check LUX mods to make code run")
  data_raw_ <- read_protti(paste(input_directory, "/msstats.csv", sep = "")) 
  data_raw  <- data_raw_ %>% 
    filter(!(run %in% paste(filter_exclude_files))) %>% # Exclude files based on filter_exclude_files
    rename("peptide_sequence_mod"= "peptide_sequence", "raw_prec_intensity"="intensity") %>%    ### introduce more meaningfull column names to avoid confusion later on
    mutate(stripped_peptide_sequence = gsub("\\[0.9840\\]|\\[57.0215\\]|\\[15.9949\\]|n\\[42.0106\\]|\\[31.9898\\]|\\[13.9793\\]|\\(UniMod:\\d+\\)", "", peptide_sequence_mod),  # stripped sequence
           peptide_sequence_mod = gsub("\\[483.2000\\]", "[483.2]", peptide_sequence_mod), 
           precursor = paste0(peptide_sequence_mod, "_", precursor_charge),  # precursor column
           raw_prec_intensity_log2 = log(raw_prec_intensity,2),
           lux_signature_psm           = ifelse(str_count(peptide_sequence_mod, "\\[13.9793\\]|\\[31.9898\\]") > 0, "yes", "no"),
           lux_signature_psm_n         = str_count(peptide_sequence_mod,        "\\[13.9793\\]|\\[31.9898\\]"),
           gLUX_SOG_483_signature_psm  = ifelse(str_count(peptide_sequence_mod, "\\[483.2\\]") > 0, "yes", "no")
    ) %>%  
    # trim Protein col
    separate(protein_name, into = c("trash", "entry", "entry_name"), sep = "\\|", remove = FALSE) %>%     
    mutate(entry = coalesce(entry, trash)) %>%
    mutate(entry = gsub("contam_", "", entry)) %>%
    dplyr::select(-trash, -entry_name) %>%  # entry_name removed because needs be added for DIA & if present already messes up colnames to x. y. 
    rename(PSM_id = protein_name)
  
  file_annotation <- data_raw[!duplicated(data_raw$run),c("run", "condition","bio_replicate")]    # file annotation variable
  
} else if (data_type == "DIA") { # has not bee adjusted since mid of 2023
  # LUX site 1:     H 14 = (151.03818)=  137.05891+13.98
  # LUX site 2:     H 32 =(169.04874)=  137.05891+32   
  paste("DIA loading according to FP_v22.0 DIA-NN output   (requires both msstats.csv & report.tsv)")
  file_annotation <- read_protti(paste(input_directory, "/msstats.csv", sep = "")) %>% 
    dplyr::select(run, condition, bio_replicate) %>% 
    distinct(run, .keep_all=TRUE)
  
  data_raw_  <- read_protti(paste(input_directory, "/report.tsv", sep = "")) # report.tsv contains protein_name, entry and RTs in case of DIA
  data_raw  <- data_raw_ %>% 
    dplyr::select(run, all_mapped_proteins, protein_group, modified_sequence, stripped_sequence, precursor_charge, precursor_quantity, rt, rt_start, rt_stop) %>% #  rt_start, rt_stop
    filter(!(run %in% paste(filter_exclude_files))) %>% # Exclude files based on filter_exclude_files
    rename("PSM_id" = "all_mapped_proteins", "entry" = "protein_group", "peptide_sequence_mod"= "modified_sequence","stripped_peptide_sequence"= "stripped_sequence" , "raw_prec_intensity"="precursor_quantity") %>%    ### introduce more meaningfull column names to avoid confusion later on
    mutate(peptide_sequence_mod = str_replace_all(peptide_sequence_mod, c(  # translate mods to DDA format >> same code works for all data types
      "N\\(UniMod:7\\)"   =  "N[0.9840]",
      "H\\(151.038177\\)" =  "H[13.9793]",   # 151.038177     H +13.9
      "H\\(169.04874\\)"  =  "H[31.9898]",   # 169.04874]     H +31.9
      "H\\(183.028005\\)" =  "H[45.9691]",   # 183.028005     H +13.9 +31.9   somehow FP stacks these mods in rare cases - so fare ignored in analysis
      #"W\\(\\)"  =  "H[31.9898]",   # ???    H +31.9
      "M\\(UniMod:35\\)"  =  "M[15.9949]",
      "C\\(UniMod:4\\)"   =  "C[57.0215]",
      "\\(UniMod:1\\)"    =  "n[42.0106]",
      "K\\(611.29496\\)"  =  "K[483.2]"
    ))) %>%
    mutate(precursor = paste0(peptide_sequence_mod, "_", precursor_charge),  # precursor column
           raw_prec_intensity_log2 = log(raw_prec_intensity,2),
           lux_signature_psm   = ifelse(str_count(peptide_sequence_mod, "\\[13.9793\\]|\\[31.9898\\]") > 0, "yes", "no"),
           lux_signature_psm_n = str_count(peptide_sequence_mod,        "\\[13.9793\\]|\\[31.9898\\]"),
           gLUX_SOG_483_signature_psm    = ifelse(str_count(peptide_sequence_mod, "\\[483.2\\]") > 0, "yes", "no")
    )
  
  data_raw  <- data_raw %>% 
    left_join(file_annotation, by = "run")
}  
# save.image("Workspace_Image.RData") # ====================================================================================================================================

## Retrieve and append Uniprot (sequence, protein_name = description of protein)
uniprot_ids <- unique(data_raw$entry) # protein in dataset
uniprot <-  # retrieve desiered info for identified proteins
  fetch_uniprot(
    uniprot_ids = uniprot_ids, 
    columns = c(  # see nomenclature:   https://www.uniprot.org/help/return_fields
      "protein_name",  "id", "gene_names", "gene_primary", "length", "sequence",
      "go_c", "go_f", "xref_string", "xref_pdb",
      "cc_interaction", "cc_disease", "cc_pharmaceutical", "cc_subcellular_location",
      "cc_ptm", "ft_lipid", "ft_carbohyd","ft_transmem", "ft_binding")) %>%
  rename(entry = accession,
         entry_name = id) %>%
  dplyr::select(-input_id)

data_raw_up <- data_raw %>%    # merge uniprot info and data table 
  left_join(
    y = uniprot[, c("entry", "entry_name", "sequence", "protein_name")], # at this point only add subset of uniprot info (otherwise unnecessarily inflated intermediate dataframes)
    by = "entry") %>%
  #mutate(sequence = ifelse(entry == "P04745" & is.na(sequence), stupid_protein_fix, sequence)) %>% # special line exclusively for P04745 which is not in latest uniprot version but corresponds to AMY1A_HUMAN / AMY1B_HUMAN / AMY1C_HUMAN
  mutate(stripped_peptide_sequence = case_when(                 # set cont sequences to NA so they don't affect sequence coverage calculation; correct column reintroduced below
    str_detect(PSM_id, str_c(filter_exclude_protein, collapse = "|")) ~ NA,
    TRUE ~ stripped_peptide_sequence))    %>%
  mutate(sequence                  = gsub("I|L", "IL", sequence), # WORKAROUND - seems necessary for DIA data FP_v22.0   DIA_LUXmod_Honly workflow
         stripped_peptide_sequence = gsub("I|L", "IL", stripped_peptide_sequence)) %>%  # Isoleucine and leucine are isomers with identical mass (131.17 Da). Standard mass spectrometry techniques used in proteomics, including those employed in DIA workflows, cannot distinguish between these two amino acids based on mass alone.
  find_peptide(                             # detect peptides in protein sequence
    protein_sequence = sequence,
    peptide_sequence = stripped_peptide_sequence  ) %>%
  assign_peptide_type(                      # tryptic, semitryptic, prior subesquent AA
    aa_before = aa_before,
    last_aa = last_aa,
    aa_after = aa_after  ) %>%
  calculate_sequence_coverage(              # sequence coverage
    protein_sequence = sequence,
    peptides = stripped_peptide_sequence  ) %>%
  dplyr::select(-sequence)
#xx <- data_raw_up[is.na(data_raw_up$aa_before)|is.na(data_raw_up$aa_after),]  # to test which proteins had no sequence matches (and thus get not assigned aa_before/after matches)

## protein sequence processing
# reintroduce stripped pep sequence for cont & non human proteins (i think i did this because unless removed above it results in error for cont& non-human proteins)
data_raw_up$stripped_peptide_sequence <- gsub("\\[0.9840\\]|\\[57.0215\\]|\\[15.9949\\]|n\\[42.0106\\]|\\[31.9898\\]|\\[13.9793\\]",
                                              "", data_raw$peptide_sequence_mod)  # create stripped pep seq column
data_raw_up$sequence_extended  <-   paste(     data_raw_up$aa_before,       data_raw_up$peptide_sequence_mod,       data_raw_up$aa_after,      sep = "")  # issue - contains NAs
data_raw_up$sequence_veneer    <-   paste("[", data_raw_up$aa_before, "].", data_raw_up$peptide_sequence_mod, ".[", data_raw_up$aa_after, "]", sep = "")  # issue - contains NAs
data_raw_up <- data_raw_up %>% dplyr::select(-aa_before, -last_aa, -aa_after)  # no longer needed after this step >> removed
# adjust modifications to format required by veneer; [..] = DDA format; UniMod
data_raw_up <- data_raw_up %>%
  mutate(sequence_veneer = str_replace_all(peptide_sequence_mod, c(  # translate mods to DDA format >> same code works for all data types
    "N\\[0.9840\\]|N\\(UniMod:7\\)" ="n",
    "M\\[15.9949\\]|\\(UniMod:35\\)"="M",   
    "C\\[57.0215\\]|\\(UniMod:4\\)" ="C",   
    "n\\[42.0106\\]|\\(UniMod:1\\)" ="" ,   
    "H\\[13.9793\\]"= "H",
    "H\\[31.9898\\]"= "H",
    "H\\[45.9691\\]"= "H"
  ))) 

## introduce surface annotations flavours >> QC plots
data_raw_up_surf <- data_raw_up %>%
  mutate(cd_antigen = "no", uniprot_2023 = "no",  tcsa_2021 = "no", surfy_2018 = "no", cspa_2015 = "no", cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high = "no")  %>% # create empty columns
  mutate(across(c(cd_antigen, uniprot_2023, tcsa_2021, surfy_2018, cspa_2015, cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high), ~ case_when(
    entry %in% surface_annotations[[cur_column()]] ~ "yes",
    grepl("cont", entry, ignore.case = TRUE) ~ "cont",
    TRUE ~ .x
  ))) %>%
  rename(meta_surfaceome = cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high) # CSPA_CD_orDouble_Surfy_TCSA_VeneerHigh (= meta_surfaceom_v3.0)   cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high = (meta_surfaceome_v1.0)

data_raw_up_surf$con_rep   <- paste(data_raw_up_surf$condition, data_raw_up_surf$bio_replicate, sep = "_")

## searched N-glyco deamidation WITHIN sequence motive  
# Extract Glyco Info   # Stringent filter (Deam in NxSTC)
data_raw_up_surf_gly <- data_raw_up_surf
data_raw_up_surf_gly$csc_signature_psm <- "no"
# consider placing lux mod search here
data_raw_up_surf_gly$csc_signature_psm[grepl('N\\[0.9840\\][^P][STCV]'            , data_raw_up_surf_gly$sequence_extended)] <- "yes"   # unmodified NxSTC consensus sequence (x /= Proline)
data_raw_up_surf_gly$csc_signature_psm[grepl('N\\[0.9840\\][M]\\[15.9949\\][STCV]', data_raw_up_surf_gly$sequence_extended)] <- "yes"   # x = methionin ox
data_raw_up_surf_gly$csc_signature_psm[grepl('N\\[0.9840\\][C]\\[57.0215\\][STCV]', data_raw_up_surf_gly$sequence_extended)] <- "yes"   # x = carbi stuff
data_raw_up_surf_gly <- data_raw_up_surf_gly %>% 
  mutate(csc_signature_psm = case_when(
    grepl("cont", PSM_id, ignore.case = TRUE) ~ "cont",    # overwrite data when cont
    TRUE ~ csc_signature_psm),                             # otherwise keep original value
    lux_signature_psm = case_when(
      grepl("cont", PSM_id, ignore.case = TRUE) ~ "cont",  # overwrite cont mod data
      TRUE ~ lux_signature_psm)
  )
## kick zero intensity data:    data_raw_up_surf_nZ = non Zero Intensity 
data_raw_up_surf_gly_nZ    <- data_raw_up_surf_gly[!is.na(data_raw_up_surf_gly$raw_prec_intensity),] %>% # kick NAs (DDA data has many)
  filter(raw_prec_intensity > 0)          # then remove 0 values. first kick out no intensity ions (initial or parallel unique filter will leave NA ions in >> useless)

# expand *signature* to whole protein group (so far it is psm level - to form accurate counts of yes no cont for *signature* expand to protein level)
data_raw_up_surf_gly_nZ <- data_raw_up_surf_gly_nZ %>%  # introduced 2024-12-04   shs_v2.18
  group_by(entry, con_rep, lux_signature_psm) %>%
  mutate(lux_signature = case_when(
    any(lux_signature_psm == "yes") ~ "yes",  
    grepl("cont", PSM_id, ignore.case = TRUE) ~ "cont",
    TRUE ~ "no"  ),
    csc_signature = case_when(
      any(csc_signature_psm == "yes") ~ "yes",
      grepl("cont", PSM_id, ignore.case = TRUE) ~ "cont",
      TRUE ~ "no"  ),
  ) %>%
  ungroup()

## Signature preview plots (created with function)
create_barplot_with_numbers <- function(data, data_level, signature, title) {
  tbl <- table(data[!duplicated(data[[data_level]]),][[signature]])
  bp <- barplot(tbl,
                main = paste0(signature, " ", title),
                xlab = "Values", ylab = "Frequency",
                col = c("red", "orange", "darkgreen"))
  text(bp, tbl, labels = tbl, pos = 3, cex = 1.6, offset = 0.5)
}
# Plot for peptide_sequence_mod
create_barplot_with_numbers(data= data_raw_up_surf_gly_nZ, data_level= "peptide_sequence_mod"     , signature= workflow_psm_signature, title= "Feature Count (peptide_sequence_mod)")
# Plot for entry
create_barplot_with_numbers(data= data_raw_up_surf_gly_nZ, data_level= "entry"                    , signature= workflow_signature     , title= "Protein Count (entry)")

## CRITICAL data filtering based on specific signature (mandatory for CSC - optional for LUX)
if (csc_filter_stringency == "consider_ONLY_nglyco_peptides" & lux_filter_stringency == "lux_signature_peptides_only" & experiment_type == "hybrid_CSC_LUX") {         # more missing values but higher confidence in datapoints
  # double signature filtering (for CSC-LUX or LUX-CSC hybrid workflow)
  data_raw_up_surf_nZ_glyF <- filter(data_raw_up_surf_gly_nZ, csc_signature == "yes" & lux_signature == "yes")             # only those peptides with direct glyco evidence are considerd for further analsis
  paste("ATTENTION - requiring CSC signature >>> !!! all other features exlcuded from quantification !!!!")
}  else if (csc_filter_stringency == "consider_ONLY_nglyco_peptides" & experiment_type == "peptideCSC") {     # more missing values but higher confidence in datapoints
  # CSC stringent filtering
  data_raw_up_surf_nZ_glyF <- filter(data_raw_up_surf_gly_nZ, csc_signature == "yes")    # only those peptides with direct glyco evidence are considerd for further analsis
  paste("ATTENTION - requiring CSC signature >>> !!! all other features exlcuded from quantification !!!!")
} else if (csc_filter_stringency == "require_just_one_ngylco_peptide_in_any_sample_per_protein" & experiment_type == "peptideCSC") {  # non-glyco peptides might distort quantification
  # CSC loose filtering
  data_raw_up_surf_nZ_glyF <- data_raw_up_surf_gly_nZ[data_raw_up_surf_gly_nZ$entry %in% unique(data_raw_up_surf_gly_nZ[data_raw_up_surf_gly_nZ$csc_signature == "yes",]$entry) ,]
  paste("ATTENTION - requiring only one glyco peptide to be identified in any sample to >>> !!! NON-gylco peptides might contribute to quantification !!!!")
} else if (lux_filter_stringency == "lux_signature_peptides_only" & experiment_type == "Lysate_ProtCSC_LUX_LUXCSC" | lux_filter_stringency == "lux_signature_peptides_only" & experiment_type == "gLUX") {  # non-glyco peptides might distort quantification
  # LUX stringent filtering
  data_raw_up_surf_nZ_glyF <- filter(data_raw_up_surf_gly_nZ, lux_signature == "yes")    # keep only lux peptides
  paste("ATTENTION - requiring LUX signature >>> !!! all other features exlcuded from quantification !!!!")
} else if (lux_filter_stringency == "gLUX_SOG_483_signature_psm") {  # non-glyco peptides might distort quantification
  # LUX stringent filtering
  data_raw_up_surf_nZ_glyF <- filter(data_raw_up_surf_gly_nZ, gLUX_SOG_483_signature_psm == "yes")    # keep only lux peptides
  paste("ATTENTION - requiring LUX signature >>> !!! all other features exlcuded from quantification !!!!")
} else { #  no signature filtering applied
  # no filtering
  data_raw_up_surf_nZ_glyF <- data_raw_up_surf_gly_nZ  
  paste("ATTENTION - signature independent analysis")
}

# clean workspace & msstats.csv file from FP is modSequence based >> filter for file specific PSMs (table has NAs that need to be removed)
rm(data_raw_up_surf,data_raw_up)
write.csv(data_raw_up_surf_gly_nZ, file = paste0(output_directory, "/data_raw_glyco_nZ_gyIFcsc.csv"))
write.csv(data_raw_up_surf_gly_nZ %>% filter(lux_signature == "yes") %>% dplyr::select(condition, bio_replicate, con_rep, entry_name, protein_name, peptide_sequence_mod, sequence_veneer, cd_antigen, uniprot_2023, tcsa_2021, surfy_2018, cspa_2015, meta_surfaceome)
          , file = paste0(output_directory, "/LUX_peptides.csv"))

## SURFACEOME QC PLOTTING ______________________________________________________________________________________________________________________________
# QC:   how many candidate proteins are already annotated (based on surface_annotation$... )
glycoAnnotDf_dummy <- data_raw_up_surf_nZ_glyF[!duplicated(data_raw_up_surf_nZ_glyF$entry),]
glycoAnnotDf <- data.frame(Category = c("CD", "UP_2023", "TCSA","SURFY" ,"CSPA" , "meta_surfaceome"),
                           Total    = replicate(6, sum(str_count(glycoAnnotDf_dummy$cd_antigen, "yes|no"))),
                           yes      = c(sum(str_count(glycoAnnotDf_dummy$cd_antigen,   "yes")),
                                        sum(str_count(glycoAnnotDf_dummy$uniprot_2023, "yes")),
                                        sum(str_count(glycoAnnotDf_dummy$tcsa_2021,    "yes")),
                                        sum(str_count(glycoAnnotDf_dummy$surfy_2018,   "yes")),
                                        sum(str_count(glycoAnnotDf_dummy$cspa_2015,    "yes")),
                                        sum(str_count(glycoAnnotDf_dummy$meta_surfaceome, "yes")))
)
glycoAnnotDf   <- glycoAnnotDf %>% mutate(no = Total-yes)
glycoAnnotDf   <- glycoAnnotDf %>% dplyr::select(Category, no, yes)
glycoAnnotDf_L <- pivot_longer(glycoAnnotDf, cols = -Category, names_to = "Annotated", values_to = "Count")
cat(paste("  >> Data overlap with surface annotations.", 
          "  >> NOTE: if peptideCSC workflow was selected, the data is signature filtered (>> don't infere enrichment purity!)", sep = "\n"))
glycoAnnotDf 
rm(glycoAnnotDf_dummy)

# set the signature info to be used in plotting
if (experiment_type == "peptideCSC") {
  data_raw_up_surf_gly_nZ$plot_signature <-  data_raw_up_surf_gly_nZ$csc_signature
} else if (experiment_type == "Lysate_ProtCSC_LUX_LUXCSC") {
  data_raw_up_surf_gly_nZ$plot_signature <-  as.character(data_raw_up_surf_gly_nZ$lux_signature)
} else if (experiment_type == "gLUX") {
  data_raw_up_surf_gly_nZ$plot_signature <-  as.character(data_raw_up_surf_gly_nZ$gLUX_SOG_483_signature_psm)
}
# CSPA PGs (1% FDR) - protein annotated in *qc_reference*
qc_plot_surface_1 <-  data_raw_up_surf_gly_nZ %>% 
  dplyr::select(con_rep, entry, !!qc_reference) %>%
  unique() %>%
  ggplot(aes_string("con_rep", fill = qc_reference))+
  geom_bar(color = "black", position = "dodge") +
  theme_classic()+
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_text(stat='count', aes(label=..count..), size = 5 ,position = "dodge", angle = 45)+
  scale_fill_brewer(palette = "BuPu")+
  labs(y = "# Identifications", x=" ", title = "PGs (1% FDR)")
# Glyco PGs (1% FDR) - proteins infered form peptides carrying deamidation mass signature within NxSTCV (=formerly N-glycosylated peptide/protein)
qc_plot_surface_2 <-  data_raw_up_surf_gly_nZ %>% 
  dplyr::select(con_rep, entry, plot_signature) %>%
  unique() %>%
  ggplot(aes(con_rep, fill = plot_signature))+
  geom_bar(color = "black", position = "dodge") +
  theme_classic()+
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_text(stat='count', aes(label=..count..), size = 5 ,position = "dodge", angle = 45)+
  scale_fill_brewer(palette = "BuPu")+
  labs(y = "# Identifications", x=" ", title = "*workflow signature* PGs (1% FDR)")
# Relative Precursor Quantities     # chose one:   cspa_2015  surfy_2018   tcsa_2021  uniprot_2023  cd_antigen
qc_plot_surface_3 <-  data_raw_up_surf_gly_nZ[data_raw_up_surf_gly_nZ$plot_signature != "cont",] %>%
  group_by(con_rep, !!sym(qc_reference)) %>%
  dplyr::mutate(TotInt = sum(raw_prec_intensity)) %>%
  dplyr::select(con_rep, !!sym(qc_reference), TotInt) %>%
  unique() %>%
  group_by(con_rep) %>%
  dplyr::mutate(SampleInt = sum(TotInt)) %>%
  ungroup() %>%
  dplyr::mutate(RelInt = TotInt/SampleInt) %>%
  ggplot(aes_string("con_rep", "RelInt", fill = qc_reference))+
  scale_fill_brewer(palette = "BuPu")+
  theme_classic()+
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_col(scale = fill, color = "black")+
  labs(y = "Relative raw_prec_intensity", x=" ", title = "Relative Raw Precursor Signal (1% FDR)", y = "PG Quantities")
## % surface signal rank plot
qc_plot_surface_4 <-  data_raw_up_surf_gly_nZ %>%
  filter(plot_signature != "cont") %>%
  group_by(con_rep, !!sym(qc_reference)) %>%
  dplyr::mutate(TotInt = sum(raw_prec_intensity)) %>%
  dplyr::select(con_rep, !!sym(qc_reference), TotInt) %>%
  unique() %>%
  group_by(con_rep) %>%
  dplyr::mutate(SampleInt = sum(TotInt)) %>%
  ungroup() %>%
  dplyr::mutate(RelInt = TotInt / SampleInt) %>%
  filter(meta_surfaceome == "yes") %>%  # Filter for meta_surfaceome = "yes" and arrange by RelInt
  arrange(desc(RelInt)) %>%
  mutate(con_rep = factor(con_rep, levels = unique(con_rep))) %>% # Create a factor for con_rep to maintain order in the plot
  ggplot(aes_string("con_rep", "RelInt", fill = qc_reference)) +
  scale_fill_brewer(palette = "BuPu") +
  theme_classic() +
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_col(color = "black") +
  labs(y = "Relative raw_prec_intensity", 
       x = " ", 
       title = "*surface_annotation* related Precursor Signal (1% FDR)")
# Relative Glyco Peptide Signal (1% FDR)
qc_plot_surface_5 <-  data_raw_up_surf_gly_nZ[data_raw_up_surf_gly_nZ$plot_signature != "cont",] %>%
  group_by(con_rep, plot_signature) %>%
  dplyr::mutate(TotInt = sum(raw_prec_intensity)) %>%
  dplyr::select(con_rep, plot_signature, TotInt) %>%
  unique() %>%
  group_by(con_rep) %>%
  dplyr::mutate(SampleInt = sum(TotInt)) %>%
  ungroup() %>%
  dplyr::mutate(RelInt = TotInt/SampleInt) %>%
  ggplot(aes(con_rep, RelInt, fill = plot_signature))+
  scale_fill_brewer(palette = "BuPu")+
  theme_classic()+
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_col(scale = fill, color = "black")+
  labs(y = "Relative raw_prec_intensity", x=" ", title = "Relative Raw *workflow signature* Peptide Signal (1% FDR)", y = "PG Quantities")
# Peptide IDs *Annotation* & CSC-Signature (Glyco)
qc_plot_surface_6 <-  data_raw_up_surf_gly_nZ %>%   
  dplyr::select(con_rep, peptide_sequence_mod, !!sym(qc_reference)) %>%
  unique() %>%
  ggplot(aes_string("con_rep", fill = qc_reference))+
  geom_bar(color = "black", position = "dodge") +
  theme_classic()+
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_text(stat='count', aes(label=..count..), size = 5 ,position = "dodge", angle = 45)+
  scale_fill_brewer(palette = "BuPu")+
  labs(y = "# Identifications", x=" ", title = "Peptide Count (1% FDR)")+
  scale_y_continuous(trans = "log10")
#
qc_plot_surface_7 <-  data_raw_up_surf_gly_nZ %>% 
  dplyr::select(con_rep, peptide_sequence_mod, plot_signature) %>%
  unique() %>%
  ggplot(aes(con_rep, fill = plot_signature))+
  geom_bar(color = "black", position = "dodge") +
  theme_classic()+
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_text(stat='count', aes(label=..count..), size = 5 ,position = "dodge", angle = 45)+
  scale_fill_brewer(palette = "BuPu")+
  labs(y = "# Identifications", x=" ", title = "*workflow signature* Peptide Count (1% FDR)")+
  scale_y_continuous(trans = "log10")
# proportion of hits in Annotation Flavors 
qc_plot_surface_8 <-  ggplot(glycoAnnotDf_L, aes(x = Category, y = Count, fill = Annotated))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x= "Annotation", y = "Overlap") +
  scale_fill_manual(values = c("no" = "gray", "yes" = "black"))+
  labs(title = "Annotation Overlap - Candidate PGs")+ 
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45))
# Target vs. off-target Signal Plot       CSC-Signature (Glyco) peptide intensity vs. non-Glyco peptide signal
qc_plot_surface_9 <- data_raw_up_surf_gly_nZ[data_raw_up_surf_gly_nZ$plot_signature != "cont",] %>% 
  dplyr::select(con_rep, plot_signature, raw_prec_intensity) %>% mutate(plot_signature = as.character(plot_signature)) %>%
  group_by(term = plot_signature) %>%
  ggplot(aes(x = con_rep, y = raw_prec_intensity, fill = plot_signature)) + 
  geom_boxplot()+
  theme_classic()+
  scale_y_continuous(trans = "log10")+
  labs(title = "Enrichment Plot: Target vs Off-Target Raw Intensity (*workflow signature* based)", x="")+
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45))
## Match CSPA with runs
# Extract CSPA Proteins for the given cell type
expected_proteins <- cbind(cell_specific_CSPA[,c(1,5)])
colnames(expected_proteins)[2] <- cell_type
expected_proteins <- expected_proteins[expected_proteins[[cell_type]]==1,]
# Initialize fraction vector
data_reps <- data_raw_up_surf_gly_nZ
reps <- unique(data_reps$con_rep)
proportions_cellType <- rep(0, times=length(reps))
counter=1
# Iterate through runs and count how many of the expected proteins are in each run
for (v in reps) {
  proportions_cellType[counter] <- table(unique(data_reps[data_reps$con_rep==v,]$entry) %in% expected_proteins$entry)[2]
  counter=counter+1  
}
names(proportions_cellType) <- reps
# Set the possible max number of matches
limit <- length(expected_proteins[[cell_type]])
# Calculate the fraction of counts relative to the limit
fraction <- pmin(proportions_cellType, limit) / limit # I totally forgot why pmin is here but it works so just leave it
remainder <- 1 - fraction
names(fraction) <- c()
names(remainder) <- c()
# Create a data frame with the condition names, fraction, and remainder
data_CSPA_fraction <- rbind(data.frame(con_rep = reps, Fraction = remainder, Match=rep("no Match", times=length(remainder))), data.frame(con_rep = reps, Fraction = fraction, Match=rep("Match", times=length(fraction))))
data_CSPA_fraction$Match <- as.factor(data_CSPA_fraction$Match)
data_CSPA_fraction$Match <- relevel(data_CSPA_fraction$Match, ref = "no Match")
# Stacked barplot 
qc_plot_surface_10 <-  ggplot(data_CSPA_fraction, aes(fill=Match, y=Fraction, x=con_rep)) + 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  labs(x= "", y = paste("% of PGs matched with\n CSPA of:", cell_type, "cells", sep = " ")) +
  theme(text = element_text(size =20), axis.text.x = element_text(hjust = 1, angle = 45))+
  scale_fill_manual(values = c("no Match" = "gray", "Match" = "black"))+
  labs(title = "") +
  geom_hline(yintercept = 0.22, color = "red", linetype = "solid", size = 1) 

# Saving QC Plots
####################################################################################################################################
# pdf
pdf(paste0(output_directory, "/_qc_surface.pdf"), onefile = TRUE, width = 3+length(unique(data_raw_up_surf_gly_nZ$con_rep))*0.75, height = 4)
for (i in 1:10) {
  current_plot <- get(paste0("qc_plot_surface_", i))
  print(current_plot)
}
dev.off()
# composite png
png(paste0(output_directory,"/_qc_surface_composite.png"), width = 30, height = 16, units = "in", res = 150)
print(plot_grid(plotlist = list(qc_plot_surface_1, qc_plot_surface_2, qc_plot_surface_3, qc_plot_surface_4, qc_plot_surface_5, qc_plot_surface_6, qc_plot_surface_7, qc_plot_surface_8, qc_plot_surface_9, qc_plot_surface_10), ncol = 4))
dev.off()
## END OF SURFACE QC PLOTS __________________________________________________________________________________________________________________________________

## DATA FILTERING - according to specifications on -------------------------------------------------------------------------------------------------------------
data_raw_pre_filtered <- data_raw_up_surf_nZ_glyF
includeFeature = rep(TRUE, nrow(data_raw_pre_filtered))
for (filter in filter_exclude_protein) {
  includeFeature = includeFeature & !grepl(filter, data_raw_pre_filtered$PSM_id,                    ignore.case=T)     } # exclude
for (filter in filter_exclude_peptideSequence) {
  includeFeature = includeFeature & !grepl(filter, data_raw_pre_filtered$stripped_peptide_sequence, ignore.case=T)     } # exclude 
for (filter in filter_include_protein) {
  includeFeature = includeFeature & grepl(filter,  data_raw_pre_filtered$PSM_id,                    ignore.case=T)  }  # include
for (filter in filter_include_peptideSequence) {
  includeFeature = includeFeature & grepl(filter,  data_raw_pre_filtered$stripped_peptide_sequence, ignore.case=T)  }  # include
# for (filter in filter_exclude_quantification) {
#   includeFeature = includeFeature & !grepl(filter, data_raw_pre_filtered$,                        ignore.case=T)  } # DISABLED 
for (filter in filter_exclude_nonproteotypic) {
  includeFeature = includeFeature & !grepl(filter, data_raw_pre_filtered$peptide_sequence_mod,      ignore.case=T)  }  # exclude
cat("Will use",sum(includeFeature),"out of",nrow(data_raw_pre_filtered),"features for Protti / MSstats analysis.")
data_raw_pre_filtered <- data_raw_pre_filtered[includeFeature, ] %>%
  dplyr::select(-PSM_id)
## END data pre-filtering __________________________________________________________________________________________________________________________________________

## clean up dataframe
data_raw_pre_filtered <- data_raw_pre_filtered %>% dplyr::select( -sequence_extended, -sequence_veneer)

## CSPA_2.0 thanks Marco & Martin __________________________________________________________________________________________________________________________________________
if (experiment_type == "peptideCSC") { 
  # Extracting per condition the Glyco PGs and merging with Jurkat CSPA (needs to be refined, maybe only do after MSstats to get more confident identifications)
  # Separate the conditions
  conditions <- unique(data_raw_pre_filtered$condition)
  # Extract all ever discovered PGs per condition and save in Output directory, the output contains which PGs were identified and in how many of the replicates, so that e.g. proteins only found in one replicate might be looked at with more doubt
  # Additionally the CSPA is filled up
  for (cond in conditions) {
    # Naming
    condition_name <- paste("condition", cond, sep = "_")
    # Subsetting CSC signature PGs
    subset_ib <- subset(data_raw_pre_filtered, subset = condition==cond & csc_signature=="yes", select = c(entry, con_rep))
    subset_ib <- table(subset_ib)
    # We don't care about how many Glyco peptides per PG so just set to one if there were any
    subset_ib[subset_ib > 0] <- 1
    # Add up the evidence
    row_sums <- rowSums(subset_ib)
    # Create named vector
    assign(condition_name, row_sums)
    
    # Search cell_specific_CSPA frame for column with matching sample time and cell state, if it does not exist yet, create it
    if (!(cond %in% names(cell_specific_CSPA))) {
      # If the column doesn't exist, create it as an empty column 
      cell_specific_CSPA[[cond]] <- 0
      print(paste(cond, "was not yet present, the corresponding column was created now", sep = " "))
    } 
    matching_cell_state <- which(names(cell_specific_CSPA) == cond)
    # With the index known the goal is to iterate over each Glyco PG in the condtition, checkîng if it is already in the list
    # if it is not, initialize it, if yes add the replicate evidence to already existing replicate evidence (depending on threshhold) 
    evidence_filter <- get(condition_name)>=CSPAv2_evidence_treshold
    addition_counter <- 0
    print(paste(length(get(condition_name)[evidence_filter]), "Glyco PGs passed evidence treshhold out of", length(get(condition_name)), "in condition", cond, sep = " "))
    for (prot in names(get(condition_name)[evidence_filter])) {
      
      if (!(prot %in% cell_specific_CSPA$entry)) {
        # If the row doesn't exist, create it as an empty row 
        cell_specific_CSPA <- rbind(cell_specific_CSPA, c(prot, rep(0, times=ncol(cell_specific_CSPA)-1)))
        cell_specific_CSPA[, 2:ncol(cell_specific_CSPA)] <- lapply(cell_specific_CSPA[, 2:ncol(cell_specific_CSPA)], as.numeric) # has to be there for the moment
        # print(paste(prot, "was not yet present, the corresponding entry was created now with condition", cond, sep = " "))
        addition_counter <- addition_counter + 1
      }
      matching_row <- which(cell_specific_CSPA$entry == prot)
      cell_specific_CSPA[matching_row, matching_cell_state] <- cell_specific_CSPA[matching_row, matching_cell_state] + 1
    }
    cell_specific_CSPA[[cond]][1] <- cell_specific_CSPA[[cond]][1] + 1
    print(paste(addition_counter, "new proteins were added to CSPA with", cond, sep = " "))
    
    # column 5 must be CSPA_2.0_*cell_type* (fyi col 6 = CSPA_1.0_*cell_type* (make all 0 if not in CSPA 1.0), col 7 = CSPA_1.0 (full)
    cell_specific_CSPA <- cell_specific_CSPA %>%
      mutate(across(5, ~ {
        sum_vals <- rowSums(cbind(cur_data()[6], cur_data()[8:ncol(cur_data())]), na.rm = TRUE)  # exclude column 7 which is entire CSPA_1.0
        ifelse(sum_vals > 0, 1, sum_vals)
      }))
    
  }
  # Finally save the Jurkat CSPA again in the directory it is situated in
  if (dry_run != T) {
    print("Updated Jurkat CSPA will be saved...")
    setwd(paste(glyco_standard_directory))
    write.csv(cell_specific_CSPA, file=paste0(cell_type, "_CSPA.csv"), row.names = FALSE) 
    setwd(paste(output_directory))  
  } else {
    print("Jurkat CSPA was not updated, change Dry_run parameter if this was no test run")
  }
}
# end CSPA_2.0 =========================================================================================================================================================
#
## Processing Pipeline Protti or MSstats ==============================================================================================================================
if (data_processing == "protti"){    #-------- Protti -------- Protti -------- Protti -------- Protti -------- Protti -------- Protti -------- Protti -------- Protti 
  # Normalize data if specified
  if (protti_normalization == "median") { 
    data_prec_norm <- data_raw_pre_filtered %>% 
      normalise(sample = run,
                intensity_log2 = raw_prec_intensity_log2,
                method = protti_normalization) %>%
      rename(normalised_prec_intensity_log2 = normalised_intensity_log2) %>% # to avoid confusion with protein quantitiy that is created later on
      mutate(normalised_prec_intensity = 2^normalised_prec_intensity_log2) # for QC plot CV calculation
  } else if (protti_normalization == "none") { 
    data_prec_norm <- data_raw_pre_filtered %>% 
      mutate(normalised_prec_intensity_log2 = raw_prec_intensity_log2) %>% # use unnormalized as normalized (keep nomenclature ...norm... for code compatibility)
      mutate(normalised_prec_intensity = 2^normalised_prec_intensity_log2) # for QC plot CV calculation
  }
  
  data_prec_norm_feature_filtered <- data_prec_norm %>% 
    filter(!is.na(normalised_prec_intensity_log2)) %>%  # comment out if desired to keep dataframe matrix complete                   
    # caluclate how many ...
    group_by(entry)            %>%                                    
    mutate(n_measurments_overall =     length(precursor[!is.na(normalised_prec_intensity_log2)]),
           n_features_overall    =     n_distinct(precursor[!is.na(normalised_prec_intensity_log2)])) %>%        # ... feature measurements per protein overall 
    ungroup() %>%      
    group_by(entry, condition) %>%                                    
    mutate(n_measurments_per_condition = n_distinct(precursor[!is.na(normalised_prec_intensity_log2)]),
           n_features_per_condition    = n_distinct(precursor[!is.na(normalised_prec_intensity_log2)])) %>%  # ... feature measurements per protein in condition 
    ungroup() %>%      
    group_by(entry, con_rep)   %>%                                    
    mutate(n_measurments_per_sample = n_distinct(precursor[!is.na(normalised_prec_intensity_log2)]),
           n_features_per_sample    = n_distinct(precursor[!is.na(normalised_prec_intensity_log2)])) %>%     # ... feature measurements per sample (con_rep)
    ungroup() %>%          
    # precursor level data FILTERING !!! feature & measurement 
    filter(
      n_measurments_overall       >= filter_min_n_measurements_overall   &
        n_measurments_per_condition >= filter_min_n_measurements_condition &
        n_measurments_per_sample    >= filter_min_n_measurements_sample    &
        n_features_overall          >= filter_min_n_features_overall       &
        n_features_per_condition    >= filter_min_n_features_per_condition &
        n_features_per_sample       >= filter_min_n_features_per_sample
    )
  
  # precursor(feature) ranking ==================================================================================================================
  data_prec_norm_feature_filtered <- data_prec_norm_feature_filtered %>% # calculate medians
    group_by(precursor) %>% 
    mutate(global_median_log2_feature_intensity    = median(normalised_prec_intensity_log2, na.rm = TRUE)) %>%
    group_by(precursor, condition) %>% 
    mutate(condition_median_log2_feature_intensity = median(normalised_prec_intensity_log2, na.rm = TRUE)) %>%
    ungroup() 
  data_prec_norm_feature_filtered <- data_prec_norm_feature_filtered %>%   # global protein rank   
    arrange(desc(global_median_log2_feature_intensity)) %>% 
    distinct(precursor, global_median_log2_feature_intensity) %>%
    mutate(rank_feature_global = row_number()) %>%
    right_join(data_prec_norm_feature_filtered, by = c("precursor", "global_median_log2_feature_intensity")) %>%
    arrange(row_number()) 
  data_prec_norm_feature_filtered <- data_prec_norm_feature_filtered %>%   # condition protein rank
    arrange(desc(condition_median_log2_feature_intensity)) %>% 
    group_by(condition) %>% #
    distinct(precursor, condition, condition_median_log2_feature_intensity) %>%
    mutate(rank_feature_condition = row_number()) %>%
    ungroup() %>%
    right_join(data_prec_norm_feature_filtered, by = c("precursor", "condition", "condition_median_log2_feature_intensity")) %>%
    arrange(row_number())
  
  # precursor cv calculation =====================================================================================================================
  data_prec_norm_feature_filtered <- data_prec_norm_feature_filtered %>%
    group_by(precursor, condition) %>%
    mutate(
      mean_precursor = mean(normalised_prec_intensity, na.rm = TRUE),
      sd_precursor = sd(normalised_prec_intensity, na.rm = TRUE),
      cv_precursor = (sd_precursor / mean_precursor) * 100 ) %>%
    dplyr::select(-sd_precursor, -mean_precursor) %>%
    ungroup()  %>%
    mutate(meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no")) %>% # for visualization
    arrange(meta_surfaceome)  # ensure meta_surfaceome"yes" on bottom of df and thus top of figures
  
  
  ## protein abundance calculation - protti developer version function (standard package function misses "min_n_peptides") ==============================================
  data_prot <-  calculate_protein_abundance(  
    data        = data_prec_norm_feature_filtered,
    min_n_peptides = protti_min_n_peptides,  # greater than condition that controls min peptide / protein filter in calculate_protein_abundance
    sample      = run,
    protein_id  = entry,
    precursor   = precursor,
    peptide     = peptide_sequence_mod,
    intensity_log2 = normalised_prec_intensity_log2,
    method      = protti_prot_quant,    
    for_plot    = FALSE, # if set TRUE the precursor and protein quantitiy will be below each other - protein intensity marked by "protein_intensity" in precursor column (filter desired with: %>% filter(precursor == "protein_intensity"))
    retain_columns = c(condition, bio_replicate, con_rep, entry, entry_name, protein_name,  
                       n_measurments_overall, n_measurments_per_condition, n_measurments_per_sample, 
                       n_features_overall, n_features_per_condition, n_features_per_sample)) %>%      # reduces dataset to protein intensity data - critical! required to account for way protti is written 
    rename(normalised_protein_intensity_log2 = normalised_prec_intensity_log2)  # needs to be removed if "filter(precursor == "protein_intensity")" line is deleted
  
  data_prot_mis <- data.frame(condition = character(), bio_replicate = character(), con_rep = character(), entry = character(), entry_name = character(),
                              protein_name = character(), n_measurments_overall = numeric(), n_measurments_per_condition = numeric(), n_measurments_per_sample = numeric(), n_features_overall = numeric(),
                              n_features_per_condition = numeric(), n_features_per_sample = numeric(), run = character(), normalised_protein_intensity_log2 = numeric(), comparison = character(),
                              missingness = character(), replicate_ids_right = numeric(), features_left = numeric(), features_right = numeric(), replicate_ids_left = numeric(),
                              features_left_right_overall = character(), n_obs_per_comparison = numeric(), replicate_ids_left_right = character()      )
  for (i in 1:length(experiment_references)) {
    paste("dtermine missingness for reference: ", experiment_references[i], "     condition: ", experiment_conditions[i])
    #paste("if ERROR message: Column`n_detect_control doesn't exist - you tried to form comparison conditionX vs conditionX (itself) or wrong naming --> check your definition of experiment_conditions & experiment_references")
    #protein missingness calculation; expands df to contain all combinations of condition-entry-comparison (NOTE: Iso data is duplicated once for each comparison)
    data_prot_missing <- data_prot %>% 
      filter(condition == experiment_references[i] | condition == experiment_conditions[i]) %>% # subset df for comparison pair
      assign_missingness(sample = run,
                         condition = condition,
                         grouping  = entry, # or precursor
                         intensity = normalised_protein_intensity_log2,
                         ref_condition = experiment_references[i],  #
                         completeness_MAR  = protti_MAR,   # default 0.7   e.g. 0.7x3 -->   > 2.1     # a numeric value that specifies the minimal degree of data completeness to be considered as MAR . Value has to be between 0 and 1, default is 0.7. It is multiplied with the number of replicates and then adjusted downward. The resulting number is the minimal number of observations for each condition to be considered as MAR. This number is always at least 1.
                         completeness_MNAR = protti_MNAR,  # default 0.2   e.g. 0.2x3 -->   < 0.6     # a numeric value that specifies the maximal degree of data completeness to be considered as MNAR. Value has to be between 0 and 1, default is 0.20. It is multiplied with the number of replicates and then adjusted downward. The resulting number is the maximal number of observations for one condition to be considered as MNAR when the other condition is complete.
                         retain_columns = c(condition, bio_replicate, con_rep, entry, entry_name, protein_name,  
                                            n_measurments_overall, n_measurments_per_condition, n_measurments_per_sample, 
                                            n_features_overall, n_features_per_condition, n_features_per_sample)
      ) 
    
    # group by relevant columns and fill in missing values based on non-NA entries within each group
    data_prot_missing <- data_prot_missing %>%
      group_by(run) %>%    # complete run specific info 
      mutate(across(c(condition, bio_replicate, con_rep), ~ifelse(is.na(.), first(.[!is.na(.)]), .))) %>%
      ungroup()  %>%
      group_by(entry) %>%  # complete entry specific info 
      mutate(across(c(entry_name, protein_name, n_measurments_overall, n_features_overall), ~ifelse(is.na(.), first(.[!is.na(.)]), .))) %>%
      ungroup()  %>% 
      group_by(entry, condition) %>%    # complete entry/condition specific info 
      mutate(across(c(n_measurments_per_condition, n_features_per_condition), ~ifelse(is.na(.), first(.[!is.na(.)]), .))) %>%
      ungroup()  %>% 
      group_by(entry, con_rep) %>%    # complete entry/condition/bio_replicate  specific info 
      mutate(across(c(n_measurments_per_sample, n_features_per_sample), ~ifelse(is.na(.), first(.[!is.na(.)]), .))) %>%
      ungroup()
    
    # calculate features left/right per comparison & replicate_ids left/right; ***for reference condition (duplicated rows for each comparison) required to re-adjust right sided features and replicate_ids
    data_prot_missing <- data_prot_missing %>% 
      group_by(entry, condition, comparison) %>%
      mutate(replicate_ids_right = sum(!is.na(n_features_per_sample) & is.numeric(n_features_per_sample))) %>%
      group_by(entry, comparison) %>%
      mutate(
        features_left      = paste(ifelse(any(condition == experiment_references[i] & !is.na(n_features_per_condition)), n_features_per_condition[condition == experiment_references[i] & !is.na(n_features_per_condition)], 0)),
        features_right     = ifelse(is.na(n_features_per_condition), 0, n_features_per_condition),
        replicate_ids_left = sum(!is.na(n_features_per_sample)[condition == experiment_references[i]] & is.numeric(n_features_per_sample[condition == experiment_references[i]])),
        # *** if reference condition row - retrieve replicate_id_right info from corresponding partner & save to dummy column to not overwrite correct values 
        features_right_dummy      = ifelse(condition == experiment_references[i],        features_right[condition == str_remove(comparison, paste0("_vs_",experiment_references[i]))],  NA),     # ***
        replicate_ids_right_dummy = ifelse(condition == experiment_references[i],   replicate_ids_right[condition == str_remove(comparison, paste0("_vs_",experiment_references[i]))],  NA)) %>% # ***for reference condition, write corrected value to dummy column
      ungroup() %>%
      mutate(features_right      = ifelse(is.na(features_right_dummy)     , features_right     , features_right_dummy     ),      # ***
             replicate_ids_right = ifelse(is.na(replicate_ids_right_dummy), replicate_ids_right, replicate_ids_right_dummy)) %>%  # ***if dummy column is NA (means this row is not ref_conditon) keep original value reference condition. Otherwise if it has a value (meanks its a rev_condition row) overwrite right_ids with correct value
      dplyr::select(-features_right_dummy, -replicate_ids_right_dummy) %>% # ***remove dummy column
      mutate(
        features_left_right_overall = paste(features_left, "_vs_", features_right, "   (", n_features_overall, " total)", sep = ""),
        n_obs_per_comparison = replicate_ids_left + replicate_ids_right,
        replicate_ids_left_right = paste(replicate_ids_left, "_vs_", replicate_ids_right, sep = "")
      )  
    # append each comparison to meta dataframe (data_prot_diff_abundance)
    data_prot_mis <- rbind(data_prot_mis, data_prot_missing)
  }
  # customize missingness for unique hits to 
  data_prot_mis <- data_prot_mis %>% 
    mutate(missingness = case_when(replicate_ids_left_right %in% unique_impute ~ "MNAR", TRUE ~ missingness))  # set unique hits that meet input criteria to MNAR so they get imputed later
  
  
  ## PROTEIN COMPARISON PAIR LEVEL FILTERING !!!!!!!! - based on minimal id requirements (replicate ids control vs. condition) ======================
  data_prot_mis_filtered <- data_prot_mis  %>% 
    filter(!replicate_ids_left_right %in% filter_prot_comparison_pair_completeness)
  cat("COMPARISON PAIR FILTERING:      removing these:      ", filter_prot_comparison_pair_completeness)
  
  
  ## IMPUTATION =======================================================================================================================================
  if(protti_imputation == "no_imputation") {
    data_prot_mis_filtered_imp         <- data_prot_mis_filtered
    data_prot_mis_filtered_imp$imputed <- "no"
    data_prot_mis_filtered_imp$imputed_prot_intensity_log2 <- data_prot_mis_filtered_imp$normalised_protein_intensity_log2
    cat(paste("", " >>>   NO PROTTI IMPUTATION PERFORMED   <<< ", " however downstream variables and colnames still named ...imputed", " ", sep = "\n"))
  } else if (protti_imputation == "ludovic" | protti_imputation == "noise") {
    paste(" >>>    IMPUTATION Method = ", protti_imputation, "   <<<", sep = "")
    data_prot_mis_filtered_imp <- impute(data        = data_prot_mis_filtered,  # output from missingness calculation
                                         sample      = run, 
                                         grouping    = entry,  # precursor level imputation possible too
                                         intensity_log2 = normalised_protein_intensity_log2,
                                         condition   = condition,
                                         comparison  = comparison,
                                         missingness = missingness,
                                         method      = protti_imputation,   # requires MAR NMAR output of assign_missingness comand
                                         skip_log2_transform_error = TRUE,  # to avoid issues for log2 values <20
                                         retain_columns = c(run, condition, bio_replicate, con_rep, entry, entry_name, protein_name, missingness, 
                                                            n_measurments_overall, n_measurments_per_condition, n_measurments_per_sample, 
                                                            n_features_overall, n_features_per_condition, n_features_per_sample,  
                                                            replicate_ids_left, replicate_ids_right, replicate_ids_left_right, n_obs_per_comparison, 
                                                            features_left, features_right, features_left_right_overall)
    ) %>%  rename(imputed_prot_intensity_log2 = imputed_intensity)
  }
  
  data_prot_mis_filtered_imp <- data_prot_mis_filtered_imp %>%
    group_by(entry, comparison) %>%
    mutate(imputed_comparison = ifelse(any(imputed == TRUE), "yes", "no")) %>%
    ungroup()
  
  ## protein ranking ----------------------------------------------------------------------
  data_prot_mis_filtered_imp <- data_prot_mis_filtered_imp %>%   # calculate medians
    group_by(entry) %>% 
    mutate(global_median_imp_log2_intensity    = median(imputed_prot_intensity_log2, na.rm = TRUE)) %>%
    group_by(entry, condition) %>% 
    mutate(condition_median_imp_log2_intensity = median(imputed_prot_intensity_log2, na.rm = TRUE)) %>%
    ungroup() 
  data_prot_mis_filtered_imp <- data_prot_mis_filtered_imp %>%   # global protein rank
    arrange(desc(global_median_imp_log2_intensity)) %>% 
    distinct(entry, global_median_imp_log2_intensity) %>%
    mutate(rank_protein_global = row_number()) %>%
    right_join(data_prot_mis_filtered_imp, by = c("entry", "global_median_imp_log2_intensity")) %>%
    arrange(row_number())
  data_prot_mis_filtered_imp <- data_prot_mis_filtered_imp %>%  # condition protein rank
    distinct(entry, condition, condition_median_imp_log2_intensity) %>%
    group_by(condition) %>% #
    arrange(desc(condition_median_imp_log2_intensity)) %>% 
    mutate(rank_protein_condition = row_number())  %>%
    ungroup() %>%
    right_join(data_prot_mis_filtered_imp, by = c("entry", "condition", "condition_median_imp_log2_intensity")) %>%
    arrange(row_number()) 
  
  # rank plotting
  data_prot_mis_filtered_imp <- data_prot_mis_filtered_imp %>% 
    mutate(meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no"),) %>%
    arrange(meta_surfaceome)
  
  ggplot(data_prot_mis_filtered_imp, 
         aes(x = rank_protein_global, y = global_median_imp_log2_intensity, color = meta_surfaceome)) +
    geom_point() +
    scale_color_manual(values = c("yes" = "#900020", "no" = "#515153")) +
    labs(title = "Global Protein Rank", x = "Rank",  y = "log2(Intensity)", color = "meta_surfaceome") +
    theme_minimal()
  plot_list <-   data_prot_mis_filtered_imp %>%
    group_by(condition) %>%
    group_split() %>%
    map(~ ggplot(., 
                 aes(x = rank_protein_global, y = global_median_imp_log2_intensity, 
                     color = meta_surfaceome,
                     text = protein_name, protein = entry_name)) +
          geom_point(alpha= 0.7) +
          scale_color_manual(values = c("yes" = "#900020", "no" = "#515153")) +
          labs(title = paste("Protein Intensity vs Rank -", .$condition[1]),
               x = "Rank",  
               y = "Intensity", 
               color = "meta_surfaceome") +
          theme_minimal()
    )  
  num_subsets <- length(plot_list)          #  grid dimensions
  num_rows    <- ceiling(num_subsets / 3)   #  grid dimensions
  grid_plot <- do.call(grid.arrange, c(plot_list, ncol = 3, nrow = num_rows))   # Arrange plots in grid
  # rank plot - grid layout
  print(grid_plot)
  
  ## CV per protein and COMPARISON (using condition would yield wrong CVs for ref condition which is duplicated!!) ======================
  data_prot_mis_filtered_imp <- data_prot_mis_filtered_imp %>%
    group_by(entry, condition, comparison) %>% # important to group also based on comparison as imputation might affect signal in one but not the other comparison (e.g. MNAR in one MAR in the other comparison)
    mutate(
      mean_protein = mean(2^imputed_prot_intensity_log2, na.rm = TRUE),
      sd_protein   = sd(2^imputed_prot_intensity_log2, na.rm = TRUE),
      cv_protein   = (sd_protein / mean_protein) * 100 ) %>%
    # dplyr::select(-sd_protein, -mean_protein) %>%
    ungroup()  %>%
    mutate(meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no")) %>% # for visualization
    arrange(meta_surfaceome)  # ensure meta_surfaceome"yes" on top of df and thus top of figures
  
  ## diff_abundance_subset testing. CALCULATE ALL COMPARISON SEPERATE OTHERWISE HIGHER P-VALUES  =================================================================================================================================
  # zero IMPUTE missing values for t-test; If NA remains as intensity the quant info will be lost in t-test ---------------------------------------------------------------------------------------------------------------
  # replicate_ids_right (numeric ranging 1:n(replicates)) refers to number of replicate ids of current row >> using this column intensities are imputed  for all entries and conditions 
  if(protti_imputation == "no_imputation") {  # zero impute to not loose one sided quant info (unique hits)
    diff_calc_input <-  data_prot_mis_filtered_imp %>%
      mutate(ifelse(replicate_ids_right == 0, 0, imputed_prot_intensity_log2)) # zero impute intensity for rows (replicates) where there has been no protein identification (each row is a protein x bio_replicate x condition )
    paste0("performing zero imputation to avoid data loss at t-test step (= removal of unique hits due to NA intensity)")
  } else if (protti_imputation == "ludovic" | protti_imputation == "noise") {  # issue is exclusive to no_imputation
    paste0("no zero imputation required")
    diff_calc_input <-  data_prot_mis_filtered_imp 
  }   
  
  # diff_abundance_subset dataframe - define empty prior loop; serves collection of diff_abundace data
  data_prot_diff_abundance <- data.frame(entry= character(), entry_name= character(), protein_name= character(), n_measurments_overall= numeric(), n_features_overall= numeric(),  
                                         condition_median_imp_log2_intensity= numeric(), rank_protein_condition= numeric(), rank_protein_global= numeric(), comparison= character(), missingness= character(), 
                                         log2FC= numeric(), CI_2.5= numeric(), CI_97.5= numeric(), avg_abundance= numeric(), t_statistic= numeric(), pvalue= numeric(), adj_pvalue= numeric(), B= numeric(), n_obs= numeric(), 
                                         features_left= numeric(), features_right= numeric(), features_left_right_overall= character(), replicate_ids_left= numeric(), replicate_ids_right= numeric(), n_obs_per_comparison= numeric(), 
                                         replicate_ids_left_right= character(), imputed= character(), condition = character(), reference = character(),
                                         imputed_comparison= character() , missing_percentage= numeric(), DF= numeric(), test_condition= character(), ttest_reference= character() 
  )
  
  for (i in unique(diff_calc_input$comparison)) { # diff loop --------------------------------------------------------------------------------------------------------------------------------------------------------------
    print(paste0("calculating diff abundance for: ", i))
    diff_abundance_subset <-  diff_calc_input[diff_calc_input$comparison == i,] %>%
      calculate_diff_abundance(
        filter_NA_missingness = FALSE,    # TRUE is default # If a reference/treatment pair has too few samples to be considered robust based on user defined cutoffs, it is annotated with NA as missingness by the assign_missingness() function. If this argument is TRUE, these NA reference/treatment pairs are filtered out.
        sample      = run,
        condition   = condition,
        grouping    = entry,              # or precursor (if changed - also change to )
        intensity_log2 = imputed_prot_intensity_log2,  # imputed_intensity or normalised_protein_intensity_log2 if no imputation desired
        missingness = missingness,
        comparison  = comparison,          # "comparison" column is output of "assign_missingness" command
        method      = protti_t_testing,    # recommendation:  >>> moderated_t-test <<<  this method is more robust than standard t-tests, especially when dealing with small sample sizes
        retain_columns = c(entry, entry_name, protein_name, replicate_ids_left_right,
                           n_measurments_overall,  
                           n_features_overall,  
                           imputed, condition_median_imp_log2_intensity, rank_protein_condition, rank_protein_global, 
                           replicate_ids_left, replicate_ids_right,  n_obs_per_comparison, 
                           features_left, features_right, features_left_right_overall, imputed_comparison
        )) %>%  
      rename(log2FC = diff,
             pvalue = pval,                               
             adj_pvalue = adj_pval) %>%
      distinct(entry, comparison, .keep_all = TRUE)  %>%   # remove values for replicates (now max 2 rows per protein - 1 for each condition)
      filter(!is.na(replicate_ids_left) & !is.na(replicate_ids_right)) %>%  # kick proteins with no replicate id on either side of comparison
      mutate(missing_percentage = 1 -(replicate_ids_left + replicate_ids_right) / max(n_obs_per_comparison),  # can be used for dataset filtering
             DF = "dummy column") # helper col for MSstats

    # prepare mapping of CV values to diff_abundance_subset data. condition and experiment_references in separate columns. CVs can be used e.g. for unique hit plot y axis spacing (high CVs seem to be associated with LUX targets see PBMC_v24 MartinGesell)
    diff_abundance_subset <- diff_abundance_subset %>%
      separate(comparison,   into = c("ttest_condition", "ttest_reference"), 
               sep = "_vs_", remove = FALSE)
    
    # append each comparison to meta dataframe (data_prot_diff_abundance)
    data_prot_diff_abundance <- rbind(data_prot_diff_abundance, diff_abundance_subset)
  }
  # export diff abundance data for downstream processing
  write.csv(data_prot_diff_abundance,    paste0(output_directory, "/_data_prot_diff_abundance.csv"))
  
  # adjust cv dataframe to make mergable
  cv_protein_merge_helper_df <- data_prot_mis_filtered_imp %>%
    dplyr::select(entry, condition, comparison, cv_protein) %>%
      separate(comparison,   into = c("ttest_condition", "ttest_reference"), 
               sep = "_vs_", remove = FALSE) %>%
    distinct() %>%  # Remove duplicate rows (replicates!)
    # work with each unique entry-comparison pair
    group_by(entry, comparison) %>%
    # summarize the data to create our desired output
    summarize(
      ttest_condition = first(ttest_condition),    # Keep the first occurrence of ttest_condition and ttest_reference (they should be the same for each group anyway)
      ttest_reference = first(ttest_reference),
      cv_protein_condition = cv_protein[condition == ttest_condition],  # For cv_condition, get the cv_protein value where condition matches ttest_condition
      cv_protein_ref = cv_protein[condition == ttest_reference],
      .groups = "drop"   # Drop the grouping to avoid any issues in further operations
    ) %>%
    dplyr::select(-ttest_condition, -ttest_reference)
    
  # include CV in diff abundance data  
  data_prot_diff_abundance <- data_prot_diff_abundance %>%
    left_join(cv_protein_merge_helper_df,
              by = c("entry" = "entry", "comparison" = "comparison"))


  ## depending on imputation done or not done there will be different sets that require Inf assignment ...
  if(protti_imputation == "no_imputation") {  # no imputation >> set Infitine and mask pvalue for all unique hits that meat filter_uni_hit_min_completeness requirement
    # set +/-Inf FCs & overwrite p-values
    data_prot_diff_abundance <- data_prot_diff_abundance %>% 
      mutate(log2FC = case_when(replicate_ids_left  == 0  &  replicate_ids_right/max(replicate_ids_right)  >= filter_uni_hit_min_completeness   ~ +Inf,    # if e.g in 2/3 or 2/4 ids in treated while 0 ids in ctrl --> unique right
                                replicate_ids_right == 0  &  replicate_ids_left/max(replicate_ids_left)    >= filter_uni_hit_min_completeness   ~ -Inf,    # if e.g in 2/3 or 2/4 ids in control while 0 ids in treated --> unique left
                                TRUE ~ log2FC)) %>%
      mutate(pvalue = ifelse(is.infinite(log2FC), NA, pvalue))  # critical to overwrite pvalues - for non-complete ided unique hits there are pvalues >> creates issue at subset_inf step in loop
  } else if (protti_imputation == "ludovic" | protti_imputation == "noise") {  # imputation >> set inf only those that have 2_vs_0 or 0_vs_2 replicate_ids_left_right (otherwise they will not be captured in inf plots (Inf intensity required to be assigned to that subset_if)
    # set FCs for proteins that are currently NA but make the filter_uni_hit_min_completeness cutoff to +/-Inf   (so they end up in subset_inf in loop)
    data_prot_diff_abundance <- data_prot_diff_abundance %>% 
      mutate(log2FC = case_when(is.na(log2FC) & replicate_ids_left  == 0  &  replicate_ids_right/max(replicate_ids_right)  >= filter_uni_hit_min_completeness   ~ +Inf,    # if e.g in 2/3 or 2/4 ids in treated while 0 ids in ctrl --> unique right
                                is.na(log2FC) & replicate_ids_right == 0  &  replicate_ids_left/max(replicate_ids_left)    >= filter_uni_hit_min_completeness   ~ -Inf,    # if e.g in 2/3 or 2/4 ids in control while 0 ids in treated --> unique left
                                TRUE ~ log2FC))
    # no need to overwrite pvalues as in if {} above - they are already NA where FC is NA
  }
  
  # assign imputation 
  data_prot_diff_abundance <- data_prot_diff_abundance %>%
    dplyr::select(-imputed) %>%   # when protein df is collapsed to comparison df the imputation column looses its value (retained but does not make sense anymore)
    group_by(comparison) %>% 
    mutate(imputed_left  = ifelse(imputed_comparison == "yes" &  replicate_ids_left  < max(replicate_ids_left ), "yes", "no"),  # if imputation & less than max replicate ids left consider as imputed left
           imputed_right = ifelse(imputed_comparison == "yes" &  replicate_ids_right < max(replicate_ids_right), "yes", "no"), # if imputation & less than max replicate ids right consider as imputed right
           imputed_left_right = case_when(imputed_left  == "yes" & imputed_right == "no" ~ "left",
                                          imputed_right == "yes" & imputed_left  == "no" ~ "right",
                                          imputed_left  == "yes" & imputed_right == "yes"  ~ "both",
                                          TRUE ~ "none" )) %>%  # to handle cases where neither is TRUE
    ungroup()
  
  # for compatibility with downstream code save results to required variable names
  data_level_feature <- data_prec_norm_feature_filtered    # NON IMPUTED; precursor level     normalized, glyco info containing, non zero, if peptideCSC workflow gylco filtered, "filter_exclude_protein" filtered data, non-imputed
  data_level_prot    <- data_prot_mis_filtered_imp %>%                # IMPUTED;     protein level       normalized, missingness, imputed according to "imputation", protein level filtered 
    filter(!is.na(imputed_prot_intensity_log2))  # kick al NA values
  
  pval_distribution_plot(                # p-val QC plot
    data     = data_prot_diff_abundance,
    grouping = protein_name,
    pval     = pvalue     )
  # END protti processing block ________________________________________________________________________________________________________________________________________________________________________________
  
  
} 
write.csv(data_level_prot,    paste0(output_directory, "/_data_prot_level.csv"))
write.csv(data_level_feature, paste0(output_directory, "/_data_feature_level.csv"))
data_median_intensity_protein <- median(data_level_prot$condition_median_imp_log2_intensity)
data_median_intensity_peptide <- median(data_level_feature$condition_median_log2_feature_intensity)

## Protti QC plotting ---------------------------------------------------------------------------------------------------------------------------------------------
if (data_processing == "protti"){          #protti QC block
  qc_plot_protti_1 <- qc_ranked_intensities(    # global protein rank
    data     = data_prot_mis_filtered_imp,                  # data_prot_mis_filtered_imp              or     data_prec_norm_feature_filtered
    sample   = con_rep,
    grouping = entry,                              # entry                          or     peptide_sequence_mod
    intensity_log2 = imputed_prot_intensity_log2,  # imputed_prot_intensity_log2    or     normalised_prec_intensity_log2
    facet    = FALSE, # logical value that specifies whether the calculation should be done group wise by sample and if the resulting plot should be faceted by sample. (default is FALSE). If facet = FALSE the median of each protein intensity will be returned.
    plot     = TRUE,
    interactive = FALSE
  ) 
  qc_plot_protti_2 <- qc_intensity_distribution(
    data        = data_prec_norm_feature_filtered,
    grouping    = peptide_sequence_mod,
    intensity_log2 = normalised_prec_intensity_log2,
    plot_style  = "histogram"
  )
  ##### missed cleavages
  trypsin_pattern <- "[KR][^P]"  # trypsin pattern
  data_prec_norm_feature_filtered$n_missed_cleavage <- str_count(data_prec_norm_feature_filtered$stripped_peptide_sequence, trypsin_pattern)
  qc_plot_protti_3 <- qc_missed_cleavages(          # missed cleavages
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    grouping  = peptide_sequence_mod,
    missed_cleavages = n_missed_cleavage,
    method    = "intensity",
    intensity = normalised_prec_intensity,
    plot      = TRUE
  )
  # proteinNameFilterExclude  = c("contam_","iRT", "_YEAST")
  contaminant_mask <- grepl(paste(filter_exclude_protein, collapse="|"), data_raw$PSM_id)
  data_raw$is_contaminant <- ifelse(contaminant_mask, "TRUE", "FALSE")
  data_raw <- data_raw %>% 
    mutate(con_rep = paste(condition, bio_replicate, sep="_")) #%>% 
  # group_by(run, entry) %>% 
  # mutate(cont_signal_per_run = sum(raw_prec_intensity)) %>%
  # arrange(desc(cont_signal_per_run)) %>%
  # ungroup()
  qc_plot_protti_4 <- qc_contaminants(
    data        = data_raw,
    sample      = con_rep,
    protein     = entry,
    is_contaminant = is_contaminant, # a logical column that indicates if the protein is a contaminant.
    intensity   = raw_prec_intensity,
    n_contaminants = 5,              # a numeric value that indicates how many contaminants should be displayed individually. The rest is combined to a group called "other". The default is 5.
    plot        = TRUE,
    interactive = FALSE
  )
  qc_plot_protti_4 <- qc_plot_protti_4 + theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 10), 
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "inches")
  )
  qc_plot_protti_5 <- qc_intensity_distribution(    # Raw Intensity Boxplot
    data     = data_prec_norm_feature_filtered,
    sample   = con_rep,
    grouping = precursor,
    intensity_log2 = raw_prec_intensity_log2,
    plot_style = "boxplot"
  )
  qc_plot_protti_6 <- qc_intensity_distribution(    # Normalized Intensity Boxplot
    data     = data_prec_norm_feature_filtered,
    sample   = con_rep,
    grouping = precursor,
    intensity_log2 = normalised_prec_intensity_log2,
    plot_style = "boxplot"
  )
  qc_plot_protti_7 <- qc_median_intensities(        # median raw signal
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    grouping  = peptide_sequence_mod,
    intensity = raw_prec_intensity_log2
  )
  qc_plot_protti_8 <- qc_ids(                       # Protein / sample plot
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    grouping  = protein_name,
    intensity = normalised_prec_intensity,
    condition = condition,
    title     = "Protein identifications per sample",
    plot      = TRUE
  )
  qc_plot_protti_9 <- qc_cvs(                       # CV plot overall
    data      = data_prec_norm_feature_filtered,
    grouping  = precursor,
    condition = condition,
    intensity = normalised_prec_intensity,
    plot = TRUE,
    plot_style = "violin"
  )
  qc_plot_protti_10 <- qc_cvs(                       # CV plot surface
    data      = data_prec_norm_feature_filtered[data_prec_norm_feature_filtered$meta_surfaceome == "yes",],
    grouping  = precursor,
    condition = condition,
    intensity = normalised_prec_intensity,
    plot = TRUE,
    plot_style = "violin"
  )
  qc_plot_protti_11 <- qc_charge_states(             # % signal contribution of charge states
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    grouping  = peptide_sequence_mod,
    charge_states = precursor_charge,
    method    = "intensity",
    intensity = normalised_prec_intensity,
    plot      = TRUE
  )
  qc_plot_protti_12 <- qc_peptide_type(
    method    = "intensity",       #  % raw signal derived from ...-tryptic peptides 
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    peptide   = stripped_peptide_sequence,
    pep_type  = pep_type,
    intensity = normalised_prec_intensity,    
    plot      = TRUE,
    interactive = FALSE
  )
  qc_plot_protti_13 <- qc_peptide_type(    
    method    = "count",          #  count of ...-tryptic peptides 
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    peptide   = stripped_peptide_sequence,
    pep_type  = pep_type,
    intensity = normalised_prec_intensity,   
    remove_na_intensities = TRUE,
    plot      = TRUE,
    interactive = FALSE
  )
  qc_plot_protti_14 <- qc_sequence_coverage(        # sequence coverage
    data               = data_prec_norm_feature_filtered,
    protein_identifier = protein_name,
    coverage           = coverage
  )
  qc_plot_protti_15 <- qc_data_completeness(
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    grouping  = peptide_sequence_mod,
    intensity = normalised_prec_intensity,
    plot      = TRUE
  )
  if(data_type == "DIA") {  # DDA data currently does not contain rt info
    qc_plot_protti_16 <-   qc_peak_width(   
      data = data_raw,
      sample = run,
      intensity = raw_prec_intensity,
      retention_time = rt,
      peak_width = NULL,
      retention_time_start = rt_start,
      retention_time_end = rt_stop,
      remove_na_intensities = TRUE,
      interactive = FALSE
    ) 
  } else {
    qc_plot_protti_16 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "currentley peak width not reported in DDA output", size = 4) +
      theme_few()    
  }
  qc_plot_protti_17 <- qc_pca(
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    grouping  = peptide_sequence_mod,
    intensity = normalised_prec_intensity_log2,
    condition = condition,
    digestion = NULL,
    plot_style = "scree"
  )
  qc_plot_protti_18 <- qc_pca(
    data      = data_prec_norm_feature_filtered,
    sample    = con_rep,
    grouping  = peptide_sequence_mod,
    intensity = normalised_prec_intensity_log2,
    condition = condition,
    components = c("PC1", "PC2"),
    plot_style = "pca"
  )
  # # qc_proteome_coverage plot takes forever - consider to exclude
  # qc_plot_protti_19 <- qc_proteome_coverage(
  #   data        = data_prec_norm_feature_filtered,
  #   sample      = con_rep,
  #   protein_id  = protein_name,
  #   organism_id = 9606,      # a numeric value that specifies a NCBI taxonomy identifier (TaxId) of the organism used. Human: 9606, S. cerevisiae: 559292, E. coli: 83333.
  #   reviewed    = TRUE,
  #   plot        = TRUE,
  #   interactive = FALSE
  # )
  qc_plot_protti_19 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "qc_proteome_coverage currently disabled", size = 4) +
    theme_few()    
  qc_plot_protti_20 <- qc_sample_correlation(
    data        = data_prec_norm_feature_filtered,
    sample      = con_rep,
    grouping    = precursor,
    intensity_log2 = normalised_prec_intensity_log2,
    condition   = condition,
    interactive = FALSE
  )
  setwd(paste(output_directory))   
  # pdf
  pdf(paste0(output_directory, "/_qc_protti.pdf"), onefile = TRUE, width = 10, height = 10)
  for (i in 1:20) {
    current_plot <- get(paste0("qc_plot_protti_", i))
    print(current_plot)
    if (i == 19) {  # new page after plot 19 - not clear why but otherwise plot 19 and 20 overlap.
      grid::grid.newpage()  
    }
  }
  dev.off()
  # png composite 
  plot_list <- list(qc_plot_protti_1, qc_plot_protti_2, qc_plot_protti_3, qc_plot_protti_4, qc_plot_protti_5, qc_plot_protti_6, qc_plot_protti_7, qc_plot_protti_8, qc_plot_protti_9, qc_plot_protti_10, 
                    qc_plot_protti_11, qc_plot_protti_12, qc_plot_protti_13, qc_plot_protti_14, qc_plot_protti_15, qc_plot_protti_16, qc_plot_protti_17, qc_plot_protti_18, qc_plot_protti_19)
  plot_list[[20]] <- ggdraw() + draw_grob(grid::grid.grabExpr(print(qc_plot_protti_20, width = 6, height = 4)))
  png("_qc_protti_composite.png", width = 30, height = 16, units = "in", res = 150)
  print(plot_grid(plotlist = plot_list, ncol = 4))
  dev.off()
  rm(plot_list)
  
  # tryptic peptides
  grouped_data <- data_prec_norm_feature_filtered %>%
    group_by(pep_type, meta_surfaceome) %>%
    summarise(count = n()) %>%
    ungroup()
  # Create the stacked bar plot with custom colors
  ggplot(grouped_data, aes(x = pep_type, y = count, fill = meta_surfaceome)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("yes" = "red", "no" = "grey")) +
    labs(title = "Precursor Count by Peptide Type and Surfaceome",
         x = "Peptide Type",
         y = "Number of Precursors",
         fill = "Surfaceome") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right")
  
  grouped_data <- data_prec_norm_feature_filtered %>%
    group_by(pep_type, meta_surfaceome) %>%
    summarise(count = n()) %>%
    group_by(pep_type) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # Create the stacked bar plot with percentages and custom colors
  ggplot(grouped_data, aes(x = pep_type, y = percentage, fill = meta_surfaceome)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("yes" = "red", "no" = "grey")) +
    labs(title = "Percentage of Surfaceome by Peptide Type",
         x = "Peptide Type",
         y = "Percentage of Total",
         fill = "Surfaceome") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right") +
    scale_y_continuous(labels = scales::percent_format(scale = 1))
  
  
  
  ## CV qc plots ===============================================================================================================================
  # qc plot for precursor level CVs  
  qc_precursor_cv <- ggplot(data_prec_norm_feature_filtered %>% distinct(precursor, run, .keep_all = TRUE),  # one datapoint per entry and comparison (otherwise replicates will appear as overlaying dots - not bad but unnecessary)
                            aes(x=condition_median_log2_feature_intensity, y=cv_precursor, color = meta_surfaceome, text = entry_name, feature = precursor)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("yes" = "#900020", "no" = "#5F5F61")) +  
    theme_minimal() +
    labs(x = "condition_median_log2_feature_intensity", y = "cv_precursor") +
    facet_wrap(~ condition, labeller = labeller(condition = label_value))  # Create separate plots for each condition with custom labels
  qc_precursor_cv_i <- ggplotly(qc_precursor_cv, tooltip = c("text", "feature")) # interactive plot
  
  # qc for protein level CVs
  qc_prot_cv <- ggplot(data_prot_mis_filtered_imp %>% distinct(entry, entry_name, comparison, .keep_all = TRUE),  # one datapoint per entry and comparison (otherwise replicates will appear as overlaying dots - not bad but unnecessary)
                       aes(x=condition_median_imp_log2_intensity, y=cv_protein, color = meta_surfaceome, text = entry_name, conditions = comparison)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("yes" = "#900020", "no" = "#5F5F61")) +  
    theme_minimal() +
    labs(x = "imputed_median_prot_intensity_log2", y = "cv_protein") +
    facet_wrap(~ condition, labeller = labeller(condition = label_value))  # Create separate plots for each condition with custom labels
  qc_prot_cv_i <- ggplotly(qc_prot_cv, tooltip = c("text", "conditions")) # interactive plot
  
  saveWidget(qc_precursor_cv_i, file= paste(experiment_name, "_qc_precursor_cv_i_", ".html", sep="")) # save interactive plot as HTML
  saveWidget(qc_prot_cv_i,      file= paste(experiment_name, "_qc_protein_cv_i_",   ".html", sep="")) # save interactive plot as HTML
  
} ## END Protti QC plotting _________________________________________________________________________________________________________________________________________
################################################################################################################################################################# 
################################################################################################################################################################# 
################################################################################################################################################################# 


# remove not needed data and safe workplace 
rm("current_plot", "data_CSPA_fraction", "data_reps", "expected_proteins", "qc_precursor_cv", "qc_precursor_cv_i", "qc_prot_cv", "qc_prot_cv_i",              # remove stuff
   "contaminant_mask", "fraction", "i", "includeFeature", "limit", "proportions_cellType", "remainder", "reps", "uniprot_ids", "create_barplot_with_numbers", # remove values
   "diff_calc_input", "glycoAnnotDf", "glycoAnnotDf_L", "grid_plot", "grouped_data", "cell_specific_CSPA", "data_raw_pre_filtered", "data_raw_up_surf_gly",          # remove intermediates
   "data_raw_up_surf_gly_nZ", "data_prec_norm", "data_prot_mis", "data_prot_mis_filtered", "data_prot_missing") 
rm(list = ls(pattern = "^qc_plot_protti|^qc_plot_surface")) # remove qc plots)

save.image(file = paste0(output_directory, "/shs_data_backup_before_loop.RData"))
# load(paste0(output_directory, "/shs_data_backup_before_loop.RData"))
 

################################################################################################################################################################# 
################################################################################################################################################################# 
###### LOOOOOOOP ################################################################################################################################################ 
#### Predefine dataframe used to aggregate Candidate-Protein (CandiProt = SigUp/SigDown/UniqueUp/UniqueDown) information across loop:
super_volcano_candi_prot_L <- data.frame(condition =character(), comparison = character(), entry = character(), entry_name = character(), protein_name= character(), log2FC=numeric(), aggreg_unique_hits=character(), adj_pvalue = numeric(), pvalue=numeric(), #q_value = numeric(), 
                                         features_left= numeric(), features_right = numeric(), features_left_right_overall= numeric(), replicate_ids_left= numeric(), replicate_ids_right= numeric(), replicate_ids_left_right= numeric(),
                                         cspa_2015= character(), cd_antigen= character(), surfy_2018= character(), meta_surfaceome =character(), poi= character(), 
                                         condition_median_imp_log2_intensity = numeric(), imputed_left_right = character(), string_interactor =character() )

loop_frame <- unique(data_prot_diff_abundance$comparison)

## for loop testing excecute:        comp_counter <- 1    
for (comp_counter in 1:length(loop_frame)) {   # ----------------------- LOOOOOOOP --------------------- LOOOOOOOP ------------------------------ LOOOOOOOP ------------------- LOOOOOOOP -------------------
  setwd(output_directory)
  
  print("protti loop")
  usedConditions <- c(gsub("_vs_.*", "", loop_frame[comp_counter]), gsub("^.*_vs_", "", loop_frame[comp_counter]))
  result_subset <- data_prot_diff_abundance[data_prot_diff_abundance$comparison == loop_frame[comp_counter],] 
  
  # adjust p value (<<-- after this point no more data filtering)
  result_subset$adj_pvalue = p.adjust(result_subset$pvalue, method="BH")
  # introduce qvalue (accodring to JS more suited for omics data than adjusted pvalue)
  # result_subset$q_value <- qvalue(result_subset$pvalue)[3] %>% unlist()
  # # selective qvalue - only for proteins crossing FC cutoff
  # somehow this does not work although pvals are in range 0:1 other 
  # result_subset$q_value_FC_subset <- NA
  # result_subset[abs(result_subset$log2FC) >= filter_log2fc_cutoff,]$q_value_FC_subset <- 
  #   qvalue(result_subset[abs(result_subset$log2FC) >= filter_log2fc_cutoff, ])[3] %>% unlist()
  # set the minimal p-values -- no "zero" p-values in output
  
  pval_distribution_plot(                # p-val QC plot
    data     = result_subset,
    grouping = protein_name,
    pval     = pvalue     )
  pval_distribution_plot(                # adj_pval QC plot
    data     = result_subset,
    grouping = protein_name,
    pval     = adj_pvalue )
  # pval_distribution_plot(                # q-val QC plot
  #   data     = result_subset,
  #   grouping = protein_name,
  #   pval     = q_value )
  
  # trim significance values 
  result_subset$pvalue      [result_subset$pvalue        < 2.220446e-16 ] = 2.220446e-16
  result_subset$adj_pvalue  [result_subset$adj_pvalue    < 2.220446e-16 ] = 2.220446e-16
  # result_subset$q_value     [result_subset$q_value       < 2.220446e-16 ] = 2.220446e-16
  #result_subset$q_value_FC_subset[result_subset$q_value_FC_subset < 2.220446e-16 ] = 2.220446e-16
  result_subset <- result_subset[order(result_subset$adj_pvalue),]  # sort by significance - good on top
  ## Volcano Plot with unique hit display
  result_subset <- result_subset[!is.na(result_subset$log2FC), ]   # filter out all proteins that were not identified - all other should have calculated FC or +/-Inf assigned
  
  for(i in poi_plotting_lists) { #  ----------------------------------------------------------------------------------------------------------------------------------------------------------  
    # define current loop pois
    poi_var = as.vector(poi_table[[i]][poi_table[[i]] != ""])
    # create directory where current poi-plots are to be saved
    if (!dir.exists(paste0(input_directory, "/_output", test_parameter, script_version, "/poi_", i))) {
      dir.create(paste0(input_directory, "/_output", test_parameter, script_version, "/poi_", i), recursive = TRUE)
    }
    setwd(          paste0(input_directory, "/_output", test_parameter, script_version, "/poi_", i)) # new directory for each poi category
    
    # Introduce surfaceome annotations & poi info
    result_subset <- result_subset %>%
      mutate(cspa_2015       = ifelse(entry %in% surface_annotations$cspa_2015      , "yes", "no"),
             surfy_2018      = ifelse(entry %in% surface_annotations$surfy_2018     , "yes", "no"),
             tcsa_2021       = ifelse(entry %in% surface_annotations$tcsa_2021      , "yes", "no"),
             meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no"),
             cd_antigen      = ifelse(entry %in% surface_annotations$cd_antigen     , "yes", "no"),
             uniprot_2023    = ifelse(entry %in% surface_annotations$uniprot_2023   , "yes", "no"),
             poi   =      ifelse(entry_name %in% poi_var                            , "yes", "no")    # pois are given as entry_names in MGs poi meta file
      )  
    # summarize surfaceome categories present in dataset
    data.frame(Category = c("cspa_2015", "surfy_2018", "tcsa_2021","uniprot_2023" ,"cd_antigen" , "meta_surfaceome", "poi"),
               Total    = replicate(7, sum(str_count(result_subset$cd_antigen, "yes|no"))),
               yes      = c(sum(str_count(result_subset$cspa_2015,    "yes")),
                            sum(str_count(result_subset$surfy_2018,   "yes")),
                            sum(str_count(result_subset$tcsa_2021,    "yes")),
                            sum(str_count(result_subset$uniprot_2023, "yes")),
                            sum(str_count(result_subset$cd_antigen,   "yes")),
                            sum(str_count(result_subset$meta_surfaceome, "yes")),
                            sum(str_count(result_subset$poi, "yes")) 
               )
    )
    
    print("---- Done:  df subsetting, SuperVolcanoData & Annotations done  ----------------------------------------------------") 
    
    # string interactor mapping
    string_interactors <- string[string$protein1_entry_name %in% string_targets[[ experiment_conditions[[comp_counter]] ]],] %>% 
      dplyr::select(protein2_entry_name) %>% 
      unique()
    result_subset <- result_subset %>% # map interactors to data 
      mutate(string_interactor = case_when(entry_name %in% string_interactors[[1]] ~ "yes", TRUE ~ "no") ) 
    # ________________________________________________________________________________________________________________________________________________________________
    
    ## dataframe for plotting purposes - no more data filtering below this point
    # aggregate SuperVolcanoPlot Information in one dataframe (that can be exported and easily checked for GO enrichment ...)
    super_volcano_data <- result_subset
    
    # specify plot dot colors for suvo
    super_volcano_data  <-  super_volcano_data %>%
      mutate(suvo_plot_color = case_when(
        super_volcano_data[[plot_fill_volcano]] == "yes" & string_interactor == "yes" ~ "surface+string",    # is both plot_fill_volcano and string_interactor
        super_volcano_data[[plot_fill_volcano]] == "yes" ~ "surface",                                        # is plot_fill_volcano
        string_interactor == "yes" ~ "string",                                                               # is string_interactor
        TRUE ~ "other"                                                                                       # is neither 
      ))
    
    # specify plot dot occupacy (incomplete dataframe make plot dots a bit transparent --> (alpha) column - adjustet_t-test does not result in missing df proteins. define here how they are displayed in plots
    super_volcano_data <- super_volcano_data %>% # recommended to readjust below part for datasets that contains replicate number > 4
      mutate(plot_alpha = case_when(replicate_ids_left_right %in% c("1_vs_1", "1_vs_2", "2_vs_1", "0_vs_1", "1_vs_0", "2_vs_0", "0_vs_2") ~ plot_no_df_occupacy,  TRUE ~ 1.0 )) %>%
      # highlight data point categories (pvalue defined at script start. specifies pvalue or adj_pvalue to be used for subsetting & "yes" tagging)
      mutate(sig_left  = ifelse(!is.na(.data[[pvalue]]) &               
                                  .data[[pvalue]] < filter_sig_cutoff &
                                  log2FC < -filter_log2fc_cutoff, "yes", 
                                ifelse(is.na(.data[[pvalue]]), "", "")),
             sig_right = ifelse(!is.na(.data[[pvalue]]) & 
                                  .data[[pvalue]] < filter_sig_cutoff & 
                                  log2FC > +filter_log2fc_cutoff, "yes", 
                                ifelse(is.na(.data[[pvalue]]), "", ""))) %>%
      # aggregate testing information into categories
      mutate(no_df_left  = ifelse(is.nan(.data[[pvalue]]) & log2FC < -filter_log2fc_cutoff, "yes", ""),     # noDF
             no_df_right = ifelse(is.nan(.data[[pvalue]]) & log2FC > +filter_log2fc_cutoff, "yes", ""),     # noDF
             unique_hit_left  = ifelse(log2FC == "-Inf", "yes", ""),  # Infinite FC
             unique_hit_right = ifelse(log2FC ==  "Inf", "yes", ""),
             aggreg_left  = ifelse(sig_left  == "yes" |  no_df_left == "yes" |  unique_hit_left == "yes" , "yes", ""),     # dump Left and Right info into individual columns for quick and dirty analysis
             aggreg_right = ifelse(sig_right == "yes" | no_df_right == "yes" | unique_hit_right == "yes" , "yes", ""),
             aggreg_unique_hits = case_when(        #
               unique_hit_left  == "yes" ~ "-", 
               unique_hit_right == "yes" ~ "+",
               TRUE ~ ""),
             left_right  = ifelse(aggreg_right  == "yes" |  aggreg_left == "yes" , "yes", "")) 
    
    super_volcano_data$protein_name <- substr(super_volcano_data$protein_name, 1, 80)  # truncate to allow nice displaying in interactive volcano plots
    # add uniprot infos to __superVolcanoData__ export df 
    super_volcano_data <- merge(super_volcano_data, proteome[,c("entry", "function_cc", "gene_ontology_molecular_function", "glycosylation")], by = "entry")
    
    # order data. 1) meta_surfaceome 
    super_volcano_data <- super_volcano_data %>% 
      arrange(meta_surfaceome) # surface hits on bottom of df >> plotted last = on top
    
    ## For Plotting dataset is separated into different cases
    subset_volcano   <- super_volcano_data[!is.na(super_volcano_data[[pvalue]]), ]    # 1) Fold-change and p value are available -> normal Volcano plot
    subset_inf       <- super_volcano_data[is.infinite(super_volcano_data$log2FC), ]  # 2) unique hits (infinite fold change) -> 1D graph below Volcano plot nad fold-change plot
    subset_no_df     <- super_volcano_data[super_volcano_data$DF == 0 & !is.na(super_volcano_data$DF),  ]  # 3) Fold-change but to few degrees of freedom for p-value calculation -> 1D FC graph below Volcano plot
    # What remains is proteins that could not be quantified 
    # Quick checksum to see all entries were distributed (following expression should come out as true)
    if (nrow(super_volcano_data) == nrow(subset_volcano) + nrow(subset_inf) + nrow(subset_no_df)){
      print("----- All quantified proteins have been assigned to a subset, go ahead ----------")
    } else {
      warning(" !!!!!!  There is an issue with proteins not being assigned to a subset  !!!!!!  ")
    }
    
    # Plotting ==================================================================================================================================================================================
    # Getting fold-change max for x-axis
    foldchangelimit <- max(abs(c(subset_volcano$log2FC, subset_no_df$log2FC)))+0.5
    
    subset_volcano <- subset_volcano %>%
      mutate(significance = as.factor(ifelse(-log10(.data[[pvalue]]) >= -log10(filter_sig_cutoff)   &    abs(log2FC) >= filter_log2fc_cutoff,   "yes",   "no")),  # Mark significant hits
             plot_label   = ifelse(.data[[plot_label_volcano]] == "yes", gsub("_HUMAN", "", entry_name), "")) # Assign plot label 
    
    # Helper frame for significance lines
    segmentation <- data.frame(x = c(-foldchangelimit, filter_log2fc_cutoff, -filter_log2fc_cutoff, filter_log2fc_cutoff),
                               y = c(-log10(filter_sig_cutoff), -log10(filter_sig_cutoff), -log10(filter_sig_cutoff), -log10(filter_sig_cutoff)),
                               xend = c(-filter_log2fc_cutoff, foldchangelimit, -filter_log2fc_cutoff, filter_log2fc_cutoff),
                               yend = c(-log10(filter_sig_cutoff), -log10(filter_sig_cutoff), max(-log10(subset_volcano[[pvalue]])), max(-log10(subset_volcano[[pvalue]]))),
                               col = rep("black", times=4),
                               linetype = rep("dashed", times=4))
    
    # Volcano Plot Subplot 1/3  -------------------------------------------------------------------------------------------------------------------------------
    plot_volcano1 <- ggplot(subset_volcano,
                            aes(x=log2FC, y=-log10(.data[[pvalue]]), label=plot_label, 
                                text= paste0("Protein: <b>", entry_name, "</b><br>", "Protein Name: <b>", protein_name, "</b><br>", "Meta Surfaceome: <b>", 
                                             meta_surfaceome, "</b><br>", "Replicate IDs: ", replicate_ids_left_right, "<br>" , "Feature IDs : ", features_left_right_overall, "<br>" ))
    ) +
      theme_few() +
      geom_segment(x=segmentation$x[1], y=segmentation$y[1], xend=segmentation$xend[1], yend=segmentation$yend[1], linetype="dashed", col="darkgrey", linewidth = 0.4) +
      geom_segment(x=segmentation$x[2], y=segmentation$y[2], xend=segmentation$xend[2], yend=segmentation$yend[2], linetype="dashed", col="darkgrey", linewidth = 0.4) +
      geom_segment(x=segmentation$x[3], y=segmentation$y[3], xend=segmentation$xend[3], yend=segmentation$yend[3], linetype="dashed", col="darkgrey", linewidth = 0.4) +
      geom_segment(x=segmentation$x[4], y=segmentation$y[4], xend=segmentation$xend[4], yend=segmentation$yend[4], linetype="dashed", col="darkgrey", linewidth = 0.4) +
      xlim(-foldchangelimit, foldchangelimit) +
      geom_point(fill = c("#cc66ff", "#cc0000", "#0000EE", "#5F5F61")[match(subset_volcano$suvo_plot_color, c("surface+string", "surface", "string", "other"))],
                 shape = ifelse(subset_volcano$imputed_comparison == "yes", 23, 21),
                 alpha = subset_volcano$plot_alpha,
                 size = 3,
                 color = "white",
                 stroke = 0.4) +
      geom_label_repel(box.padding = 2, size=6, max.overlaps = 100000) +
      labs(fill="meta_surfaceome") +
      theme(axis.title.x = element_blank(),
            legend.position = "top",
            text = element_text(size = 14),
            axis.text = element_text(size = 12)) +
      ylab(paste("-log10(", pvalue, ")", sep ="")) +
      geom_label(label= usedConditions[2],
                 x=-foldchangelimit+0.25, y=0,
                 label.padding = unit(0.8, "lines"),
                 label.size = 1.5,
                 color = "white",
                 fill="black",
                 size = 6) +
      geom_label(label= usedConditions[1],
                 x=foldchangelimit-0.25, y=0,
                 label.padding = unit(0.8, "lines"),
                 label.size = 1.5,
                 color = "white",
                 fill="black",
                 size = 6)
    
    # Subplot Unique IDs ----------------------------------------------------------------------------------------------------------
    # collect protein info, subset df, calculate plot info columns
    prot_level_data_uni <- data_level_prot[data_level_prot$entry %in% unique(subset_inf$entry), ] %>%  # retrieve condition quant data for unique hits
      filter(condition== usedConditions[1] | condition== usedConditions[2]) %>%
      filter(!is.na(imputed_prot_intensity_log2) & !is.na(run))   # remove proteins that were not quantified in current conditions
    # plot_label and data sort by plot_label
    subset_inf <- subset_inf %>% 
      mutate(plot_label   = ifelse(.data[[plot_label_volcano]] == "yes", gsub("_HUMAN", "", entry_name), ""))
    # Split unique hit dataset in - and + infinity
    subset_inf_neg <- subset_inf[subset_inf$log2FC==-Inf,]
    subset_inf_pos <- subset_inf[subset_inf$log2FC==Inf,]
    # Identify axis limits and global intensity median
    intensity_limits <- c(min(subset_inf$condition_median_imp_log2_intensity)-0.25, max(subset_inf$condition_median_imp_log2_intensity)+0.25)
    intensity_median <- median(result_subset$condition_median_imp_log2_intensity)
    cv_limits        <- c(0, log(max(c(subset_inf$cv_protein_condition, subset_inf$cv_protein_ref), na.rm = TRUE)+0.25 , 10))
    
    # Subplot 2/3 neg Unique IDs ---------------------------------------------------------------------------------------------------
    plot_negative_unique_hits <- ggplot(subset_inf_neg,
                                        aes(x = condition_median_imp_log2_intensity, y = log(cv_protein_ref, 10), label=plot_label, text = paste0("Protein: <b>", entry_name, "</b><br>", "Protein Name: <b>", protein_name, "</b><br>", "Meta Surfaceome: <b>", meta_surfaceome, "</b><br>", "Replicate IDs: ", replicate_ids_left_right, "<br>", "Feature IDs : ", features_left_right_overall, "<br>"))
    ) +
      theme_few() +
      geom_point(fill = c("#cc66ff", "#cc0000", "#0000EE", "#5F5F61")[match(subset_inf_neg$suvo_plot_color, c("surface+string", "surface", "string", "other"))],
                 shape = ifelse(subset_inf_neg$imputed_comparison == "yes", 23, 21),
                 alpha = subset_inf_neg$plot_alpha,
                 size = 3,
                 color = "white",
                 stroke = 0.4) +
      geom_segment(x = -0.3, xend = 0.3, y = intensity_median, yend = intensity_median, col = "darkgrey", linetype = "solid", linewidth = 2) +
      geom_label_repel(box.padding = 2, size=6, max.overlaps = 100000) +
      ylim(cv_limits) +
      xlim(intensity_limits) +
      theme(legend.position = "none",
            text = element_text(size = 14),
            axis.text = element_text(size = 12)) +
      ylab("log10(cv_protein)") +
      xlab("log2(Intensity)") +
      geom_segment(x=data_median_intensity_protein, y=-1, xend=data_median_intensity_protein, yend=100, linetype="dashed", col="darkgrey", linewidth = 0.4)  # add median signal line for orientation 
    
    plot_positive_unique_hits <- ggplot(subset_inf_pos,
                                        aes(x=condition_median_imp_log2_intensity, y=log(cv_protein_condition, 10), label=plot_label, text = paste0("Protein: <b>", entry_name, "</b><br>", "Protein Name: <b>", protein_name, "</b><br>", "Meta Surfaceome: <b>", meta_surfaceome, "</b><br>", "Replicate IDs: ", replicate_ids_left_right, "<br>" , "Feature IDs : ", features_left_right_overall, "<br>" ))
    ) +
      theme_few() +
      geom_point(fill = c("#cc66ff", "#cc0000", "#0000EE", "#5F5F61")[match(subset_inf_pos$suvo_plot_color, c("surface+string", "surface", "string", "other"))],
                 shape = ifelse(subset_inf_pos$imputed_comparison == "yes", 23, 21),
                 alpha = subset_inf_pos$plot_alpha,
                 size = 3,
                 color = "white",
                 stroke = 0.4) +
      geom_segment(x = -0.3, xend = 0.3, y = intensity_median, yend = intensity_median, col = "darkgrey", linetype = "solid", linewidth = 2) +
      geom_label_repel(box.padding = 2, size=6, max.overlaps = 100000) +
      ylim(cv_limits) +
      xlim(intensity_limits) +
      theme(legend.position = "none",
            text = element_text(size = 14),
            axis.text = element_text(size = 12)) +
      ylab("") +
      xlab("log2(Intensity)") +
      geom_segment(x=data_median_intensity_protein, y=-1, xend=data_median_intensity_protein, yend=100, linetype="dashed", col="darkgrey", linewidth = 0.4)  # add median signal line for orientation 
    
    
    # static super-volcano-plot  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    unique_arrange         <- ggarrange(plot_negative_unique_hits, plot_positive_unique_hits)
    composite              <- ggarrange(plot_volcano1,             unique_arrange, heights = c(1.5, 0.5), nrow= 2, align = "v")
    ggsave(file=paste("_", experiment_name, "_super-volcano-plot_",  loop_frame[comp_counter], ".png", sep=""), plot=composite, width=10, height=15, dpi = 300)
  } # end poi loop ----------------------------------------------------------------------------------------------------------------------------------------------------------
  setwd(output_directory)
  
  # interactive super-volcano-plot -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # handling missing plots -  those would lead to ggplotly() and subplot() crash ...
  if (nrow(plot_negative_unique_hits$data) == 0) {
    plot_negative_unique_hits <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "no left sided unique hits", size = 4) +
      theme_few()    }
  if (nrow(plot_positive_unique_hits$data) == 0) {
    plot_positive_unique_hits <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "no right sided unique hits", size = 4) +
      theme_few()    }     
  # ggplotly
  plot_interactive_volcano              <- ggplotly(plot_volcano1            , tooltip="text")       %>%  layout(hoverlabel = list(align = "left")) 
  # style(marker.symbol = c("circle", "diamond"))
  plot_interactive_negative_unique_hits <- ggplotly(plot_negative_unique_hits, tooltip="text")       %>%  layout(hoverlabel = list(align = "left")) 
  plot_interactive_positive_unique_hits <- ggplotly(plot_positive_unique_hits, tooltip="text")       %>%  layout(hoverlabel = list(align = "left")) 
  # Combine interactive plots into one (nice here is upper volcano is larger. Issue axis are lost ...)
  unique_arrange_interactive        <- subplot(plot_interactive_negative_unique_hits, plot_interactive_positive_unique_hits, nrows = 1, widths  = c(0.5,0.5)  , titleX = TRUE, titleY = TRUE) 
  composite_interactive             <- subplot(plot_interactive_volcano             , unique_arrange_interactive           , nrows = 2, heights = c(2/3,1/3)  , titleX = TRUE, titleY = TRUE) %>%
    layout(showlegend = FALSE, title = 'Super Volcano')
  #composite_interactive
  saveWidget(composite_interactive, file=paste0("_", experiment_name, "_super-volcano-plot_interactive_", loop_frame[comp_counter], ".html") ) # save interactive plot as HTML
  
  # export testing result in toomuchvolcano compatible format (colnames adjusted)
  write.csv(super_volcano_data[order(super_volcano_data$left_right, decreasing = TRUE),] %>% rename("Protein" = "entry", "Label" = "comparison"),
            file=paste(experiment_name,"__superVolcanoData__",loop_frame[comp_counter],".csv", sep=""), row.names=F)
  
  
  ## super_volcano_candi_prot_L - append data with each loop  =====================================================================================================================================================================================
  # *redundant for protti wf where ref condition is specified ... *   Aggregate super_volcano_data Data over loops 
  if(usedConditions[1] == paste(unique(result_subset$ttest_reference))) {
    super_volcano_condition <- usedConditions[2]     # set what is left of "_vs_" as condition of interest
  } else {
    super_volcano_condition <- usedConditions[1]    # set what is right of "_vs_" as condition of interest
  }
  
  super_volcano_data$condition  <- paste(super_volcano_condition)
  super_volcano_candi_prot_L <- rbind(super_volcano_candi_prot_L, 
                                      super_volcano_data[super_volcano_data$sig_right == "yes" | super_volcano_data$sig_left == "yes" | super_volcano_data$aggreg_unique_hits %in% c("+","-"),
                                                         c("condition", "comparison", "entry", "entry_name", "log2FC", "aggreg_unique_hits", "adj_pvalue", "pvalue", #"q_value",
                                                           "features_left", "features_right",  "features_left_right_overall", "replicate_ids_left", "replicate_ids_right",  "replicate_ids_left_right",
                                                           "cspa_2015", "cd_antigen", "surfy_2018", "meta_surfaceome", "poi", "condition_median_imp_log2_intensity", "protein_name", "imputed_left_right", "string_interactor")]
  )
  
  ## median protein intensity XY scatter plot  ==================================================================================================================================================================================================
  # filter for Groups of interest + subsets df + calculate no_df_right / Protein-Condtion combo + add imputation info
  scatter_plot_data <- data_level_prot[data_level_prot$comparison == paste0(usedConditions[1], "_vs_", usedConditions[2]), c("condition", "entry", "entry_name", "protein_name", "condition_median_imp_log2_intensity", 
                                                                                                                             "replicate_ids_left", "replicate_ids_right", "replicate_ids_left_right", "features_left_right_overall", "imputed_comparison")] %>% # retrieve protein level data for current comparison
    distinct(entry, condition, .keep_all = TRUE) %>% # kick values of condition replicates
    # surface and poi annotation 
    mutate(cd_antigen = "no", uniprot_2023 = "no",  tcsa_2021 = "no", surfy_2018 = "no", cspa_2015 = "no", cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high = "no")  %>% # create empty columns
    mutate(across(c(cd_antigen, uniprot_2023, tcsa_2021, surfy_2018, cspa_2015, cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high), ~ case_when(
      entry %in% surface_annotations[[cur_column()]] ~ "yes",
      grepl("cont", entry, ignore.case = TRUE) ~ "cont",
      TRUE ~ .x
    ))) %>%
    rename(meta_surfaceome = cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high)  %>%
    mutate(poi =  ifelse(entry_name %in% poi_var, "yes", "no"))   %>%
    # transform to plotting format
    pivot_wider(id_cols = c("entry", "entry_name", "protein_name", "replicate_ids_left_right", "features_left_right_overall", "imputed_comparison", 
                            "cspa_2015",  "surfy_2018", "tcsa_2021", "uniprot_2023", "meta_surfaceome", "cd_antigen", "poi"), 
                names_from = "condition", 
                values_from = "condition_median_imp_log2_intensity") %>%
    filter(if_any(all_of(usedConditions), ~!is.na(.)))  %>% # kick full NA columns (data_level_prot contains all combinations of entry and condition accross experiment - so double sided NAs)
    arrange(meta_surfaceome) # surface proteins on top
  
  ## annotate string (& based on that color) and plot alpha
  # string_interactors info introduced;    fyi string_interactors defined in result_subset prior volcano plot section dependent on target (= loop iteration)
  scatter_plot_data  <-  scatter_plot_data %>%
    mutate(string_interactor = case_when(entry_name %in% string_interactors[[1]] ~ "yes", TRUE ~ "no") ) %>%    # introduce string
    # plot dot colors
    mutate(plot_color = case_when( # specify colors
      scatter_plot_data[[plot_fill_volcano]] == "yes" & string_interactor == "yes" ~ "surface+string",    # is both plot_fill_volcano and string_interactor
      scatter_plot_data[[plot_fill_volcano]] == "yes" ~ "surface",                                        # is plot_fill_volcano
      string_interactor == "yes" ~ "string",                                                               # is string_interactor
      TRUE ~ "other"                                                                                       # is neither 
    )) %>%
    # plot alpha 
    mutate(plot_alpha = case_when(replicate_ids_left_right %in% c("1_vs_1", "1_vs_2", "2_vs_1", "0_vs_1", "1_vs_0", "2_vs_0", "0_vs_2") ~ plot_no_df_occupacy,  TRUE ~ 1.0 )) 
  
  
  xyMax <-  max( c(scatter_plot_data[[usedConditions[1]]], scatter_plot_data[[usedConditions[2]]]) ,  na.rm = TRUE) 
  xyMin <-  min( c(scatter_plot_data[[usedConditions[1]]], scatter_plot_data[[usedConditions[2]]]) ,  na.rm = TRUE) 
  # replace missing values by min for nice plot range
  scatter_plot_data <- scatter_plot_data %>%   
    mutate(
      !!sym(usedConditions[1]) := case_when(
        !is.na(!!sym(usedConditions[2])) & is.na(!!sym(usedConditions[1])) ~ xyMin,     #  2 is not NA but 1 is NA --> set 1 to min
        TRUE ~ !!sym(usedConditions[1])      ),
      !!sym(usedConditions[2]) := case_when(
        !is.na(!!sym(usedConditions[1])) & is.na(!!sym(usedConditions[2])) ~ xyMin,     #  1 is not NA but 2 is NA --> set 2 to min 
        TRUE ~ !!sym(usedConditions[2])      ),
      filter_col = case_when(
        is.na(!!sym(usedConditions[1])) & is.na(!!sym(usedConditions[2])) ~ NA_real_,   #  both are NA --> filter out 
        TRUE ~ 1      )) %>%
    filter(!is.na(filter_col)) %>%
    dplyr::select(-filter_col)
  
  # condition scatter plot ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  plot_condition_scatter <- ggplot(scatter_plot_data, 
                                   aes(x = .data[[usedConditions[1]]], 
                                       y = .data[[usedConditions[2]]], 
                                       fill = plot_color, 
                                       shape = imputed_comparison, 
                                       text = paste0("Protein:                   <b>", entry_name,          "</b><br>",
                                                     "Protein Name:        <b>", protein_name,              "</b><br>",
                                                     "Meta Surfaceome:  <b>", meta_surfaceome,        "</b><br>",
                                                     "Replicate IDs:  ", replicate_ids_left_right,    "<br>",
                                                     "Feature IDs :   ", features_left_right_overall, "<br>")       )) +
    theme_few() +
    geom_point(aes(alpha = plot_alpha),
               size = 2, 
               color = "white",
               stroke = 0.2) +
    scale_shape_manual(values = c("no" = 21, "yes" = 23))+ #c("left" = 25, "right" = 24, "no" = 21, "both" = 23)) +  # Ensure all points use a fillable shape
    scale_fill_manual(values = c("surface+string" = "#cc66ff", 
                                 "surface"        = "#cc0000", 
                                 "string"         = "#0000EE", 
                                 "other"          = "#5F5F61"),
                      guide = guide_legend(override.aes = list(shape = 21))) +  # Force legend to show filled shapes
    scale_alpha_identity() +  # Use the alpha values directly from the data
    guides(shape = "none") + # Remove the shape legend if not needed
    scale_x_continuous(limits = c(xyMin-0.1, xyMax))+
    scale_y_continuous(limits = c(xyMin-0.1, xyMax))+
    geom_abline(intercept = 0                    , slope = 1, color = "black", linetype = "solid" , size = 0.3) +
    geom_abline(intercept =  filter_log2fc_cutoff, slope = 1, color = "black", linetype = "dashed", size = 0.3) +  
    geom_abline(intercept = -filter_log2fc_cutoff, slope = 1, color = "black", linetype = "dashed", size = 0.3) 
  
  # interactive plottin, saving, clean up scatter variables
  plot_interactive_condition_scatter <- ggplotly(plot_condition_scatter, tooltip= "text") %>%  # Interactive  
    layout(hoverlabel = list(align = "left"))
  saveWidget(plot_interactive_condition_scatter, file=paste(experiment_name, "_ConditionScatter_Interactive_", loop_frame[comp_counter], ".html", sep="")) # save interactive plot as HTML
  rm(composite_interactive, scatter_plot_data, xyMax, xyMin, plot_condition_scatter, plot_interactive_condition_scatter) # clean environment
  
  # clean workspace --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  print(paste("Completed evaluation of", loop_frame[comp_counter], sep = " "))
  rm(list = c("quant_data_subset", "result_subset", #"super_volcano_data",     # data
              "subset_volcano", "subset_inf" , "subset_no_df",             # Plot subsets
              "prot_level_data_uni", "unique_arrange",               # Plot subsets
              "plot_volcano1", "plot_NaN", "plot_negative_unique_hits", "plot_positive_unique_hits", "plot_combined_volcano1", "composite", # Plots
              "plot_interactive_volcano", "plot_interactive_NaN",                 # Interactive plots
              "plot_interactive_negative_unique_hits", "plot_interactive_positive_unique_hits",   # Interactive plots
              "plot_combined_interactive_volcano",
              "n_features_per_condition", "features", "n_features_per_condition.exp", "measurements", "n_measurments_per_condition", 
              "super_volcano_condition"    ))
  
} # condition loop end 
write.csv(super_volcano_candi_prot_L, file= "_super_volcano_candi_prot_L.csv") # write result files
xxx <- super_volcano_candi_prot_L       # backup variable for super_volcano_candi_prot_L df
#save.image("Workspace_Image.RData")     # save workspace for next time (can be loaded from global environment)
rm(list = c(paste0("qc_plot_protti_", c(1:20)), 
            paste0("qc_plot_surface_", c(1:9))      ))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#  gLUX specifc analysis 
######### barplots for surface scatters all relative to same control (quick and dirty analysis) ____________________________________________________________________________________
# plotting unique hits
# plot_data <- data_level_prot %>%
#   filter(replicate_ids_left_right == "0_vs_1" | replicate_ids_left_right == "0_vs_2"  | replicate_ids_left_right == "0_vs_3") %>%     # commented out for modification filter figure
#   dplyr::select(condition, entry, meta_surfaceome) %>%
#   group_by(condition, meta_surfaceome) %>%
#   distinct() %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   mutate(meta_surfaceome = factor(meta_surfaceome, levels = c("yes", "no")))
# # barplot
# ggplot(plot_data, aes(x = condition, y = count, fill = meta_surfaceome)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_fill_manual(values = c("yes" = "firebrick", "no" = "black")) +
#   labs(title = "Unique Hits by Condition",
#        x = "Condition",
#        y = "Count",
#        fill = "meta_surfaceome") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Volcano enrichment summary barplot sigupunique (p-value filtered)
plot_data <- super_volcano_candi_prot_L %>%
  filter(log2FC == Inf) %>%
  dplyr::select(comparison, entry, meta_surfaceome) %>%
  group_by(comparison, meta_surfaceome) %>%
  distinct() %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(meta_surfaceome = factor(meta_surfaceome, levels = c("yes", "no")))
# barplot
ggplot(plot_data, aes(x = comparison, y = count, fill = meta_surfaceome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("yes" = "firebrick", "no" = "black")) +
  labs(title = "Unique Hit Plot (~volcano plot bottom right)",
       x = "Condition",
       y = "Count",
       fill = "meta_surfaceome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Scatterplot enrichment summary barplot (up+unique)
# sampe for diff abundance data to also include significant fold change proteins in barplot
condition_pairs <- unique(super_volcano_candi_prot_L$comparison)
plot_data <- data_prot_diff_abundance %>%
  filter(comparison %in% condition_pairs) %>% # only use desired comparisons (for t-test regular and "all" comparison being generated)
  filter(log2FC == Inf | log2FC > filter_log2fc_cutoff) %>%
  mutate(meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no")) %>% # for visualization
  dplyr::select(comparison, entry, meta_surfaceome) %>%
  group_by(comparison, meta_surfaceome) %>%
  distinct() %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(meta_surfaceome = factor(meta_surfaceome, levels = c("yes", "no")))
# barplot
ggplot(plot_data, aes(x = comparison, y = count, fill = meta_surfaceome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("yes" = "firebrick", "no" = "black")) +
  labs(title = "FC filtered protein plot (~scatter plot)",
       x = "Condition",
       y = "Count",
       fill = "meta_surfaceome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Volcano enrichment summary barplot sigupunique (p-value filtered)
plot_data <- super_volcano_candi_prot_L %>%
  filter(log2FC == Inf |
           (log2FC > filter_log2fc_cutoff & !is.na(pvalue) & pvalue <= filter_sig_cutoff)) %>%
  dplyr::select(comparison, entry, meta_surfaceome) %>%
  group_by(comparison, meta_surfaceome) %>%
  distinct() %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(meta_surfaceome = factor(meta_surfaceome, levels = c("yes", "no")))
# barplot
ggplot(plot_data, aes(x = comparison, y = count, fill = meta_surfaceome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("yes" = "firebrick", "no" = "black")) +
  labs(title = "FC & p-value filtered protein plot (~volcano plot)",
       x = "Condition",
       y = "Count",
       fill = "meta_surfaceome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################
##  below part contains +/- sophisticated random ideas for data visualizations  ###################################################################################################################################################################
###################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################

# retrieve POI traces accross conditions for clustering and trace plotting ########################################################################################################################################################################
prot_level_data_protein_traces <-  data_level_prot[, c("entry", "entry_name", "condition", "condition_median_imp_log2_intensity")] %>%
  group_by(entry, condition) %>%
  mutate(poi = ifelse(entry_name %in% poi_var, "yes", "no")) %>%
  filter(poi == "yes")%>%
  ungroup() %>%
  distinct(condition, entry, entry_name, condition_median_imp_log2_intensity, poi) %>%
  pivot_wider(id_cols = c(entry, entry_name, poi), names_from = condition, values_from = condition_median_imp_log2_intensity)  

# clustering of protein abundance trends -----------------------------------------------------------------------------------------------------------------------------------------------------------
library(cluster)
prot_level_data_protein_traces[is.na(prot_level_data_protein_traces)] <- 0   
data_scaled                 <- scale(prot_level_data_protein_traces[, -c(1:3)])   # Scale the data (optional, but recommended for k-means clustering)
cluster_proteins <- function(k) {                       # Function to perform k-means clustering
  set.seed(123)                                         # For reproducibility
  clusters <- kmeans(data_scaled, centers = k, nstart = 25)
  return(clusters)
}
clustering <- cluster_proteins(k = 4)        # specify number of clusters here
prot_level_data_clustered <- prot_level_data_protein_traces %>%
  mutate(cluster = clustering$cluster)# transfer clustering info into wide df

prot_level_data_clustered <- prot_level_data_clustered %>%
  pivot_longer(cols = setdiff(colnames(prot_level_data_protein_traces), c("entry","entry_name", "poi", "cluster")),  # this will fail as soon as other col than poi selected
               names_to = "condition",
               values_to = "condition_median_imp_log2_intensity")

plot_colors <- c("#FF6347", "#00CED1", "#9370DB", "#FFD700", "#8B008B", "#00FF7F",
                 "#FF8C00", "#9400D3", "#00FFFF", "#FF1493", "#00BFFF", "#8A2BE2",
                 "#FFA07A", "#00FA9A", "#EE82EE", "#FFD39B", "#006400", "#BDB76B",
                 "#8B0000", "#556B2F", "#FF69B4", "#CD5C5C", "#6B8E23", "#9932CC",
                 "#8FBC8F", "#E9967A", "#483D8B", "#2F4F4F", "#800000", "#000080")
# library(viridis) # when plotting hundrets of proteins
# plot_colors <- viridis(1000)

# cluster wise plotting of proteins -----------------------------------------------------------------------------------------------------------------------------------------------------------------
for (counter in 1:max(prot_level_data_clustered$cluster)) { 
  plot_protein_traces <- ggplot(prot_level_data_clustered[prot_level_data_clustered$cluster == counter, ], aes(x = factor(condition, levels = setdiff(colnames(prot_level_data_protein_traces), c("entry","entry_name", "poi", "cluster"))), # c("NAct", "Act_0", "Act_0_5h", "Act_1h", "Act_2h", "Act_4h", "Act_24h", "Act_72h"))
                                                                                                               y = condition_median_imp_log2_intensity, color = entry_name, group = entry)) +
    geom_line(size = 1) +
    scale_color_manual(values = plot_colors)+
    theme_minimal()+
    ylim(min(prot_level_data_clustered[prot_level_data_clustered$condition_median_imp_log2_intensity >0, ]$condition_median_imp_log2_intensity)-0.1, max(prot_level_data_clustered$condition_median_imp_log2_intensity)+0.1)+
    labs(title = "Intensity along Activation",
         x = "condition",
         y = "medianLog2Int",
         color = "entry")+
    geom_point(shape = 18, size = 4)+
    geom_text(aes(label = entry_name), nudge_y = 0.1, check_overlap = TRUE)
  
  print(plot_protein_traces)
}
paste(length(unique(prot_level_data_clustered$entry)), " POIs proteins were identified", sep="")




## surfaceome category plot for candiproteins (sigup+unique hits)-------------------------------------------
# Define the categories and their subcategories in the desired order
categories <- list(
  "T-cell"    = c("tcr_signaling_plus", "tcr_signaling_minus", "t_act_plus", "t_act_minus"),
  "Receptor"  = c("receptor_complex", "receptor_clustering", "receptor_transactivation", "receptor_internalization", "receptor_recycling", "receptor_metabolism", "mhc_complex"),
  "Signaling" = c("signal_transduction", "kinase", "phosphatase", "gpcr", "antigen_binding", "growth_factor_binding", "cytokine_receptor", "immune_response"),
  "Adhesion"  = c("cell_adhesion", "integrin_binding"),
  "Transport" = c("transporter_activity", "ion_channel"),
  "Others"    = c("hydrolase_activity", "redox_homeostasis", "oxidoreductase_activity", "apoptotic_process")
)
# Flatten the list to get the order of subcategories
subcategory_order <- unlist(categories)
# Create a mapping of subcategory order to category
subcategory_to_category <- data.frame(
  subcategory_order = seq_along(subcategory_order),
  subcategory = subcategory_order,
  category = rep(names(categories), sapply(categories, length))
)
subcategory_to_category$key <- rownames(subcategory_to_category)
# Create a mapping dataframe
category_mapping <- data.frame(
  subcategory = subcategory_order,
  category = rep(names(categories), sapply(categories, length))
)

# retrieve sigup and unique data subset & add uniprot derived go-term info
svcp_L_sigupuni <- super_volcano_candi_prot_L %>% 
  filter(log2FC > filter_log2fc_cutoff, # filter FC (Inf means pvalue is NA !!! thats why below ....
         (is.na(!!sym(pvalue)) & log2FC == Inf) | !!sym(pvalue) < filter_sig_cutoff)   %>%   # ... it is important to account for NA cases for Inf FC to not kick these IDs
  dplyr::select("comparison", "entry", "entry_name", "meta_surfaceome", "replicate_ids_left_right") %>%
  left_join(proteome %>% dplyr::select("entry", "gene_ontology_i_ds"), #  add go term info
            by = "entry") %>%
  mutate(   # annotate surfaceome categories of interest (extend list if needed and communicate suggestion to Martin)
    tcr_signaling_plus       = as.integer(str_detect(gene_ontology_i_ds, "GO:0050852")),
    tcr_signaling_minus      = as.integer(str_detect(gene_ontology_i_ds, "GO:0050860")),
    t_act_plus               = as.integer(str_detect(gene_ontology_i_ds, "GO:0042110")),
    t_act_minus              = as.integer(str_detect(gene_ontology_i_ds, "GO:0050868")),
    
    receptor_complex         = as.integer(str_detect(gene_ontology_i_ds, "GO:0043235")),
    receptor_clustering      = as.integer(str_detect(gene_ontology_i_ds, "GO:0043113")),
    receptor_transactivation = as.integer(str_detect(gene_ontology_i_ds, "GO:0035624")),
    receptor_internalization = as.integer(str_detect(gene_ontology_i_ds, "GO:0031623")),
    receptor_recycling       = as.integer(str_detect(gene_ontology_i_ds, "GO:0001881")),
    receptor_metabolism      = as.integer(str_detect(gene_ontology_i_ds, "GO:0043112")),
    
    signal_transduction      = as.integer(str_detect(gene_ontology_i_ds, "GO:0007165")),
    kinase                   = as.integer(str_detect(gene_ontology_i_ds, "GO:0004672")),
    phosphatase              = as.integer(str_detect(gene_ontology_i_ds, "GO:0016791")),
    gpcr                     = as.integer(str_detect(gene_ontology_i_ds, "GO:0004930")),
    antigen_binding          = as.integer(str_detect(gene_ontology_i_ds, "GO:0003823")),
    mhc_complex              = as.integer(str_detect(gene_ontology_i_ds, "GO:0042611")),
    growth_factor_binding    = as.integer(str_detect(gene_ontology_i_ds, "GO:0019838")),
    cytokine_receptor        = as.integer(str_detect(gene_ontology_i_ds, "GO:0004896")),
    immune_response          = as.integer(str_detect(gene_ontology_i_ds, "GO:0006955")),
    
    cell_adhesion            = as.integer(str_detect(gene_ontology_i_ds, "GO:0007155")),
    integrin_binding         = as.integer(str_detect(gene_ontology_i_ds, "GO:0005178")),
    
    transporter_activity     = as.integer(str_detect(gene_ontology_i_ds, "GO:0005215")),
    ion_channel              = as.integer(str_detect(gene_ontology_i_ds, "GO:0005216")),
    hydrolase_activity       = as.integer(str_detect(gene_ontology_i_ds, "GO:0016787")),
    oxidoreductase_activity  = as.integer(str_detect(gene_ontology_i_ds, "GO:0016491")),
    apoptotic_process        = as.integer(str_detect(gene_ontology_i_ds, "GO:0006915")),
    redox_homeostasis        = as.integer(str_detect(gene_ontology_i_ds, "GO:0045454"))# immune_response         = ifelse(str_detect(gene_ontology_i_ds, "GO:0006955"), 1, 0),
  )

sums_by_condition <- svcp_L_sigupuni %>%    # category summary per condition & join with annotation
  group_by(comparison) %>%
  summarise(across(all_of(subcategory_order), sum)) %>%   # sum across columns defined in variable
  pivot_longer(cols = -comparison, names_to = "subcategory", values_to = "sum") %>%
  mutate(category_mapping = row_number()) %>%
  left_join(subcategory_to_category, by = c("subcategory" = "key")) %>% # join with annotation 
  dplyr::select(-subcategory_order, -subcategory) %>%
  rename("subcategory" = "subcategory.y") %>%
  mutate(
    category    = factor(category, levels = names(categories)),
    subcategory = factor(subcategory, levels = subcategory_order)
  ) 

# category bar plotting
colors <- c("#377eb8", "#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#6a3d9a")

svcp_category_plot <- ggplot(sums_by_condition, aes(x = subcategory, y = sum, fill = comparison)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(title = "sum of categories by condition",
       y = "sum",
       fill = "comparison") +
  facet_grid(. ~ category, scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(angle = 90, hjust = 0))

ggsave(file=paste0("_", experiment_name, "_svcp_category_plot_", ".png"), plot=svcp_category_plot, width = 20, height = 10)
svcp_category_plot

## upset plot & category heatmap ########################################################################################################################################################################     
# groups in 
svcp_L_sigupuni <- svcp_L_sigupuni %>%
  arrange(condition) %>%
  group_by(entry) %>%
  mutate(conditions_listed = paste(unique(condition), collapse = "_")) %>%
  ungroup()

# create binary matrix in wide format
upset_data <- svcp_L_sigupuni[,c("condition","entry", "entry_name")] %>%
  mutate(value = 1) %>%
  mutate(entry_name = stringr::str_replace_all(entry_name, "_HUMAN", "")) %>%
  pivot_wider(names_from = condition, values_from = value, values_fill = 0)

# create UpSet plot
upset(as.data.frame(upset_data), 
      sets = unique(svcp_L_sigupuni$condition),
      order.by = "freq", 
      nsets = length(conditions),  # Show all conditions
      point.size = 3, 
      line.size = 1,
      mainbar.y.label = "Number of Entries",
      sets.x.label = "Entries per Condition",
      text.scale = c(2, 2, 1.6, 1.6, 2, 1.6))  # Adjust text size
# aggregate info for categories per overlap-subset and calculate % of subset made up by category  xyz
heatmap_data <- svcp_L_sigupuni %>%
  group_by(conditions_listed) %>%
  summarize(across(c(immune_response, cytokine_receptor, mhc_protein, cell_adhesion, integrin_binding, 
                     signal_transduction, protein_kinase, gpcr_activity, transporter_activity, ion_channel, 
                     hydrolase_activity, oxidoreductase_activity, antigen_binding, t_cell_receptor, apoptotic_process, 
                     growth_factor_binding, redox_homeostasis), 
                   ~mean(. == 1, na.rm = TRUE) * 100),
            count = n()) %>%   # . ==1 yields TRUE / FALSE vector dependent on col value. when forming mean of logical vector the fraction of total list which is TRUE is the result; *100 to obtain %
  arrange(desc(count)) %>%
  dplyr::select(-count) %>%
  pivot_longer(cols = -conditions_listed, 
               names_to = "category", 
               values_to = "value") %>%
  mutate(
    conditions_listed = factor(conditions_listed, levels = names(sort(table(svcp_L_sigupuni$conditions_listed), decreasing = TRUE))),
    category          = factor(category         , levels = rev(subcategory_order))  # Reverse for bottom-to-top ordering
  )

# Create the heatmap
ggplot(heatmap_data, aes(x = conditions_listed, y = category, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, max(heatmap_data$value))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Enrichment Overlap Group", y = "GO-Term Category", fill = "% of Group") +
  coord_fixed()

data <- svcp_L_sigupuni %>%
  dplyr::select(meta_surfaceome, conditions_listed) %>%
  table()

# convert data to a data frame
data <- as.data.frame(table(svcp_L_sigupuni$meta_surfaceome, svcp_L_sigupuni$conditions_listed))
colnames(data) <- c("meta_surfaceome", "conditions_listed", "count")

# Create the stacked bar plot
ggplot(data, aes(y = conditions_listed, x = count, fill = meta_surfaceome)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("no" = "grey", "yes" = "black")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(x = "Count", 
       y = "Enrichment Overlap Group", 
       fill = "Meta Surfaceome",
       title = "Distribution of Meta Surfaceome by Enrichment Overlap-/Subset") +
  coord_flip()  # This flips the coordinates to make it a horizontal bar plot


## sample wise heatmap of sig-up and unique IDs of each condition ########################################################################################################################################################################     
# create keys for data retrieval
data_level_prot <- data_level_prot %>%
  mutate(key = paste(condition, entry, sep = "_"))
svcp_L_sigupuni <- svcp_L_sigupuni %>%
  mutate(key = paste(condition, entry, sep = "_"))
# filter data_level_prot for sig up unique proteins in corresponding condition >> sample wise quant info for sig up unique 
prot_level_data_sigupuni <- data_level_prot %>% 
  filter(key %in% svcp_L_sigupuni$key) %>%
  dplyr::select(-key) %>%  # Remove the temporary key column
  dplyr::select(condition, bio_replicate, con_rep, entry, entry_name, protein_name, imputed_prot_intensity_log2, features_left_right_overall, replicate_ids_left_right)

library(tidyverse)
library(pheatmap)
# Prepare the data for the heatmap
heatmap_data <- prot_level_data_sigupuni %>%
  dplyr::select(con_rep, entry_name, imputed_prot_intensity_log2) %>%
  pivot_wider(names_from = con_rep, values_from = imputed_prot_intensity_log2) %>%
  mutate(across(-entry_name, ~replace_na(., 0))) %>%  # Impute NA as 0
  column_to_rownames("entry_name")

# Create annotation dataframe for entry_name
annotation_row <- data.frame(
  entry_name = rownames(heatmap_data),
  row.names  = rownames(heatmap_data)
)
# Create annotation dataframe for con_rep
annotation_col <- data.frame(
  condition = sub("_.*", "", colnames(heatmap_data)),
  replicate = sub(".*_", "", colnames(heatmap_data)),
  row.names = colnames(heatmap_data)
)
# Create a custom color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# Find the midpoint (corresponding to 0)
mid <- max(abs(min(heatmap_data, na.rm = TRUE)), abs(max(heatmap_data, na.rm = TRUE)))
breaks <- seq(-mid, mid, length.out = 101)
# Create the heatmap
plot_heatmap_prot_sigupuni <- pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  #  annotation_row = annotation_row,
  annotation_col = annotation_col,
  main = "Heatmap of Protein Intensity (Significantly Up-regulated and Unique)",
  fontsize_row = 8,
  fontsize_col = 8,
  scale = "row",  # Scale rows for better visualization
  color = color_palette,
  breaks = seq(-2, 2, length.out = 101)
)
# Save the heatmap as a PDF
# pdf("protein_intensity_heatmap_sigupuni.pdf", width = 12, height = 10)
# print(p)
# dev.off()


## condition wise heatmap of sig-up and unique IDs of each condition ########################################################################################################################################################################     
# Create the heatmap
pheatmap(upset_data %>% 
           dplyr::select(-entry) %>%
           column_to_rownames("entry_name"), 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Condition wise Clustered Heatmap of SigUp & Unique Hits",
         color = colorRampPalette(c("grey", "white", "black"))(50),
         fontsize_row = 6,
         fontsize_col = 10)




?order
## To do ###################################################################################################################################################################################################################################################
# Labelling in Volcano is still a bit off, need to optimize that
# add functionality that the plotting criteria can be changed (ATM just for Glyco Motif)
# Peptide based quantification and plotting -> find out how to do that in MSStats
