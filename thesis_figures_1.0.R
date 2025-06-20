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
library(UpSetR)
library(ggupset)

# *****************************************************************************************************************************************************************
# general plotting theme from Amanda ******************************************************************************************************************************
# my plot theme to fit most journal figure requirements:
# black axis ticks and borders
# 5-7pt or 6-8pt text for a 2-panel figure
# Line width 0.5-1.5 pt
# Colour blind friendly colour scheme (would recommend defining colours in source functions and use throughout) e.g. no red and green used together
# Gray fills between 10-80%
# Figure resolution at least 300 dpi ( I like to export in pdf not svg or tiff if using ggsave so figure is small and vectorized)
plot_theme <- function() {
  theme_classic() %+replace%
    theme(
      axis.text     = element_text(size = 16, face = "plain", family = "Arial"),
      axis.text.x   = element_text(family = "Arial"),
      axis.title    = element_text(size = 16, face = "plain", family = "Arial"),
      axis.title.x  = element_text(family = "Arial"),
      title         = element_text(size = 16, face = "plain", family = "Arial"),
      strip.text.x  = element_text(size = 16, family = "Arial"),
      strip.text.y  = element_text(size = 16, family = "Arial"),
      legend.text   = element_text(size = 16, family = "Arial"),
      legend.title  = element_text(size = 16, face = "plain", family = "Arial"),
      panel.border  = element_rect(colour = "black", fill = NA, size = 1),
      plot.margin   = unit(c(t = 1, r = 1, b = 0.1, l = 1), "cm"),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
}
# plot_theme <- function() { # Martin
# theme_classic() %+replace%
#   theme(
#     axis.text = element_text(size = 12, face = "plain", family = "Arial"),
#     axis.text.x = element_text(family = "Arial"),
#     axis.title = element_text(size = 12, face = "bold", family = "Arial"),
#     axis.title.x = element_text(family = "Arial"),
#     title = element_text(size = 12, face = "plain", family = "Arial"),
#     strip.text.x = element_text(size = 12, family = "Arial"),
#     strip.text.y = element_text(size = 12, family = "Arial"),
#     legend.text = element_text(size = 12, family = "Arial"),
#     legend.title = element_text(size = 12, face = "bold", family = "Arial"),
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     plot.margin = unit(c(t = 1, r = 1, b = 0.1, l = 1), "cm"),
#     plot.background = element_rect(fill = "transparent", color = NA)
#   )
# } 
# ***************************************************************************************************************************************************************
# *****************************************************************************************************************************************************************

##################################################################################################################################################################################################################
## GENERAL RESOURCES ####################################################################################################################################################################################################
##################################################################################################################################################################################################################
CSPA <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_per_cell_type.csv") %>%
  filter(protein_count > 0) %>%
  rename(entry = id_link) %>%
  pull(entry) %>%
  unique()
paste("CSPA comprises", length(CSPA), "proteins")
CSPA_Jurkat <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_per_cell_type.csv") %>%
  filter(jurkat > 0) %>%
  rename(entry = id_link) %>%
  pull(entry) %>%
  unique()
paste("Jurkat CSPA comprises", length(CSPA_Jurkat), "proteins")
#
surface_annotations = read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/surface.annotations.csv", header = TRUE, sep = ",")
# 
meta_surfaceome   <- surface_annotations %>% 
  filter(cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high != "") %>%
  pull(cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high)
paste("Meta surfaceome comprises", length(meta_surfaceome), "proteins")
#
proteome_upsp_202501 <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/human_upsp_202501.csv")
paste("Human proteome (upsp 2025-01) comprises", length(unique(proteome_upsp_202501$entry)), "proteins")
#

##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################




##################################################################################################################################################################################################################
# Chapter 2 ######################################################################################################################################################################################################
# Figure 2.3.1: CSC LLOQ panT  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
v25_LLOQ_CSC <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_CSC_panT_lowinput/_output_1-2_ludo_adjp_0.7string_shs2.25/_data_prot_level.csv") %>%
  mutate(condition = gsub("0_5", "0.5", condition)) 

v25_LLOQ_CSC_pep <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_CSC_panT_lowinput/_output_1-2_ludo_adjp_0.7string_shs2.25/data_raw_glyco_nZ_gyIFcsc.csv") %>%
  mutate(condition = gsub("0_5", "0.5", condition))   

# Fig2.3.1_b: CSC Peptides ......................................................................................................................
Fig2.3.1_b <- v25_LLOQ_CSC_pep %>%
  dplyr::select(condition, bio_replicate, peptide_sequence_mod, csc_signature) %>%
  distinct() %>%
  group_by(condition, csc_signature, bio_replicate) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(condition, csc_signature) %>%
  summarise(
    mean_count = mean(count),
    sd = sd(count),
    se = sd / sqrt(n()),
    median_count = median(count)  # New line added
  ) %>%
  #
  ggplot(aes(x = condition, y = mean_count, fill = csc_signature)) +
  geom_col(position = position_dodge(0.9), width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_count - se, ymax = mean_count + se),
    position = position_dodge(0.9),
    width = 0.25,
    color = "black"  ) +
  scale_fill_manual(values = c("#e15759", "darkgrey", "#4e79a7")) +
  labs(
    title = "Peptides",
    x     = "panT cell input [e6]",
    y     = "Unique Peptides",
    fill  = "CSC signature"  ) +
  plot_theme()

Fig2.3.1_b
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_b.png",
  plot   = Fig2.3.1_b,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)


# Fig2.3.1_c: CSC Proteins ......................................................................................................................
Fig2.3.1_c  <- v25_LLOQ_CSC_pep %>%
  dplyr::select(condition, bio_replicate, entry, csc_signature) %>%
  distinct() %>%
  group_by(condition, csc_signature, bio_replicate) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(condition, csc_signature) %>%
  summarise(
    mean_count = mean(count),
    sd = sd(count),
    se = sd / sqrt(n()),
    median_count = median(count)  # New line added
  ) %>%
  #
  ggplot(aes(x = condition, y = mean_count, fill = csc_signature)) +
  geom_col(position = position_dodge(0.9), width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_count - se, ymax = mean_count + se),
    position = position_dodge(0.9),
    width = 0.25,
    color = "black"  ) +
  scale_fill_manual(values = c("#e15759", "darkgrey", "#4e79a7")) +
  labs(
    title = "Proteins",
    x     = "panT cell input [e6]",
    y     = "Unique Proteins",
    fill  = "CSC signature"  ) +
  plot_theme()

Fig2.3.1_c
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_c.png",
  plot   = Fig2.3.1_c,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)


# Fig2.3.1_d: CSC Precursor Signal ......................................................................................................................
Fig2.3.1_d <- v25_LLOQ_CSC_pep %>%
  filter(csc_signature != "cont") %>%
  group_by(con_rep, csc_signature) %>%
  mutate(TotInt = sum(raw_prec_intensity)) %>%
  select(condition, bio_replicate, con_rep, csc_signature, TotInt) %>%
  unique() %>%
  group_by(con_rep) %>%
  mutate(SampleInt = sum(TotInt)) %>%
  ungroup() %>%
  mutate(RelInt = TotInt / SampleInt) %>%
  filter(csc_signature == "yes") %>%
  group_by(condition) %>%
  summarise(
    mean_RelInt = mean(RelInt)*100,  # to 100%
    RelInt_sd = sd(RelInt)*100,      # to 100% 
    count = n(),
    RelInt_se = RelInt_sd / sqrt(count)
  ) %>%
  ggplot(aes(x = condition, y = mean_RelInt)) +
  geom_col(position = position_dodge(0.9), width = 0.8, fill = "#4e79a7") +
  geom_errorbar(
    aes(ymin = mean_RelInt - RelInt_se, ymax = mean_RelInt + RelInt_se),
    position = position_dodge(0.9),
    width = 0.25,
    color = "black"
  ) +
  labs(
    title = "CSC signal contribution",
    x     = "panT cell input [e6]",
    y     = "CSC feature signal [%]",
    fill  = "CSC signature"  ) +
  plot_theme()
  
Fig2.3.1_d
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_d.png",
  plot = Fig2.3.1_d,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_e: CSC Protein ID data completeness ......................................................................................................................
Fig2.3.1_e <- v25_LLOQ_CSC_pep %>%
  filter(csc_signature == "yes", raw_prec_intensity > 0) %>%
  select(condition, entry, bio_replicate) %>%
  distinct() %>%
  group_by(condition, entry) %>%
  mutate(detection_per_condition = n_distinct(bio_replicate)) %>%
  ungroup() %>%
  distinct(condition, entry, detection_per_condition) %>% # filter out replicates
  group_by(condition, detection_per_condition) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(data_completeness = round(detection_per_condition/3 * 100, 0)) %>%
  #
  ggplot(aes(x = condition, y = count, fill = factor(data_completeness))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Detection completeness",
       x     = "Experimental condition",
       y     = "Protein count",
       fill  = "Detection %") +
  scale_fill_manual(values = c("black", "#7F7F7F", "#CCCCCC")) +
  plot_theme() +
  theme(legend.position = "right")

Fig2.3.1_e
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_e.png",
  plot = Fig2.3.1_e,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  
  
# Fig2.3.1_X  Fig2.3.1_f glyco sites/peptides, glyco glyco-peptides/protein ......................................................................................................................
v25_LLOQ_CSC_pep_csc <- v25_LLOQ_CSC_pep %>%
  select(entry, entry_name, condition, csc_signature_psm, peptide_sequence_mod) %>%
  distinct() %>%
  filter(csc_signature_psm == "yes") %>%
  mutate(csc_signature_count = str_count(peptide_sequence_mod, fixed("N[0.9840]"))) %>%   # no stringent filter required because signature prefiltered in line above
  group_by(entry, condition) %>%
  mutate(csc_peptides_per_protein = n_distinct(peptide_sequence_mod)) %>%
  ungroup()

table(v25_LLOQ_CSC_pep_csc$csc_signature_count) # for most the numer is 1 (1746) only for 13,3 % it is 2 (232)
summary(v25_LLOQ_CSC_pep_csc$csc_peptides_per_protein)


Fig2.3.1_X <- ggplot(v25_LLOQ_CSC_pep_csc, aes(x = condition, y = csc_signature_count, fill = factor(csc_signature_count))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Signatures per peptide",
    x     = "panT cell input [e6]",
    y     = "Number of signatures",
    fill  = "Signature count") +
  scale_fill_manual(values = c("black", "#7F7F7F")) +
  plot_theme()  # Use built-in theme instead of plot_theme() until verified

Fig2.3.1_X
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_X.png",
  plot = Fig2.3.1_X,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  
#
#
Fig2.3.1_f <- ggplot(v25_LLOQ_CSC_pep_csc, aes(x = condition, y = csc_peptides_per_protein)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "#4e79a7") +
  labs(
    title = "CSC peptides / protein",
    x     = "panT cell input [e6]",
    y     = "Number of CSC peptides"
  ) +
  coord_cartesian(ylim = c(0, 30)) +  # Limits y-axis to 0–40
  plot_theme()  # Use built-in theme instead of plot_theme() until verified

Fig2.3.1_f
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_f.png",
  plot = Fig2.3.1_f,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_g overlap/upset Jurkat CSPA & pT CSC ......................................................................................................................
# write.csv(data.frame(Values = CSPA_Jurkat), "/Users/mgesell/Downloads/CSPA_Jurkat.csv", row.names = FALSE)
# Prepare data in required format
panT_csc_ids <- v25_LLOQ_CSC %>% 
  filter(csc_signature_psm == "yes") %>% 
  pull(entry) %>%
  unique()

upset_data <- bind_rows(
  tibble(entry = CSPA_Jurkat  , set = "Jurkat CSPA"),
  tibble(entry = panT_csc_ids , set = "panT cells") ,
  tibble(entry = CSPA         , set = "CSPA contained") 
  # tibble(entry = intersect(intersect(CSPA, CSPA_Jurkat), panT_csc_ids), # to not display full CSPA but show which proteins already annotate in CSPA  
  #                               set = "CSPA subset"),
  ) %>% 
  group_by(entry) %>% 
  summarise(sets = list(set)) %>%  # Critical: list column of set memberships
  ungroup() %>%
  filter(!(map_lgl(sets, ~ length(.) == 1 && all(. == "CSPA contained"))))

Fig2.3.1_g <- upset_data %>% 
  ggplot(aes(x = sets)) +
  geom_bar(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_upset(
    sets = c("Jurkat CSPA", "CSPA contained", "panT cells"),
    name = "",
    n_intersections = 20
  ) +
  labs(
    y     = "Intersection",
    title = "Protein Overlap"
  ) +
  plot_theme() +
  theme(axis.text.y = element_text(size = 14))

Fig2.3.1_g
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_g.png",
  plot = Fig2.3.1_g,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  
# fyi
cat("\n",
    "Novel on panT & CSPA supported ", length(setdiff(intersect(panT_csc_ids, CSPA), CSPA_Jurkat))  , "\n",
    "Overlap panT-CSPA              ", length(intersect(intersect(CSPA_Jurkat, panT_csc_ids), CSPA)), "\n",
    "Jurkat but not panT            ", length(setdiff(intersect(CSPA_Jurkat, CSPA), panT_csc_ids))  , "\n",
    "Novel panT                     ", length(setdiff(setdiff(panT_csc_ids, CSPA), CSPA_Jurkat))
)
# info to write in thesis
paste("Unique overall CSC PEPTIDES", v25_LLOQ_CSC_pep %>% filter(csc_signature == "yes") %>% pull(peptide_sequence_mod) %>% unique() %>% length() )
paste("Unique overall CSC PROTEINS", v25_LLOQ_CSC_pep %>% filter(csc_signature == "yes") %>% pull(entry) %>% unique() %>% length() )
paste("Consider reporting median CV values at least in text")

# =============================================================================================================================================================
# Figure 2.3.1_i: TCR-LUX LLOQ panT  ..........................................................................................................................................................................................................................................
v25_LLOQ_LUX_data_prot_diff_abundance <- read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_LUX_panT_lowinput_TCR/_output_1-2_ludo_adjp_0.7string_shs2.25/_data_prot_diff_abundance.csv") %>%
  arrange(comparison ) %>%
  group_by(entry_name) %>%
  mutate(overlap = paste(comparison, collapse = "_&_")) %>%
  ungroup()
# sig up data
v25_LLOQ_LUX_sigup <- v25_LLOQ_LUX_data_prot_diff_abundance %>%
  filter(log2FC >1, adj_pvalue <= 0.05) %>% # sig up 
  select(entry, entry_name, comparison, overlap) %>%
  mutate(comparison = str_replace_all(comparison,
                                      c("TCR_5e5_vs_Iso_5e5" = "0.5",
                                        "TCR_1e6_vs_Iso_1e6" = "1",
                                        "TCR_5e6_vs_Iso_5e6" = "5"
                                      )),
         comparison = factor(comparison, levels = c("0.5", "1", "5")),
         LUX_id = 1  )  %>%
  mutate(meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no"))

# core known TCR community
abTCR_chains <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv") %>% 
  filter(abTCR_chains_cd3_mgmanual != "") %>% 
  pull(abTCR_chains_cd3_mgmanual)
TCR_LUX_pois <- c(c("CD4_HUMAN", "CD8A_HUMAN" , "CD8B_HUMAN"), abTCR_chains)

# Figure 2.3.1_i: Volcano plots ..........................................................................................................................................................................................................................................
# consider to only display red grey blue purple annotation (surface / non-surface, string-interactors, string&surface )
filter_log2fc_cutoff            = log(2, 2)   # Cutoff for Volcano plotting
pvalue                = "adj_pvalue"          # select "pvalue" or "adj_pvalue" depending on how stringend you want to be.   recommendation for CSC is adj_pvalue    
filter_sig_cutoff               = 0.05  
#
plot_no_df_occupacy         =  0.7      # make lower confidence ids less bright
plot_label_volcano          =  "poi"               # default poi   
plot_label_inf              =  plot_label_volcano  # default plot_label_volcano
plot_fill_volcano           =  "meta_surfaceome"   # default meta_surfaceome  
plot_fill_label_inf         =  plot_fill_volcano   # default plot_fill_volcano                 
plot_fill_condition_scatter =  "meta_surfaceome" # x-y median condition signal plot (color surface or poi recommended)
#
count_var = 1
loop_frame <- unique(v25_LLOQ_LUX_data_prot_diff_abundance$comparison)
# comp_counter = 1
for (comp_counter in 1:length(loop_frame)) {   # ----------------------- LOOOOOOOP --------------------- LOOOOOOOP ------------------------------ LOOOOOOOP ------------------- LOOOOOOOP -------------------
  usedConditions <- c(gsub("_vs_.*", "", loop_frame[comp_counter]), gsub("^.*_vs_", "", loop_frame[comp_counter]))
  result_subset <- v25_LLOQ_LUX_data_prot_diff_abundance[v25_LLOQ_LUX_data_prot_diff_abundance$comparison == loop_frame[comp_counter],] 
  # adjust p value (<<-- after this point no more data filtering)
  result_subset$adj_pvalue = p.adjust(result_subset$pvalue, method="BH")
  result_subset$pvalue      [result_subset$pvalue        < 2.220446e-16 ] = 2.220446e-16
  result_subset$adj_pvalue  [result_subset$adj_pvalue    < 2.220446e-16 ] = 2.220446e-16
  result_subset <- result_subset[order(result_subset$adj_pvalue),]  # sort by significance - good on top
  result_subset <- result_subset[!is.na(result_subset$log2FC), ]   # filter out all proteins that were not identified - all other should have calculated FC or +/-Inf assigned
  #
  poi_var <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv" , header = TRUE, sep = ",") %>%
    filter(abTCR_chains_cd3_mgmanual != "") %>%
    pull(abTCR_chains_cd3_mgmanual)
  #
  result_subset <- result_subset %>%
    mutate(cspa_2015       = ifelse(entry %in% surface_annotations$cspa_2015      , "yes", "no"),
           surfy_2018      = ifelse(entry %in% surface_annotations$surfy_2018     , "yes", "no"),
           tcsa_2021       = ifelse(entry %in% surface_annotations$tcsa_2021      , "yes", "no"),
           meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no"),
           cd_antigen      = ifelse(entry %in% surface_annotations$cd_antigen     , "yes", "no"),
           uniprot_2023    = ifelse(entry %in% surface_annotations$uniprot_2023   , "yes", "no"),
           poi   =      ifelse(entry_name %in% poi_var                            , "yes", "no")    # pois are given as entry_names in MGs poi meta file
    )  

  # # string interactor mapping
  # string_interactors <- string[string$protein1_entry_name %in% string_targets[[ experiment_conditions[[comp_counter]] ]],] %>% 
  #   dplyr::select(protein2_entry_name) %>% 
  #   unique()
  # result_subset <- result_subset %>% # map interactors to data 
  #   mutate(string_interactor = case_when(entry_name %in% string_interactors[[1]] ~ "yes", TRUE ~ "no") ) 

  super_volcano_data <- result_subset
  # specify plot dot colors for suvo
  super_volcano_data  <-  super_volcano_data %>%
    mutate(suvo_plot_color = case_when(
      # super_volcano_data[[plot_fill_volcano]] == "yes" & string_interactor == "yes" ~ "surface+string",    # is both plot_fill_volcano and string_interactor
      super_volcano_data[[plot_fill_volcano]] == "yes" ~ "surface",                                        # is plot_fill_volcano
      # string_interactor == "yes" ~ "string",                                                               # is string_interactor
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
  super_volcano_data <- merge(super_volcano_data, proteome_upsp_202501[,c("entry", "function_cc", "gene_ontology_molecular_function", "glycosylation")], by = "entry")
  
  # order data. 1) meta_surfaceome 
  super_volcano_data <- super_volcano_data %>% 
    arrange(meta_surfaceome) # surface hits on bottom of df >> plotted last = on top
  
  ## For Plotting dataset is separated into different cases
  subset_volcano   <- super_volcano_data[!is.na(super_volcano_data[[pvalue]]), ]    # 1) Fold-chang

  foldchangelimit_max <- max(v25_LLOQ_LUX_data_prot_diff_abundance$log2FC)
  foldchangelimit_min <- min(v25_LLOQ_LUX_data_prot_diff_abundance$log2FC)
  
  
  subset_volcano <- subset_volcano %>%
    mutate(significance = as.factor(ifelse(-log10(.data[[pvalue]]) >= -log10(filter_sig_cutoff)   &    abs(log2FC) >= filter_log2fc_cutoff,   "yes",   "no")),  # Mark significant hits
           plot_label   = ifelse(.data[[plot_label_volcano]] == "yes", gsub("_HUMAN", "", entry_name), "")) # Assign plot label 
  
  # Helper frame for significance lines
  segmentation <- data.frame(x = c(foldchangelimit_min, filter_log2fc_cutoff, -filter_log2fc_cutoff, filter_log2fc_cutoff),
                             y = c(-log10(filter_sig_cutoff), -log10(filter_sig_cutoff), -log10(filter_sig_cutoff), -log10(filter_sig_cutoff)),
                             xend = c(-filter_log2fc_cutoff, foldchangelimit_max, -filter_log2fc_cutoff, filter_log2fc_cutoff),
                             yend = c(-log10(filter_sig_cutoff), -log10(filter_sig_cutoff), max(-log10(subset_volcano[[pvalue]])), max(-log10(subset_volcano[[pvalue]]))),
                             col = rep("black", times=4),
                             linetype = rep("dashed", times=4))
  
  # Volcano Plot 
  plot_volcano <- ggplot(subset_volcano,
                          aes(x=log2FC, y=-log10(.data[[pvalue]]),# label=plot_label, 
                              text= paste0("Protein: <b>", entry_name, "</b><br>", "Protein Name: <b>", protein_name, "</b><br>", "Meta Surfaceome: <b>", 
                                           meta_surfaceome, "</b><br>", "Replicate IDs: ", replicate_ids_left_right, "<br>" , "Feature IDs : ", features_left_right_overall, "<br>" ))
  ) +
    geom_segment(x=segmentation$x[1], y=segmentation$y[1], xend=segmentation$xend[1], yend=segmentation$yend[1], linetype="dashed", col="darkgrey", linewidth = 0.4) +
    geom_segment(x=segmentation$x[2], y=segmentation$y[2], xend=segmentation$xend[2], yend=segmentation$yend[2], linetype="dashed", col="darkgrey", linewidth = 0.4) +
    geom_segment(x=segmentation$x[3], y=segmentation$y[3], xend=segmentation$xend[3], yend=segmentation$yend[3], linetype="dashed", col="darkgrey", linewidth = 0.4) +
    geom_segment(x=segmentation$x[4], y=segmentation$y[4], xend=segmentation$xend[4], yend=segmentation$yend[4], linetype="dashed", col="darkgrey", linewidth = 0.4) +
    xlim(foldchangelimit_min, foldchangelimit_max) +
    geom_point(fill = c("#cc66ff", "#cc0000", "#0000EE", "#5F5F61")[match(subset_volcano$suvo_plot_color, c("surface+string", "surface", "string", "other"))],
               shape = ifelse(subset_volcano$imputed_comparison == "yes", 23, 21),
               alpha = subset_volcano$plot_alpha,
               size = 1.5,
               color = "white",
               stroke = 0.4) +
    # geom_label_repel(box.padding = 2, size=6, max.overlaps = 100000) +
    labs(fill="meta_surfaceome") +
    plot_theme() +
    # theme(axis.title.x = element_blank(),
    #       legend.position = "top",
    #       text = element_text(size = 14),
    #       axis.text = element_text(size = 12)) +
    xlab(paste("log2(fold change)"   , sep ="")) +
    ylab(paste("-log10(adj. p-value)", sep ="")) +
    geom_label(label= "Iso",
               x     = foldchangelimit_min+0.7, y=0,
               label.padding = unit(0.2, "lines"),
               label.size = 0.5,
               color = "white",
               fill  = "black",
               size  = 12/.pt) +
    geom_label(label= "TCR",
               x     =  foldchangelimit_max-1.2, y=0,
               label.padding = unit(0.2, "lines"),
               label.size = 0.5,
               color = "white",
               fill  = "black",
               size  = 12/.pt) 
    
  #  plot_volcano
  ggsave(
    filename = paste0("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_h", count_var, ".png"),
    plot = plot_volcano,
    width  = 8.00,  # 10.66,  # document is 16 cm wide             (before 12cm used)
    height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
    units  = "cm",
    dpi    = 300    # default for good quality
  )  
  
  count_var = count_var + 1
  
}

# Figure 2.3.1_j: Surfaceome contriburtion plots (Volcano sig-up) ..........................................................................................................................................................................................................................................
surfaceome_sigup_table <- v25_LLOQ_LUX_sigup %>%
  mutate(meta_surfaceome = case_when(entry_name %in% abTCR_chains[!abTCR_chains %in% c("CD3D_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN")] ~ "abTCR", TRUE ~ meta_surfaceome)) %>%
  # mutate(comparison = gsub("0.5", "0,5", comparison)) %>%
  group_by(comparison, meta_surfaceome) %>%
  tally() %>%
  tidyr::pivot_wider(
    names_from = meta_surfaceome,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(Surfaceome  = (yes    /sum(yes,no,abTCR))*100,
         Other       = (no     /sum(yes,no,abTCR))*100,
         abTCR       = (abTCR  /sum(yes,no,abTCR))*100,
         Surface_and_CD3_and_TCR = (sum(yes,abTCR)/sum(yes,no,abTCR))*100,
  )
#
ggplot(surfaceome_sigup_table %>%
         pivot_longer(cols = c(no, yes), names_to = "response", values_to = "count"),
       aes(x = comparison, y = count, fill = response)) +
  geom_bar(stat = "identity") +
  labs(title = "Surfaceome", x = "Comparison", y = "Count") +
  scale_fill_manual(values = c("no" = "#5F5F61", "yes" = "#cc0000")) +
  plot_theme()

Fig2.3.1_j <- surfaceome_sigup_table %>%
  pivot_longer(cols = c(Other, Surfaceome, abTCR), names_to = "Annotation", values_to = "count") %>%
  mutate(Annotation = factor(Annotation, levels = c("Other", "abTCR", "Surfaceome"))) %>%
  ggplot(aes(x = comparison, y = count, fill = Annotation)) +
  geom_bar(stat = "identity") +
  labs(title = "Enriched proteins", x = "panT cell input [e6]", y = "Enriched proteins [%]") +
  scale_fill_manual(values = c("Other" = "darkgrey", "Surfaceome" = "#cc0000", "abTCR" = "#0000EE")) +
  plot_theme()

Fig2.3.1_j
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_j.png",
  plot = Fig2.3.1_j,
  width  = 10.26,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Figure 2.3.1_k: sig up Upset (overlap) ..........................................................................................................................................................................................................................................
# Figure 2.3.1_l: sig up Upset for QC proteins (TCR chains, CD3, CD4, CD8 id overlap) ..........................................................................................................................................................................................................................................
upset_data <- bind_rows(
  # tibble(entry = v25_LLOQ_LUX_sigup %>% filter(comparison == "0.5")  %>% pull(entry) %>% unique(), 
  #        set = "0.5"),
  # tibble(entry = v25_LLOQ_LUX_sigup %>% filter(comparison == "1")  %>% pull(entry) %>% unique(), 
  #        set = "1") ,
  # tibble(entry = v25_LLOQ_LUX_sigup %>% filter(comparison == "5")  %>% pull(entry) %>% unique(), 
  #        set = "5")
  # ) %>% 
  tibble(entry = v25_LLOQ_LUX_sigup %>% filter(entry_name %in% abTCR_chains[!abTCR_chains %in% c("CD3D_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN")], comparison == "0.5")  %>% pull(entry) %>% unique(), 
           set = "0.5"),
    tibble(entry = v25_LLOQ_LUX_sigup %>% filter(entry_name %in% abTCR_chains[!abTCR_chains %in% c("CD3D_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN")], comparison == "1")  %>% pull(entry) %>% unique(), 
           set = "1") ,
    tibble(entry = v25_LLOQ_LUX_sigup %>% filter(entry_name %in% abTCR_chains[!abTCR_chains %in% c("CD3D_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN")], comparison == "5")  %>% pull(entry) %>% unique(), 
           set = "5")
  ) %>% 
  group_by(entry) %>% 
  summarise(sets = list(set)) %>%  # Critical: list column of set memberships
  ungroup() 

Fig2.3.1_k <- upset_data %>% 
  ggplot(aes(x = sets)) +
  geom_bar(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_upset(
    sets = c("0.5", "1", "5"),
    name = "",
    n_intersections = 20
  ) +
  labs(
    y     = "Intersection",
    title = "Enrichment Overlap"
  ) +
  plot_theme() +
  theme(axis.text.y = element_text(size = 14))

Fig2.3.1_k
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_k.png",
  plot = Fig2.3.1_k,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  


# Figure 2.3.1_l: Overlap/upset LUX-sig-up and CSC IDs ..........................................................................................................................................................................................................................................
upset_data <- bind_rows(
   tibble(entry = v25_LLOQ_LUX_sigup %>% pull(entry) %>% unique(), 
          set = "LUX"),
   tibble(entry = v25_LLOQ_CSC %>% pull(entry) %>% unique(), 
          set = "CSC"),
 ) %>% 
   group_by(entry) %>% 
   summarise(sets = list(set)) %>%  # Critical: list column of set memberships
   ungroup() 
 
 Fig2.3.1_l <- upset_data %>% 
   ggplot(aes(x = sets)) +
   geom_bar(fill = "black", color = "white", linewidth = 0.3) +
   scale_x_upset(
     sets = c("LUX", "CSC"),
     name = "",
     n_intersections = 20
   ) +
   labs(
     y     = "Intersection",
     title = "LUX to CSC Overlap"
   ) +
   plot_theme() +
   theme(axis.text.y = element_text(size = 14))

Fig2.3.1_l
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_l.png",
  plot = Fig2.3.1_l,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Figure 2.3.1_m: LUX sig-up vs. CSC signal ..........................................................................................................................................................................................................................................
v25_LLOQ_LUX_prot_nonzero <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_LUX_panT_lowinput_TCR/_output_1-2_ludo_adjp_0.7string_shs2.25/_data_prot_level.csv") %>%
  filter(normalised_protein_intensity_log2 > 0)

csc_subset_LUX_05 <- v25_LLOQ_CSC %>% 
  dplyr::filter(condition == "0.5", 
                normalised_protein_intensity_log2 > 0,
                entry %in% (v25_LLOQ_LUX_prot_nonzero %>% dplyr::filter(condition == "TCR_5e5") %>% pull(entry) %>% unique())) %>%
  mutate(comparison = "0.5",
         sig_up = case_when(entry %in% (v25_LLOQ_LUX_sigup %>% dplyr::filter(comparison == "0.5") %>% pull(entry) %>% unique()) 
                            ~ "Yes", TRUE ~ "No"))
csc_subset_LUX_1  <- v25_LLOQ_CSC %>% 
  dplyr::filter(condition == "1", 
                normalised_protein_intensity_log2 > 0,
                entry %in% (v25_LLOQ_LUX_prot_nonzero %>% dplyr::filter(condition == "TCR_1e6") %>% pull(entry) %>% unique())) %>%
  mutate(comparison = "1",
         sig_up = case_when(entry %in% (v25_LLOQ_LUX_sigup %>% dplyr::filter(comparison == "1") %>% pull(entry) %>% unique()) 
                            ~ "Yes", TRUE ~ "No"))
csc_subset_LUX_5  <- v25_LLOQ_CSC %>% 
  dplyr::filter(condition == "5",
                normalised_protein_intensity_log2 > 0,
                entry %in% (v25_LLOQ_LUX_prot_nonzero %>% dplyr::filter(condition == "TCR_5e6") %>% pull(entry) %>% unique())) %>%
  mutate(comparison = "5",
         sig_up = case_when(entry %in% (v25_LLOQ_LUX_sigup %>% dplyr::filter(comparison == "5") %>% pull(entry) %>% unique()) 
                            ~ "Yes", TRUE ~ "No"))

csc_subset_LUX <- rbind(csc_subset_LUX_05, csc_subset_LUX_1, csc_subset_LUX_5) %>%
  dplyr::select(condition, normalised_protein_intensity_log2, sig_up) %>%
  distinct()

Fig2.3.1_m <- ggplot(csc_subset_LUX, aes(x = condition, y = normalised_protein_intensity_log2, fill = sig_up)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "CSC signal for LUX IDs",
       x = "panT cell input [e6]",
       y = "log2(Norm. protein intensity) ",
       fill = "Sig-enriched") +
  scale_fill_manual(values = c("Yes" = "black", "No" = "darkgrey")) + # Customize colors as you like
  plot_theme()

Fig2.3.1_m
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_m.png",
  plot = Fig2.3.1_m,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  












######## end of functional section ________________________________________________


