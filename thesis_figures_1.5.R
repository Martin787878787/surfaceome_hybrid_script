#Load libraries required ----------------------------------------------------------------------------------------------------------
rm(list = ls())           # Purge workspace
set.seed(123)
script_version = "_1.5" # version stamp on output directory
# v1.5   enrichment based on confidence gradient IDs (naive CD4+ gets all panT t-test IDs added due to higher confidence)
# Core data wrangling and manipulation
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(magrittr)
# Plotting
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggthemes)
library(viridis)
library(gridExtra)
library(cowplot)
library(ggplotify)
library(plotly)
library(htmlwidgets)
# Overlap and set visualization
library(ggupset)
# Proteomics and differential analysis
library(protti)

sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/shs_functions", pattern = "\\.R$", full.names = TRUE), source)
sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/thesis_figures_functions", pattern = "\\.R$", full.names = TRUE), source)

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
# ***************************************************************************************************************************************************************
# *****************************************************************************************************************************************************************

##################################################################################################################################################################################################################
## GENERAL RESOURCES ####################################################################################################################################################################################################
##################################################################################################################################################################################################################
CSPA <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_per_cell_type.csv") %>%
  filter(protein_count > 0) %>%
  dplyr::rename(entry = id_link) %>%
  pull(entry) %>%
  unique()
paste("CSPA comprises", length(CSPA), "proteins")
CSPA_Jurkat <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/CSPA_per_cell_type.csv") %>%
  filter(jurkat > 0) %>%
  dplyr::rename(entry = id_link) %>%
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


##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
#   ____ _                 _              ____  
#  / ___| |__   __ _ _ __ | |_ ___ _ __  |___ \ 
# | |   | '_ \ / _` | '_ \| __/ _ \ '__|   __) |
# | |___| | | | (_| | |_) | ||  __/ |     / __/ 
# \____|_| |_|\__,_| .__/ \__\___|_|    |_____|
#                  |_|                         
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################

# Figure 2.3.1: CSC LLOQ panT  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
v25_LLOQ_CSC <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_CSC_panT_lowinput/_output_1-2_ludo_adjp_0.7string_shs2.27/_data_prot_level.csv") %>%
  mutate(condition = gsub("0_5", "0.5", condition)) 

v25_LLOQ_CSC_pep <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_CSC_panT_lowinput/_output_1-2_ludo_adjp_0.7string_shs2.27/_data_feature_level.csv") %>%
  mutate(condition = gsub("0_5", "0.5", condition)) %>%
  mutate(N_deam_psm_n = str_count(peptide_sequence_mod,   "N\\[0.9840\\]"))

table(v25_LLOQ_CSC_pep %>% dplyr::filter(csc_signature_psm == "yes", !is.na(ng_sites)) %>%  dplyr::select(entry, csc_signature_psm_n, peptide_sequence_mod, ng_sites) %>% distinct() %>% pull(csc_signature_psm_n))

# QC of deamidation signatures - how many in NxSTC signature - should be majority otherwise likely chem deam possible too
NsiteQC <- v25_LLOQ_CSC_pep %>% dplyr::select(peptide_sequence_mod, N_deam_psm_n, csc_signature_psm_n) %>% 
  distinct() %>% # keep only unique modified peptides (no charge states)
  dplyr::select(-peptide_sequence_mod) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Count")
  
ggplot(NsiteQC %>%
         group_by(Type) %>%
         summarise(Sum = sum(Count)), aes(x = Type, y = Sum, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("N_deam_psm_n" = "#ff7f0e", "csc_signature_psm_n" = "black")) +
  plot_theme()

bad_csc_sites   <- v25_LLOQ_CSC_pep %>% filter(csc_signature_psm == "no" , N_deam_psm_n >0) %>% dplyr::select(entry_name, semi_stripped_peptide_sequence, meta_surfaceome, cspa_2015) %>% distinct()
weird_csc_sites <- v25_LLOQ_CSC_pep %>% filter(csc_signature_psm == "yes", N_deam_psm_n >1) %>% # require at least 2 sites (CSC and non-csc site)
  dplyr::select(entry_name, semi_stripped_peptide_sequence, meta_surfaceome, cspa_2015) %>% distinct()
table(bad_csc_sites$meta_surfaceome)
table(weird_csc_sites$meta_surfaceome)
table(bad_csc_sites$cspa_2015)
table(weird_csc_sites$cspa_2015)
##### ________________________________________________________________________________________________________________________________________________________________________________________________
######################################################################################################################################################################################################

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
    fill  = "CSC\nsignature"  ) +
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

paste("Cummulative (replicates) CSC peptides 0.5: ", length(v25_LLOQ_CSC_pep %>% filter(condition == "0.5"  , csc_signature_psm == "yes") %>% pull(peptide_sequence_mod) %>% unique()))
paste("Cummulative (replicates) CSC peptides 1: "  , length(v25_LLOQ_CSC_pep %>% filter(condition == "1"    , csc_signature_psm == "yes") %>% pull(peptide_sequence_mod) %>% unique()))
paste("Cummulative (replicates) CSC peptides 5: "  , length(v25_LLOQ_CSC_pep %>% filter(condition == "5"    , csc_signature_psm == "yes") %>% pull(peptide_sequence_mod) %>% unique()))

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
    fill  = "CSC\nsignature"  ) +
  plot_theme()

paste("Cummulative (replicates) CSC proteins 0.5: ", length(v25_LLOQ_CSC_pep %>% filter(condition == "0.5", csc_signature_psm == "yes") %>% pull(entry) %>% unique()))
paste("Cummulative (replicates) CSC proteins 1: "  , length(v25_LLOQ_CSC_pep %>% filter(condition == "1"    , csc_signature_psm == "yes") %>% pull(entry) %>% unique()))
paste("Cummulative (replicates) CSC proteins 5: "  , length(v25_LLOQ_CSC_pep %>% filter(condition == "5"    , csc_signature_psm == "yes") %>% pull(entry) %>% unique()))


Fig2.3.1_c 
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_c.png",
  plot   = Fig2.3.1_c,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)

# Fig2.3.1_X  Fig2.3.1_d glyco sites/peptides, glyco glyco-peptides/protein ......................................................................................................................
v25_LLOQ_CSC_pep_csc <- v25_LLOQ_CSC_pep %>%
  dplyr::select(entry, entry_name, condition, csc_signature_psm, peptide_sequence_mod) %>%
  distinct() %>%
  filter(csc_signature_psm == "yes") %>%
  mutate(csc_signature_count = str_count(peptide_sequence_mod, fixed("N[0.9840]"))) %>%   # no stringent filter required because signature prefiltered in line above
  group_by(entry, condition) %>%
  mutate(csc_peptides_per_protein = n_distinct(peptide_sequence_mod)) %>%
  ungroup()

table(v25_LLOQ_CSC_pep_csc$csc_signature_count) # for most the numer is 1 (1746) only for 13,3 % it is 2 (232)
table(v25_LLOQ_CSC_pep_csc %>% filter(condition == "5") %>% pull(csc_signature_count)) # for most the numer is 1 (1746) only for 13,3 % it is 2 (232)
summary(v25_LLOQ_CSC_pep_csc$csc_peptides_per_protein)
paste("median sites / protein for 5 e6 panT input:", median(v25_LLOQ_CSC_pep_csc %>% filter(condition == "5")   %>% pull(csc_peptides_per_protein)))
paste("median sites / protein for 1 e6 panT input:", median(v25_LLOQ_CSC_pep_csc %>% filter(condition == "1")   %>% pull(csc_peptides_per_protein)))
paste("median sites / protein for 5 e5 panT input:", median(v25_LLOQ_CSC_pep_csc %>% filter(condition == "0.5") %>% pull(csc_peptides_per_protein)))


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

# Fig2.3.1_d glyco-peptides/protein ......................................................................................................................
Fig2.3.1_d <- ggplot(v25_LLOQ_CSC_pep_csc, aes(x = condition, y = csc_peptides_per_protein)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "#4e79a7") +
  labs(
    title = "Peptides/protein",
    x     = "panT cell input [e6]",
    y     = "Number of CSC peptides"
  ) +
  coord_cartesian(ylim = c(0, 30)) +  # Limits y-axis to 0–40
  plot_theme()  # Use built-in theme instead of plot_theme() until verified

Fig2.3.1_d
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_d.png",
  plot = Fig2.3.1_d,
  width  = 7.5,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_e protein CVs ......................................................................................................................
Fig2.3.1_e <- ggplot(v25_LLOQ_CSC, aes(x = condition, y = cv_protein, fill = condition)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "Protein CVs",
    x     = "panT cell input [e6]",
    y     = "Protein CV [%]"  ) +
  # scale_fill_manual(values = c("#cc3300", "#ff9900", "#009933")) +
  scale_fill_manual(values = c("#ffffff", "#808080", "#4d4d4d")) +
  plot_theme() +
  theme(legend.position = "none") # dont show legend - already x axis contains the info

Fig2.3.1_e
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_e.png",
  plot = Fig2.3.1_e,
  width  = 8,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

paste("median CV 5e6 cells CSC:  ", round(v25_LLOQ_CSC %>% filter(condition == "5") %>% dplyr::select(entry, cv_protein) %>% distinct() %>% pull(cv_protein) %>% median(na.rm = TRUE),1))
paste("median CV 1e6 cells CSC:  ", round(v25_LLOQ_CSC %>% filter(condition == "1") %>% dplyr::select(entry, cv_protein) %>% distinct() %>% pull(cv_protein) %>% median(na.rm = TRUE),1))
paste("median CV 5e5 cells CSC:  ", round(v25_LLOQ_CSC %>% filter(condition == "0.5") %>% dplyr::select(entry, cv_protein) %>% distinct() %>% pull(cv_protein) %>% median(na.rm = TRUE),1))

# Fig2.3.1_f: CSC Precursor Signal ......................................................................................................................
Fig2.3.1_f <- v25_LLOQ_CSC_pep %>%
  filter(csc_signature != "cont") %>%
  group_by(con_rep, csc_signature) %>%
  mutate(TotInt = sum(raw_prec_intensity)) %>%
  dplyr::select(condition, bio_replicate, con_rep, csc_signature, TotInt) %>%
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
    title = "MS signal",
    x     = "panT cell input [e6]",
    y     = "CSC feature signal [%]",
    fill  = "CSC signature"  ) +
  plot_theme()
  
Fig2.3.1_f
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_f.png",
  plot = Fig2.3.1_f,
  width  = 7.5,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_g: CSC Protein ID data completeness ......................................................................................................................
Fig2.3.1_g <- v25_LLOQ_CSC_pep %>%
  filter(csc_signature == "yes", raw_prec_intensity > 0) %>%
  dplyr::select(condition, entry, bio_replicate) %>%
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
  labs(title = "Data completeness",
       x     = "panT cell input [e6]",
       y     = "Protein count",
       fill  = "Detection [%]") +
  scale_fill_manual(values = c("black", "#7F7F7F", "#CCCCCC")) +
  plot_theme() +
  theme(legend.position = "right")

Fig2.3.1_g
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_g.png",
  plot = Fig2.3.1_g,
  width  = 10,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_h overlap/upset Jurkat CSPA & pT CSC ......................................................................................................................
# write.csv(data.frame(Values = CSPA_Jurkat), "/Users/mgesell/Downloads/CSPA_Jurkat.csv", row.names = FALSE)
# Prepare data in required format
panT_csc_ids <- v25_LLOQ_CSC_pep %>% 
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

Fig2.3.1_h <- upset_data %>% 
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

Fig2.3.1_h
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_h.png",
  plot = Fig2.3.1_h,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  
# fyi
cat("\n",
    "Novel on panT & CSPA supported ", length(setdiff(intersect(panT_csc_ids, CSPA), CSPA_Jurkat))  , "\n",
    "Overlap panT-CSPA              ", length(intersect(panT_csc_ids, CSPA)), "\n",
    "Jurkat but not panT            ", length(setdiff(CSPA_Jurkat, panT_csc_ids))  , "\n",
    "Novel panT                     ", length(setdiff(panT_csc_ids, CSPA))
)

# info to write in thesis
paste("Unique overall CSC PEPTIDES", v25_LLOQ_CSC_pep %>% filter(csc_signature == "yes") %>% pull(peptide_sequence_mod) %>% unique() %>% length() )
paste("Unique overall CSC PROTEINS", v25_LLOQ_CSC_pep %>% filter(csc_signature == "yes") %>% pull(entry) %>% unique() %>% length() )
paste("Unique 5 e6 CSC PROTEINS", v25_LLOQ_CSC_pep %>% filter(csc_signature == "yes", condition == "5") %>% pull(entry) %>% unique() %>% length() )
paste("Consider reporting median CV values at least in text")


# Fig2.3.1_i map CSC sites to glyco_upsp_202501_list  ......................................................................................................................
evaluate_csc_sites_result <- evaluate_csc_sites(csc_output__data_raw_up_surf_gly_nZ = v25_LLOQ_CSC_pep)
csc_sites_plus_upsp_annotation       <- evaluate_csc_sites_result[[1]]
novel_Ng_proteins_surface_annotated  <- evaluate_csc_sites_result[[2]] %>% left_join(proteome_upsp_202501 %>% dplyr::select(entry, entry_name), by = "entry")
cat(novel_Ng_proteins_surface_annotated$entry %>% unique(), sep = "\n")
#
query_protein = "P09486"
query_protein_peptides <- extract_csc_peptides_for_protti(csc_output__data_raw_up_surf_gly_nZ = v25_LLOQ_CSC_pep,
                                                          csc_sites_plus_upsp_annotation      = csc_sites_plus_upsp_annotation, 
                                                          query_protein                       = query_protein)
#
Ng_summary <- csc_sites_plus_upsp_annotation %>%
  dplyr::select(entry, site, summary) %>%
  distinct() %>%
  mutate(summary = gsub("NOVEL site", "novel", summary)) %>%
  dplyr::rename(CSC_site = site) %>%
  count(summary) %>% 
  mutate(summary = factor(summary, levels = c("known", "validated", "novel"))) 
Ng_summary
#
Fig2.3.1_i <- Ng_summary %>%
  #
  ggplot(aes(x = NA, y = n, fill = summary)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("novel" = "#0000cc", "validated" = "#4e79a7", "known" = "grey")) +
  labs(y = "Count", x = "", fill = "Summary") +
  labs(
    title = "N-glyco sites",
    x     = "Overall",
    y     = "Glyco site count",
    fill  = "Category"  ) +
  plot_theme()+   # Hide x-axis text)
  theme(axis.text.x = element_blank())
Fig2.3.1_i
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_i.png",
  plot = Fig2.3.1_i,
  width  = 9,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  


# ===================================================================================================================================================================================================================================================================================
# Figure 2.3.1_k: TCR-LUX LLOQ panT  ..........................................................................................................................................................................................................................................
v25_LLOQ_LUX_data_prot_diff_abundance <- read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_LUX_panT_lowinput_TCR/_output_1-2_ludo_adjp_0.7string_shs2.25/_data_prot_diff_abundance.csv") %>%
  arrange(comparison ) %>%
  group_by(entry_name) %>%
  mutate(overlap = paste(comparison, collapse = "_&_")) %>%
  ungroup()
# sig up data
v25_LLOQ_LUX_sigup <- v25_LLOQ_LUX_data_prot_diff_abundance %>%
  filter(log2FC >1, adj_pvalue <= 0.05) %>% # sig up 
  dplyr::select(entry, entry_name, comparison, overlap) %>%
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
filter_log2fc_cutoff        = log(2, 2)   # Cutoff for Volcano plotting
pvalue                      = "adj_pvalue"          # select "pvalue" or "adj_pvalue" depending on how stringend you want to be.   recommendation for CSC is adj_pvalue    
filter_sig_cutoff           = 0.05  
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
  
    foldchangelimit_max <- max(v25_LLOQ_LUX_data_prot_diff_abundance$log2FC, na.rm = TRUE)
    foldchangelimit_min <- min(v25_LLOQ_LUX_data_prot_diff_abundance$log2FC, na.rm = TRUE)
    
    
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
      filename = paste0("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_k", count_var, ".png"),
      plot = plot_volcano,
      width  = 8.00,  # 10.66,  # document is 16 cm wide             (before 12cm used)
      height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
      units  = "cm",
      dpi    = 300    # default for good quality
    )  
    
    count_var = count_var + 1
}

# Figure Fig2.3.1_l: Surfaceome contriburtion plots (Volcano sig-up) ..........................................................................................................................................................................................................................................
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
Fig2.3.1_l <- surfaceome_sigup_table %>%
  pivot_longer(cols = c(Other, Surfaceome, abTCR), names_to = "Annotation", values_to = "count") %>%
  mutate(Annotation = factor(Annotation, levels = c("Other", "abTCR", "Surfaceome"))) %>%
  ggplot(aes(x = comparison, y = count, fill = Annotation)) +
  geom_bar(stat = "identity") +
  labs(title = "LUX annotation", 
       x     = "panT cell input [e6]", 
       y     = "Significantly\nenriched proteins [%]",
       fill  = "Category") +
  scale_fill_manual(values = c("Other" = "darkgrey", "Surfaceome" = "#cc0000", "abTCR" = "#0000EE")) +
  plot_theme()

Fig2.3.1_l
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_l.png",
  plot = Fig2.3.1_l,
  width  = 11,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Figure Fig2.3.1_: sig up Upset (overlap) ..........................................................................................................................................................................................................................................

# Figure Fig2.3.1_m: sig up Upset for QC proteins (TCR chains, CD3, CD4, CD8 id overlap) ..........................................................................................................................................................................................................................................
upset_data <- bind_rows(
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

Fig2.3.1_m <- upset_data %>% 
  ggplot(aes(x = sets)) +
  geom_bar(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_upset(
    sets = c("0.5", "1", "5"),
    name = "",
    n_intersections = 20
  ) +
  labs(
    y     = "Intersection",
    title = "abTCR chains"
  ) +
  plot_theme() +
  theme(axis.text.y = element_text(size = 14))

Fig2.3.1_m
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_m.png",
  plot = Fig2.3.1_m,
  width  = 7,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_n: sig up Upset for QC proteins (TCR chains, CD3, CD4, CD8 id overlap) ..........................................................................................................................................................................................................................................
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
  filter(protein1_entry_name %in% string_filter   |    protein2_entry_name %in% string_filter  )

# string interactor mapping
string_interactors <- string %>% filter(protein1_entry_name %in% string_targets$TCR | protein2_entry_name %in% string_targets$TCR) %>% 
  dplyr::select(protein1_entry_name, protein2_entry_name) %>% 
  unlist() %>%
  unique()

upset_data <- bind_rows(
  tibble(entry = v25_LLOQ_LUX_sigup %>% filter(entry_name %in% string_interactors, comparison == "0.5")  %>% pull(entry) %>% unique(), 
         set = "0.5"),
  tibble(entry = v25_LLOQ_LUX_sigup %>% filter(entry_name %in% string_interactors, comparison == "1")  %>% pull(entry) %>% unique(), 
         set = "1") ,
  tibble(entry = v25_LLOQ_LUX_sigup %>% filter(entry_name %in% string_interactors, comparison == "5")  %>% pull(entry) %>% unique(), 
         set = "5")
) %>% 
  group_by(entry) %>% 
  summarise(sets = list(set)) %>%  # Critical: list column of set memberships
  ungroup() 

Fig2.3.1_n <- upset_data %>% 
  ggplot(aes(x = sets)) +
  geom_bar(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_upset(
    sets = c("0.5", "1", "5"),
    name = "",
    n_intersections = 20
  ) +
  labs(
    y     = "Intersection",
    title = "STRING Interactors"
  ) +
  plot_theme() +
  theme(axis.text.y = element_text(size = 14))

Fig2.3.1_n
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_n.png",
  plot = Fig2.3.1_n,
  width  = 9,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_o: Overlap/upset LUX-sig-up and CSC IDs ..........................................................................................................................................................................................................................................
upset_data <- bind_rows(
   tibble(entry = v25_LLOQ_LUX_sigup %>% filter(comparison == "1") %>% pull(entry) %>% unique(), 
          set = "LUX"),
   tibble(entry = v25_LLOQ_CSC %>% pull(entry) %>% unique(), 
          set = "CSC"),
 ) %>% 
   group_by(entry) %>% 
   summarise(sets = list(set)) %>%  # Critical: list column of set memberships
   ungroup() 
 
 Fig2.3.1_o <- upset_data %>% 
   ggplot(aes(x = sets)) +
   geom_bar(fill = "black", color = "white", linewidth = 0.3) +
   scale_x_upset(
     sets = c("LUX", "CSC"),
     name = "",
     n_intersections = 20
   ) +
   labs(
     y     = "Intersection",
     title = "1 e6 LUX sig-up\nvs. CSC signal"
   ) +
   plot_theme() +
   theme(axis.text.y = element_text(size = 14))

Fig2.3.1_o
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_o.png",
  plot = Fig2.3.1_o,
  width  = 8,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# Fig2.3.1_p: LUX sig-up vs. CSC signal ..........................................................................................................................................................................................................................................
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

Fig2.3.1_p <- ggplot(csc_subset_LUX, aes(x = condition, y = normalised_protein_intensity_log2, fill = sig_up)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "CSC signal for LUX IDs",
       x     = "panT cell input [e6]",
       y     = "log2(Normalized\nCSC protein intensity) ",
       fill  = "Significantly\nenriched") +
  scale_fill_manual(values = c("Yes" = "black", "No" = "darkgrey")) + # Customize colors as you like
  plot_theme()

Fig2.3.1_p
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig2.3.1_p.png",
  plot = Fig2.3.1_p,
  width  = 10.66,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  


##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
#     ____ _                  _            _____ 
#   / ___| |__   __ _ _ __  | |_ ___ _ __ |___ / 
#  | |   | '_ \ / _` | '_ \| __/ _ \ '__|   |_ \ 
#  | |___| | | | (_| | |_) | ||  __/ |    ___) |
#  \____|_| |_|\__,_| .__/ \__\___|_|    |____/ 
#                   |_|                         
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
# LUX =========================================================================================================================================================================================================================
v31_LUX_data_prot   <- rbind( read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.27_4pTss/_data_prot_level.csv"),
                              read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.27_ss_meta/_data_prot_level.csv"),
                              read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.27_full_meta/_data_prot_level.csv")) %>%  
  as.data.frame() %>%
  dplyr::select(-X) # determine comparison overlap

v31_LUX_data_prot_diff_abundance   <- rbind( read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.27_4pTss/_data_prot_diff_abundance.csv"),
                         read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.27_ss_meta/_data_prot_diff_abundance.csv"),
                         read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.27_full_meta/_data_prot_diff_abundance.csv")) %>%  
  as.data.frame() %>%
  dplyr::select(-X) # determine comparison overlap

# bind and process individually processed comparisons
v31_LUX_data_prot_diff_abundance <- v31_LUX_data_prot_diff_abundance %>%
  mutate(plot_heading = case_when(ttest_condition == "nCD4_TCR"    ~ "Naive CD4+",
                                  ttest_condition == "nCD8_TCR"    ~ "Naive CD8+",
                                  ttest_condition == "nnCD4_TCR"   ~ "Memory CD4+",
                                  ttest_condition == "nnCD8_TCR"   ~ "Memory CD8+",
                                  ttest_condition == "CD4meta_TCR" ~ "Meta CD4+",
                                  ttest_condition == "CD8meta_TCR" ~ "Meta CD8+",
                                  ttest_condition == "nMeta_TCR"   ~ "Meta Naive",
                                  ttest_condition == "nnMeta_TCR"  ~ "Meta Memory",
                                  ttest_condition == "metaTCR"     ~ "Meta panT",
                                  TRUE ~ "error"),
         plot_right   = case_when(ttest_condition == "nCD4_TCR"    ~ "TCR",
                                  ttest_condition == "nCD8_TCR"    ~ "TCR",
                                  ttest_condition == "nnCD4_TCR"   ~ "TCR",
                                  ttest_condition == "nnCD8_TCR"   ~ "TCR",
                                  ttest_condition == "CD4meta_TCR" ~ "TCR",
                                  ttest_condition == "CD8meta_TCR" ~ "TCR",
                                  ttest_condition == "nMeta_TCR"   ~ "TCR",
                                  ttest_condition == "nnMeta_TCR"  ~ "TCR",
                                  ttest_condition == "metaTCR"     ~ "TCR",
                                  TRUE ~ "error"),
         plot_left   = case_when(ttest_condition  == "nCD4_TCR"    ~ "Iso",
                                  ttest_condition == "nCD8_TCR"    ~ "Iso",
                                  ttest_condition == "nnCD4_TCR"   ~ "Iso",
                                  ttest_condition == "nnCD8_TCR"   ~ "Iso",
                                  ttest_condition == "CD4meta_TCR" ~ "Iso",
                                  ttest_condition == "CD8meta_TCR" ~ "Iso",
                                  ttest_condition == "nMeta_TCR"   ~ "Iso",
                                  ttest_condition == "nnMeta_TCR"  ~ "Iso",
                                  ttest_condition == "metaTCR"     ~ "Iso",
                                  TRUE ~ "error")
         )

# sig up data
v31_LUX_data_prot_diff_abundance_sigup <- v31_LUX_data_prot_diff_abundance %>%
  filter(log2FC >1, adj_pvalue <= 0.05) %>% # sig up 
  dplyr::select(entry, entry_name, plot_heading, comparison) %>% # comparison
  mutate(LUX_id = 1  )  %>%
  mutate(meta_surfaceome = ifelse(entry %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "yes", "no")) %>%
  arrange(plot_heading) %>%  
  group_by(entry_name) %>% 
  mutate(overlap           = paste(sort(unique(plot_heading)), collapse = "||"),
         proximity_id_cout = str_count(overlap, "\\|\\|")+1) %>%  
  ungroup() 
message(length(unique(v31_LUX_data_prot_diff_abundance_sigup %>% filter(plot_heading == "Meta panT") %>% pull(entry) )), " unique sig-enriched proteins in panT core")
message(length(unique(v31_LUX_data_prot_diff_abundance_sigup %>% filter(plot_heading %in% c("Meta CD4+", "Meta CD8+", "Meta Naive", "Meta Memory", "Meta panT")) %>% pull(entry) )), " unique sig-enriched proteins in the 5 meta subsets")
message(length(unique(v31_LUX_data_prot_diff_abundance_sigup %>% filter(plot_heading %in% c("Naive CD4+", "Naive CD8+", "Memory CD4+", "Memory CD8+"))           %>% pull(entry) )), " unique sig-enriched proteins in the 4 pT subsets"  )
message(length(unique(v31_LUX_data_prot_diff_abundance_sigup$entry)), " unique sig-enriched proteins accros all comparisons")

# core known TCR community
abTCR_chains <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv") %>% 
  filter(abTCR_chains_cd3_mgmanual != "") %>% 
  pull(abTCR_chains_cd3_mgmanual)
TCR_LUX_pois <- c(c("CD4_HUMAN", "CD8A_HUMAN" , "CD8B_HUMAN"), abTCR_chains)

# Figure 3.3.1_a: abTCR-LUX  ..........................................................................................................................................................................................................................................
# consider to only display red grey blue purple annotation (surface / non-surface, string-interactors, string&surface )
filter_log2fc_cutoff        = log(2, 2)   # Cutoff for Volcano plotting
pvalue                      = "adj_pvalue"          # select "pvalue" or "adj_pvalue" depending on how stringend you want to be.   recommendation for CSC is adj_pvalue    
filter_sig_cutoff           = 0.05  
#
plot_no_df_occupacy         =  0.7      # make lower confidence ids less bright
plot_label_volcano          =  "poi"               # default poi   
plot_label_inf              =  plot_label_volcano  # default plot_label_volcano
plot_fill_volcano           =  "meta_surfaceome"   # default meta_surfaceome  
plot_fill_label_inf         =  plot_fill_volcano   # default plot_fill_volcano                 
plot_fill_condition_scatter =  "meta_surfaceome" # x-y median condition signal plot (color surface or poi recommended)
#
count_var = 1
loop_frame <- unique(v31_LUX_data_prot_diff_abundance$comparison)
regenerate_volcanos = FALSE
if (regenerate_volcanos == TRUE) { # comp_counter = 1
  for (comp_counter in 1:length(loop_frame)) {   # ----------------------- LOOOOOOOP --------------------- LOOOOOOOP ------------------------------ LOOOOOOOP ------------------- LOOOOOOOP -------------------
    # usedConditions <- c(gsub("_vs_.*", "", loop_frame[comp_counter]), gsub("^.*_vs_", "", loop_frame[comp_counter]))
    result_subset <- v31_LUX_data_prot_diff_abundance[v31_LUX_data_prot_diff_abundance$comparison == loop_frame[comp_counter],] 
    # adjust p value (<<-- after this point no more data filtering)
    result_subset$adj_pvalue = p.adjust(result_subset$pvalue, method="BH")
    result_subset$pvalue      [result_subset$pvalue        < 2.220446e-16 ] = 2.220446e-16
    result_subset$adj_pvalue  [result_subset$adj_pvalue    < 2.220446e-16 ] = 2.220446e-16
    result_subset <- result_subset[order(result_subset$adj_pvalue),]  # sort by significance - good on top
    result_subset <- result_subset[!is.na(result_subset$log2FC), ]   # filter out all proteins that were not identified - all other should have calculated FC or +/-Inf assigned
    # used to define poi col below in result_subset
    poi_var <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv" , header = TRUE, sep = ",") %>%
      filter(abTCR_chains_cd3_mgmanual != "") %>%
      pull(abTCR_chains_cd3_mgmanual)
    # annoated surface and poi
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
    
    foldchangelimit_max <- max(v31_LUX_data_prot_diff_abundance$log2FC, na.rm = TRUE)
    foldchangelimit_min <- min(v31_LUX_data_prot_diff_abundance$log2FC, na.rm = TRUE)
  
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
      geom_label(label = unique(subset_volcano$plot_heading),
                 x      = -Inf,           # Places at left edge of plot
                 y      =  Inf,           # Places at top edge of plot
                 hjust  = 0,              # Aligns box to the left of (x, y)
                 vjust  = 1,              # Aligns box to the top of (x, y)
                 label.padding = unit(0.2, "lines"),
                 label.size = 0.5,
                 color = "white",
                 fill  = "black",
                 size  = 12/.pt) +
      geom_label(label = unique(subset_volcano$plot_left),
                 x     = foldchangelimit_min+0.7, y=0,
                 label.padding = unit(0.2, "lines"),
                 label.size = 0.5,
                 color = "black",
                 fill  = "grey",
                 size  = 12/.pt) +
      geom_label(label= unique(subset_volcano$plot_right),
                 x     =  foldchangelimit_max-1.2, y=0,
                 label.padding = unit(0.2, "lines"),
                 label.size = 0.5,
                 color = "black",
                 fill  = "grey",
                 size  = 12/.pt) 
    
    #  plot_volcano
    ggsave(
      filename = paste0("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.1_a", count_var, ".png"),
      plot = plot_volcano,
      width  = 8.00,  # 10.66,  # document is 16 cm wide             (before 12cm used)
      height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
      units  = "cm",
      dpi    = 300    # default for good quality
    )  
    
    count_var = count_var + 1
  }
}
# __________________________________________________________________________________________________________________________________________________________________________________________________________________
# Fig3.3.1_b: abTCR-LUX  ..........................................................................................................................................................................................................................................

upset_data <- bind_rows( # ""  ""  "" "" ""   ""   ""   ""  ""
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Meta panT")  %>% pull(entry) %>% unique(),      set = "Meta panT"),
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Meta CD4+")  %>% pull(entry) %>% unique(),      set = "Meta CD4+"),
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Meta CD8+")  %>% pull(entry) %>% unique(),      set = "Meta CD8+"),
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Meta Naive") %>% pull(entry) %>% unique(),      set = "Meta Naive"),
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Meta Memory")%>% pull(entry) %>% unique(),      set = "Meta Memory"),
  ) %>% 
  group_by(entry) %>% 
  summarise(sets = list(set)) %>%  # Critical: list column of set memberships
  ungroup() 

Fig3.3.1_b <- upset_data %>% 
  ggplot(aes(x = sets)) +
  geom_bar(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_upset(
    sets = unique(v31_LUX_data_prot_diff_abundance_sigup$plot_heading),
    name = "",
    n_intersections = 13  ) +
  labs(y     = "Intersection",
       title = "Protein Community Intersection") +
  plot_theme() # +
  # theme(axis.text.y = element_text(size = 14))

Fig3.3.1_b
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.1_b.png",
  plot = Fig3.3.1_b,
  width  = 15,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  


######################################################################################################################################################################################################################### ########################################################
#### enrichment analysis - go, complex, family, domain, ... ######################################################################################################################################################################################################################
# protein family enrichment analysis =====================================================================================================================================================
sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/thesis_figures_functions", pattern = "\\.R$", full.names = TRUE), source)
proteome_20250804 <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/human_upsp_202508.csv", header = TRUE) %>% setNames(gsub("\\.", "_", tolower(names(.)))) %>% setNames(gsub("__", "_", names(.))) %>% setNames(gsub("_$", "",  names(.)))

# extend the dataset so each subset also contains panT information --> full set and set-panT proteins for correct determination of enrichment 
v31_LUX_data_prot_diff_abundance_sigup_upsp <- v31_LUX_data_prot_diff_abundance_sigup %>%
  mutate(plot_heading = gsub("Meta ", "", plot_heading)) %>%
  dplyr::select(entry, plot_heading) %>%
  left_join(proteome_20250804 %>% dplyr::select(entry, protein_families, domain_cc, domain_ft, drugbank, drugcentral), by = "entry")

# Figure 3.3.xxx: GO term bubble plot  ..........................................................................................................................................................................................................................................
library(gprofiler2)      # GO gost 
library(org.Hs.eg.db)    # GO gost
library(biomaRt)         # GO gost

# perform GO enrichment for all 
gost_LUX <- NULL
for (i in unique(v31_LUX_data_prot_diff_abundance_sigup_upsp$plot_heading)) {
  go_result <- go_gost(query_list    = unique(v31_LUX_data_prot_diff_abundance_sigup_upsp %>% filter(plot_heading == i) %>% pull(entry)),
                       set = i,
                       max_term_size = 100,
                       # # NULL = all GO branches;       a vector of data sources to use. Currently, these include GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP. Please see the g:GOSt web tool for the comprehensive list and details on incorporated data sources.
                       go_source     = c("GO:BP", "GO:MF", "GO:CC", # basic
                                         "CORUM",                   # complexes
                                         "KEGG","REAC", "WP",       # signaling 
                                         "HPA", "HP", "HPA" #,      # phenotype (human phenotype atlas) 
                                         # "TF", "MIRNA"            # regulatory DNA motifs
                                         )   
                       )                    
  gost_LUX <- rbind(gost_LUX,go_result)
}
# issue with some complexes from CORUM that AB+C complexes are cosidered to have 2 components - which is wrong
# >> manually correct (cant think of universal fix right now)
gost_LUX <- gost_LUX %>%
  mutate(ID_count    = str_count(intersection_gene,",")+1)
# combine GO sets to one long df
gost_LUX <- gost_LUX %>%
  group_by(term_name) %>%
  mutate(overlap = paste(sort(unique(comparison)), collapse = "_")) %>%
  ungroup() %>%
  mutate(term_name =  str_replace_all(term_name,  # remove redundant complex terms
                                      c("alpha-beta"    = "ab", # catch before downstream modifies ab T cell activaiton and similar
                                        "ITG"           = "IT",
                                        "IT"            = "ITG",
                                        "integrin alpha"= "ITGA",
                                        "integrin beta" = "ITGB",
                                        "alpha"         = "ITGA" , 
                                        "-beta"         = "-ITGB" , 
                                        "interleukin "  = "IL",
                                        "ab T cell"     = "alpha-beta T cell" # correct
                                      )),
         term_name = str_replace(term_name, "ITGcomplex", "integrin complex")) %>%
  distinct(term_name, comparison, recall, .keep_all = TRUE)  %>% # remove redundant complex terms 
  mutate(term_size   = case_when(term_name == "ITGA5-ITGB1-SPP1 complex"   ~ 3,
                                 term_name == "ITGAV-ITGB1-SPP1 complex"   ~ 3,
                                 term_name == "ITGA4-ITGB1-CD81 complex"   ~ 3,
                                 term_name == "ITGA1-ITGB1-COL6A3 complex" ~ 3,
                                 term_name == "ITGA1-ITGB1-PTPN2 complex"  ~ 3,
                                 term_name == "ITGA4-ITGB1 complex"  ~ 2,
                                 TRUE ~term_size),
         recall       = ID_count/term_size) %>%         
  dplyr::filter(!term_name %in% c("ITGAv-ITGB1 complex",    # duplicated name
                                  "ITGA5-ITGB4 complex",    # wrong annotated in corum - acutally ITGA5-ITGB1 (contain in other term category - so removed corum term)
                                  "ITGA3-ITGB1 complex")  ) # redundant through ITGA3-ITGB1-BSG complex 

## plot GO result =====================================================================================================================================================
# Define order and elements to that are to be plotted together 
group_filters <- list( 
  c("panT", "CD4+" , "CD8+" , "Naive",  "Memory", "Naive CD4+", "Memory CD4+", "Naive CD8+" , "Memory CD8+" ),# NULL,
  c("panT", "CD4+" , "CD8+" , "Naive",  "Memory"),
  c("Naive CD4+", "Memory CD4+", "Naive CD8+" , "Memory CD8+")               )
## full recall plot ---------------------------------------------------------------
# Define common plot parameters
common_params <- list( 
  data  = gost_LUX, #%>% dplyr::filter(grepl("CORUM", evidence_codes)),  
  min_recall = 0.7,    max_p_value = 0.05,    grouping = "comparison",    term_column = "term_name",    distance_column = "recall",    distance_method = "euclidean", # distance parameters
  x_var = "comparison",    y_var = "term_name",    size_var = "recall",    fill_var = "p_value",    
  title_var = "g:GOSt Recall > 0.7" )
# Loop through group filters and create plots
go_LUX_result_recall <- lapply(group_filters, function(filter) {
  # Call the plotting function, which now returns a list: list(plot, dataframe)
  result_list <- do.call(plot_dual_distance_bubble, c(common_params, list(group_filter = filter)))
  plot <- result_list[[1]]  # Extract plot object
  dataframe <- result_list[]  # Extract dataframe (optional: use/store as needed)
  # Adjust plot theme if filter is NULL
  if (is.null(filter)) {
    plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  print(plot)  # display plot
  return(result_list) # store both plot and dataframe
})
go_LUX_result_recall[[1]][[1]] # plot for group_filters[1]
go_LUX_recall_table        <- go_LUX_result_recall[[1]][[2]] # <<<<<<< result table for group_filters[1] list element 2
go_LUX_term_proteins_maxID <- go_LUX_result_recall[[1]][[2]] %>%
  dplyr::select(term_name, term_size, recall, evidence_codes, intersection_gene) %>% 
  group_by(term_name) %>%
  slice_max(recall, n = 1) %>% # keep top 1 independent (line above) of term 
  ungroup() %>%
  distinct() %>%
  arrange(desc(row_number()))  # inverts order of df so easy to compare output table with plot

write.csv(go_LUX_term_proteins_maxID, "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/enriched_go_protein_maxID.csv"   , row.names = FALSE)

# export plot of interest 
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.x_go_recall.png",
  plot   = go_LUX_result_recall[[1]][[1]],
  width  = 38,  # document is 16 cm wide             (before 12cm used)
  height = 21,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# ## exclude recall = 1 go plot --------------------------------------------
# LUX_go_top10other <- gost_LUX %>%  # extract top 10 recall values with recall < 1
#   dplyr::select(term_name, recall) %>%
#   distinct() %>%
#   filter(recall < 0.8) %>%            # Only terms with recall < ...
#   arrange(desc(recall)) %>%
#   slice_head(n = 10) %>%
#   pull(term_name)
# 
# common_params <- list( 
#   data  = gost_LUX %>% dplyr::filter(term_name %in% LUX_go_top10other),                # Select top 10 rows,   # <<<
#   min_recall = 0,                                   # <<<
#   max_p_value = 0.05,    grouping = "comparison",    term_column = "term_name",    distance_column = "recall",    distance_method = "euclidean", # distance parameters
#   x_var = "comparison",    y_var = "term_name",    size_var = "recall",    fill_var = "p_value",    title_var = "g:GOSt - Recall 1" )
# # Loop through group filters and create plots
# go_LUX_result_go <- lapply(group_filters, function(filter) {
#   # Call the plotting function, which now returns a list: list(plot, dataframe)
#   result_list <- do.call(plot_dual_distance_bubble, c(common_params, list(group_filter = filter)))
#   plot <- result_list[[1]]  # Extract plot object
#   dataframe <- result_list[]  # Extract dataframe (optional: use/store as needed)
#   # Adjust plot theme if filter is NULL
#   if (is.null(filter)) {
#     plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   }
#   print(plot)  # display plot
#   return(result_list) # store both plot and dataframe
# })
# go_LUX_result_go[[1]][[1]] # plot for group_filters[1]
# go_LUX_go_table <- go_LUX_result_go[[2]][[2]] # <<<<<<< result table for group_filters[1]
# # export plot of interest 
# ggsave(
#   filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.x_go_top10other.png",
#   plot   = go_LUX_result_go[[1]][[1]],
#   width  = 32,  # document is 16 cm wide             (before 12cm used)
#   height = 24,  # 4/3 width/high ratio is common      (before 8 cm used)
#   units  = "cm",
#   dpi    = 300    # default for good quality
# )  

# ## top n p value plot ---------------------------------------------------------------
# top_n_pvalue_terms_go = 10
# # Define common plot parameters
# common_params <- list( 
#   data        = gost_LUX,      min_recall = 0,    # critial recall 0 !!!
#   # dynamically define p_value cutoff to obtain top n catgories (global minimum values count --> change of comparisons included affects displayed terms)
#   max_p_value = gost_LUX %>% arrange(p_value) %>% distinct(term_name, .keep_all = TRUE) %>% slice_head(n = top_n_pvalue_terms_go) %>% pull(p_value) %>% max(),    
#   grouping = "comparison",    term_column = "term_name",    distance_column = "recall",    distance_method = "euclidean", # distance parameters
#   x_var = "comparison",    y_var = "term_name",    size_var = "recall",    fill_var = "p_value",    title_var = "g:GOSt - pValue filtered" )
# # Loop through group filters and create plots
# go_LUX_result_p <- lapply(group_filters, function(filter) {
#   plot <- do.call(plot_dual_distance_bubble, c(common_params, list(group_filter = filter)))
#   if (is.null(filter)) {     plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))   }
#   print(plot)  # display plot
#   return(plot) # store plot 
# })
# go_LUX_result_p[[1]]
# go_LUX_p_table <- go_LUX_result_p[[2]] # <<<<<<< result table
# # export plot of interest 
# ggsave(
#   filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.x_go_pvalue.png",
#   plot   = go_LUX_result_p[[1]],
#   width  = 32,  # document is 16 cm wide             (before 12cm used)
#   height = 12,  # 4/3 width/high ratio is common      (before 8 cm used)
#   units  = "cm",
#   dpi    = 300    # default for good quality
# )  # end of gost baset enrichemtn analysis __________________________________________________________________________________________________________________


## protein family enrichment =====================================================================================================================================================
# qid list result #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#- #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
protein_entry_vector <- read.csv("/Users/mgesell/Downloads/targets.csv", header = TRUE) %>% dplyr::filter(entry != "" & CPI.Target == "Yes"
                                                                                                          ) %>% pull(entry)
enriched_go <- go_gost(query_list    = protein_entry_vector,
                       set = "Protein Families: LUX_targets",   max_term_size = 100,
                       # # NULL = all GO branches;       a vector of data sources to use. Currently, these include GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP. Please see the g:GOSt web tool for the comprehensive list and details on incorporated data sources.
                       go_source     = c("GO:BP", "GO:MF", "GO:CC", # basic
                                         "CORUM",                   # complexes
                                         "KEGG","REAC", "WP",       # signaling 
                                         "HPA", "HP", "HPA" #,      # phenotype (human phenotype atlas) 
                                         # "TF", "MIRNA"            # regulatory DNA motifs
                       )) 
plot_dual_distance_bubble(data = enriched_go ,  min_recall = 0.7,    max_p_value = 0.05,   
  group_filter = unique(enriched_go$comparison),
  grouping = "comparison",   term_column = "term_name",    distance_column = "recall",    
  distance_method = "euclidean", # distance parameters
  x_var = "comparison"  ,    y_var = "term_name",    size_var = "recall",    
  fill_var = "p_value"    ,    title_var = "g:GOSt" # plot parameters
)


(enriched_go)  c(common_params, list(group_filter = filter))
enriched_protein_families <- pF_pD_enrichment_analysis(data_subset_name     = "Protein Families: LUX_targets",
                                                       protein_entry_vector = protein_entry_vector,
                                                       mode = "protein_family")
enriched_protein_domains <- pF_pD_enrichment_analysis(data_subset_name      = "Protein Domains: LUX_targets",
                                                       protein_entry_vector = protein_entry_vector,
                                                       mode = "domain_ft")
enriched_protein_families[[2]]
enriched_protein_domains[[1]]
 #___________________________________________________________________________________________________________________________

group_filters <- list(c("panT", "CD4+" , "CD8+" , "Naive",  "Memory", "Naive CD4+", "Memory CD4+", "Naive CD8+" , "Memory CD8+" ))
protein_fam_meta_entrichment_table  <- NULL # empty df
protein_fam_meta_proteins           <- NULL # empty df
for (i in 1:length(group_filters[[1]])) {
  protein_entry_vector <- v31_LUX_data_prot_diff_abundance_sigup_upsp %>% 
    dplyr::filter(plot_heading == group_filters[[1]][i]) %>% # comment this line out to get full dataset view
    pull(entry) %>% unique()
  enriched_protein_families <- pF_pD_enrichment_analysis(data_subset_name     = group_filters[[1]][i] ,
                                                         protein_entry_vector = protein_entry_vector,
                                                         mode = "protein_family")
  print(enriched_protein_families[[1]])   # volcano plot
  
  protein_fam_meta_entrichment_table <- rbind(protein_fam_meta_entrichment_table, enriched_protein_families[[2]])# term table
  protein_fam_meta_proteins          <- rbind(protein_fam_meta_proteins         , enriched_protein_families[[3]]) # term proteins
}
# isolate per term all proteins that were detected
pf_overall_proteins_per_term <- protein_fam_meta_proteins %>%
  group_by(term) %>%
  summarise(overall_proteins_per_term = n_distinct(entry)) %>%
  ungroup()
pf_proteins_per_term <- protein_fam_meta_proteins %>%
  dplyr::select(entry, term) %>%
  distinct() %>%
  left_join(dplyr::select(proteome_upsp_202501, entry, gene_names_primary), by = "entry") %>%
  dplyr::select(term, gene_names_primary) %>%
  group_by(term) %>%
  summarise(gene_name = paste(unique(gene_names_primary), collapse = ", "))
# merge
pf_meta_genes_maxID <- pf_overall_proteins_per_term %>%
  left_join(pf_proteins_per_term, by = "term") 

protein_fam_meta_entrichment_table <- protein_fam_meta_entrichment_table %>% arrange(desc(recall))
write.csv(protein_fam_meta_entrichment_table, "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/enriched_pf_table.csv"   , row.names = FALSE)
write.csv(protein_fam_meta_proteins         , "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/enriched_pf_proteins.csv", row.names = FALSE)
write.csv(pf_meta_genes_maxID               , "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/enriched_pf_genes.csv", row.names = FALSE)
#    protein_fam_meta_proteins  # <<<<<<<<<<<
#    protein_fam_meta_proteins %>% dplyr::filter(term == "ITAM") %>% pull(entry_name) %>% unique()
#
result_pf_enrich <- plot_dual_distance_bubble(data = protein_fam_meta_entrichment_table,  min_recall = 0,    max_p_value = 0.05,   
                          group_filter = group_filters[[1]],
                          grouping = "plot_heading",   term_column = "term",    distance_column = "recall",    
                          distance_method = "euclidean", # distance parameters
                          x_var = "plot_heading"  ,    y_var = "term",    size_var = "recall",    
                          fill_var = "p_value"    ,    title_var = "Protein familiy enrichment" # plot parameters
)
result_pf_enrich[[1]]
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.x_pf.png",
  plot   = result_pf_enrich[[1]],
  width  = 32,  # document is 16 cm wide             (before 12cm used)
  height = 10,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  # end of protein familiy enrichment analysis __________________________________________________________________________________________________________________

# protein domain enrichment ..................................................................................................
group_filters <- list(c("panT", "CD4+" , "CD8+" , "Naive",  "Memory", "Naive CD4+", "Memory CD4+", "Naive CD8+" , "Memory CD8+" ))
protein_domain_meta_entrichment_table  <- NULL # empty df
protein_domain_meta_proteins           <- NULL # empty df
for (i in 1:length(group_filters[[1]])) {
  protein_entry_vector <- v31_LUX_data_prot_diff_abundance_sigup_upsp %>% 
    dplyr::filter(plot_heading == group_filters[[1]][i]) %>% # comment this line out to get full dataset view
    pull(entry) %>% unique()
  enriched_protein_domains <- pF_pD_enrichment_analysis(data_subset_name     = group_filters[[1]][i],
                                                        protein_entry_vector = protein_entry_vector,
                                                        mode = "domain_ft")
  print(enriched_protein_domains[[1]])   # volcano plot
  
  protein_domain_meta_entrichment_table <- rbind(protein_domain_meta_entrichment_table, enriched_protein_domains[[2]]) # term table
  protein_domain_meta_proteins          <- rbind(protein_domain_meta_proteins         , enriched_protein_domains[[3]]) # term proteins
}
# isolate per term all proteins that were detected
domain_overall_proteins_per_term <- protein_domain_meta_proteins %>%
  group_by(term) %>%
  summarise(overall_proteins_per_term = n_distinct(entry)) %>%
  ungroup()
domain_count_per_condition <- protein_domain_meta_proteins %>%
  group_by(plot_heading, term) %>%
  summarise(domain_count_per_condition = n()) %>%
  group_by(term) %>% 
  slice_max(domain_count_per_condition, n = 1) %>% # keep top 1 independent (line above) of term 
  ungroup() %>%
  dplyr::select(-plot_heading) %>%
  distinct()
domain_proteins_per_term <- protein_domain_meta_proteins %>%
  dplyr::select(entry, term) %>%
  distinct() %>%
  left_join(dplyr::select(proteome_upsp_202501, entry, gene_names_primary), by = "entry") %>%
  dplyr::select(term, gene_names_primary) %>%
  group_by(term) %>%
  summarise(gene_name = paste(unique(gene_names_primary), collapse = ", "))
# merge
protein_domain_meta_genes_maxID <- domain_overall_proteins_per_term %>%
  left_join(domain_count_per_condition, by = "term") %>%
  left_join(domain_proteins_per_term  , by = "term") 
##  
write.csv(protein_domain_meta_entrichment_table, "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/enriched_pd_table.csv"   , row.names = FALSE)
write.csv(protein_domain_meta_proteins         , "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/enriched_pd_proteins.csv", row.names = FALSE)
write.csv(protein_domain_meta_genes_maxID      , "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/enriched_pd_genes.csv", row.names = FALSE)

protein_domain_meta_entrichment_table <- protein_domain_meta_entrichment_table %>% arrange(desc(recall))
# protein_domain_meta_proteins  # <<<<<<<<<<<
protein_domain_meta_proteins %>% dplyr::filter(term == "Ig-like V-type") %>% pull(entry_name) %>% unique()
#
result_pd_enrich <- plot_dual_distance_bubble(data = protein_domain_meta_entrichment_table,  min_recall = 0,    max_p_value = 0.05,   
                          group_filter = group_filters[[1]],
                          grouping = "plot_heading",   term_column = "term",    distance_column = "recall",    
                          distance_method = "euclidean", # distance parameters
                          x_var = "plot_heading"  ,    y_var = "term",    size_var = "recall",    
                          fill_var = "p_value"    ,    title_var = "Protein domain enrichment" # plot parameters
)
result_pd_enrich[[1]]
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.x_pd.png",
  plot   = result_pd_enrich[[1]],
  width  = 32,  # document is 16 cm wide             (before 12cm used)
  height = 12,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good qualityrotein_fam_meta_entrichment_table <- protein_domain_meta_entrichment_table %>% arrange(desc(recall))
)
# end of protein domain enrichment analysis ______________________________________________________________________________________________________________________________________________



















# *************************************************************************************************************************************************
# *************************************************************************************************************************************************
# central processing step: IDs are transfered from high confidence tests (12 vs 12; 6 vs 6 replicates) to lower confidence test sets (3 vs 3)
# that means eventhough e.g. not sig-up in naiveCD4+ subsets any ID from panT meta subsets transfered to that subsets (ID expansion following confidence gradient)
# based on enhanced statistical power panT t-test (12 replicated and at least 9/12 replicate IDs required for imputation) is more credible than 3vs3 replicate tests from individual subsets
# --> use panT sig-enriched proteins as core abTCR community --> expand (copy) IDs into all subsets. same for other meta IDs. distributed to downstream of 
# *************************************************************************************************************************************************
# *************************************************************************************************************************************************

# # extend dataset subsets - unique(CD8 is panT + CD8)
# df_list <- list(
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "Naive" , "CD4+", "Naive CD4+" ))   %>%    mutate(plot_heading_ssEXT = "Naive CD4+"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "Naive" , "CD8+", "Naive CD8+" ))   %>%    mutate(plot_heading_ssEXT = "Naive CD8+"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "Memory", "CD4+", "Memory CD4+"))   %>%    mutate(plot_heading_ssEXT = "Memory CD4+"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "Memory", "CD8+", "Memory CD8+"))   %>%    mutate(plot_heading_ssEXT = "Memory CD8+"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "Naive" , "Naive CD4+" , "Naive CD8+" ))   %>%    mutate(plot_heading_ssEXT = "Naive"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "Memory", "Memory CD4+", "Memory CD8+"))   %>%    mutate(plot_heading_ssEXT = "Memory"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "CD4+"  , "Naive CD4+" , "Memory CD4+"))   %>%    mutate(plot_heading_ssEXT = "CD4+"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT", "CD8+"  , "Naive CD8+" , "Memory CD8+"))   %>%    mutate(plot_heading_ssEXT = "CD8+"),
#   v31_LUX_data_prot_diff_abundance_sigup_upsp %>%    filter(plot_heading %in% c("panT"                                 ))   %>%    mutate(plot_heading_ssEXT = "panT")
# )
# v31_LUX_data_prot_diff_abundance_sigup_ssEXTENDED <- bind_rows(df_list) %>%  distinct() %>%
#   left_join(proteome_20250804 %>% dplyr::select(entry, entry_name, involvement_in_disease, domain_cc, domain_ft, topological_domain, protein_families),  by = "entry") 
# 


###############################################################################################################################################
# CYTOSCAPE -----------------------------------------------------------------------------------------------------------------------------------
# cytoscape information mapping, sig-up proteins vs. annotations
## regenerate table when changing GO list
sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/thesis_figures_functions", pattern = "\\.R$", full.names = TRUE), source)
GO_descendants_df <- go_children_table_generator2(my_terms = c("GO:0036398", "GO:0036399", "GO:0050852", "GO:0050862", "GO:0050860", "GO:0042110", "GO:0050870", "GO:0050868",
                                         "GO:0043235", "GO:0043113", "GO:0031623", "GO:0001881", "GO:0043112", "GO:0038023", "GO:0003824", "GO:0030674",
                                         "GO:0000151", "GO:0061630", "GO:0005215", "GO:0008233", "GO:0042105", "GO:0005085", "GO:0005096", "GO:0003924",
                                         "GO:0007165", "GO:0004672", "GO:0004721", "GO:0004930", "GO:0003823", "GO:0042611", "GO:0019838", "GO:0004896",
                                         "GO:0006955", "GO:0007155", "GO:0005178", "GO:0005215", "GO:0005216", "GO:0016787", "GO:0016491", "GO:0006915",
                                         "GO:0045454", "GO:0031295", "GO:0140359", "GO:0023052", "GO:0035591", "GO:0097197", "GO:0002684", "GO:0002683", 
                                         "GO:0007166"
                                         ), only_new = TRUE) # define whether only new terms should be searched and appended to df. otherwise full serach - which takes long
## load resources ...........
GO_descendants_df     <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/GO_term_children/GO_term_children.csv", header = TRUE , check.names = FALSE)

protein_family_df     <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/protein_families/human_upsp_proteinfamilies_2025_07_23.csv", header = TRUE) %>% 
  dplyr::rename("protein_families" = "Protein.families", "entry" = "Entry")

string_uniprot_mapping <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/string_uniprot_ENSP_mapping/stringTOupsp_human_20250801.csv", header = TRUE) %>%
  dplyr::select(Entry, STRING) %>% dplyr::rename("entry" = "Entry", "ENSPid_string" = "STRING") %>%
  mutate(ENSPid_string = gsub(";", "", ENSPid_string))

proteome_domain <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/protein_domains/domains_upsp_20250801.csv", header = TRUE) %>%
  dplyr::select(Entry, Domain..CC., Domain..FT.) %>%
  dplyr::rename("entry" = "Entry", "domains_cc" = "Domain..CC.", "domains" = "Domain..FT.") %>%
  mutate(ITIM = ifelse(str_detect(domains, "ITIM") | str_detect(domains_cc, "ITIM"), "Yes", "No"),
         ITAM = ifelse(str_detect(domains, "ITAM") | str_detect(domains_cc, "ITAM"), "Yes", "No"),
         ITSM = ifelse(str_detect(domains, "ITSM") | str_detect(domains_cc, "ITSM"), "Yes", "No"),
         ITxM = ifelse(ITIM == "Yes" | ITAM == "Yes" |ITSM == "Yes", "Yes", "No")) 

poi_table       <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv", header = TRUE)

immune_diseases <- str_c(  c("autoimmunity", "immunodeficiency", "impaired immune", "B cell deficiency", "T cell deficiency", "auto-inflammatory",
                             "immune dysregulation", "rheuma", "diabetes", "lupus", "graves", "addison", "Sjögren", "psoriasis", "bowel", "arteritis", "vasculitides", "celiac", "allergic", 
                             "leukemia", "lymphopenia", "hypogammaglobulinemia", "lymphoma", "T cell disorders", "Kaposi sarcoma", "Multiple sclerosis", "Hashimoto", "Myasthenia", "Pernicious", "Antiphospholipid",
                             "Dermatomyositis", "Scleroderma", "cholangitis", "Wiskott-Aldrich", "DiGeorge", "granulomatous", "IgA deficiency", "Granulomatosis", "polyangiitis", "Churg-Strauss", "arteritis", "HIV", "Sarcoidosis", "polyglandular"
                              ), 
                         collapse = "|")
proteome_immune_disease <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/protein_disease/upsp_disease_2025-08-04.csv", header = TRUE) %>%
  dplyr::rename("entry" = "Entry", "entry_name" = "Entry.Name", "disease" = "Involvement.in.disease") %>%
  mutate(immune_disease = case_when(str_detect(disease, regex(immune_diseases, ignore_case = TRUE)) ~ "Yes", TRUE ~"No"))

proteome_drug_20250804 <- proteome_20250804 %>%
  dplyr::select(entry,  drugbank, drugcentral) %>%
  mutate(drug_target    = case_when(drugbank != "" |  drugcentral != "" ~ "Yes", TRUE ~"No"))
table(proteome_drug_20250804$drug_target)

APMS_TCR_ppis <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/publications/TCR_co-ip_literature.csv", header = TRUE) 
# [1] "APMS_Romain2022_CD3z_FC2"           "APMS_Romain2022_CD3z_FC10"          "APMS_Romain2022_ZAP70_FC2"          "APMS_Romain2022_ZAP70_FC10"        
# [5] "APMS_Romain2019_Lck"                "combined"                           "x"                                  "APMS_Romain2022_CD3z_anytime_FC10" 
# [9] "APMS_Romain2022_CD3z_anytime_FC2"   "APMS_Romain2022_ZAP70_anytime_FC10" "APMS_Romain2022_ZAP70_anytime_FC2"  "combined_anytime_and_Lck_2019"   

#  string 
library(rbioapi)             # install.packages("rbioapi") # required for string interactor mapping with confidence score x and interaction criterium physica/funcitonal(full)

# string_targets_chapter3  <- c("CD3D_HUMAN",           "CD3E_HUMAN",           "CD3G_HUMAN",           "CD3Z_HUMAN",            "CD4_HUMAN",           "CD8A_HUMAN" ,          "CD8B_HUMAN")  # define target & or fixed complex partners as well as condition
string_ids          <-        c("9606.ENSP00000300692", "9606.ENSP00000354566", "9606.ENSP00000431445", "9606.ENSP00000354782", "9606.ENSP00000011653", "9606.ENSP00000386559", "9606.ENSP00000331172") # manually retrieved. CD8B directly from string beause not uniprot annotated string id

string_AB_physical_400 <- rba_string_interaction_partners(
  ids            = string_ids,
  species        = 9606,
  required_score = 700,         # set to 400 for medium, 700 for high confidence
  network_type   = "physical" ) %>%    # for direct/physical; use "functional" for total
  # map ENSP (.$stringId_B) and gene_names_primary (.$preferredName_B) back to uniprot format
  left_join(string_uniprot_mapping, by = c("stringId_B"      = "ENSPid_string")) %>% # double mapping part 1
  dplyr::rename("entry_ENS" = "entry") %>%
  left_join(proteome_upsp_202501 %>% dplyr::select(entry, gene_names_primary)  , by = c("preferredName_B" = "gene_names_primary")) %>% # double mapping part 2
  dplyr::rename("entry_primgenename" = "entry") %>%
  dplyr::mutate(entry = coalesce(entry_ENS, entry_primgenename), 
                match = entry_ENS == entry_primgenename) %>% # QC column showing issues of conflicting mappings (so for none found)
  pull(entry) %>% unique()

string_AB_full_400 <- rba_string_interaction_partners(
  ids            = string_ids,
  species        = 9606,
  required_score = 400,         # set to 400 for medium, 700 for high confidence
  network_type   = "functional" ) %>%    # for direct/physical; use "functional" for total
  # map ENSP (.$stringId_B) and gene_names_primary (.$preferredName_B) back to uniprot format
  left_join(string_uniprot_mapping, by = c("stringId_B"      = "ENSPid_string")) %>% # double mapping part 1
  dplyr::rename("entry_ENS" = "entry") %>%
  left_join(proteome_upsp_202501 %>% dplyr::select(entry, gene_names_primary)  , by = c("preferredName_B" = "gene_names_primary")) %>% # double mapping part 2
  dplyr::rename("entry_primgenename" = "entry") %>%
  dplyr::mutate(entry = coalesce(entry_ENS, entry_primgenename), 
                match = entry_ENS == entry_primgenename) %>% # QC column showing issues of conflicting mappings (so for none found)
  pull(entry) %>% unique()

string_AB_full_700 <- rba_string_interaction_partners(
  ids            = string_ids,
  species        = 9606,
  required_score = 700,         # set to 400 for medium, 700 for high confidence
  network_type   = "functional" ) %>%    # for direct/physical; use "functional" for total
  # map ENSP (.$stringId_B) and gene_names_primary (.$preferredName_B) back to uniprot format
    left_join(string_uniprot_mapping, by = c("stringId_B"      = "ENSPid_string")) %>% # double mapping part 1
  dplyr::rename("entry_ENS" = "entry") %>%
  left_join(proteome_upsp_202501 %>% dplyr::select(entry, gene_names_primary)  , by = c("preferredName_B" = "gene_names_primary")) %>% # double mapping part 2
  dplyr::rename("entry_primgenename" = "entry") %>%
  dplyr::mutate(entry = coalesce(entry_ENS, entry_primgenename), 
                match = entry_ENS == entry_primgenename) %>% # QC column showing issues of conflicting mappings (so for none found)
  pull(entry) %>% unique()


#######################
# helper vecotrs to define unique sets 
meta_all_but_naive_list    <-  v31_LUX_data_prot_diff_abundance_sigup %>% filter(plot_heading %in% c("Meta panT",               "Meta Memory", "Meta CD4+", "Meta CD8+")) %>% pull(entry) %>% unique()
meta_all_but_memory_list   <-  v31_LUX_data_prot_diff_abundance_sigup %>% filter(plot_heading %in% c("Meta panT", "Meta Naive",                "Meta CD4+", "Meta CD8+")) %>% pull(entry) %>% unique()
meta_all_but_CD4_list      <-  v31_LUX_data_prot_diff_abundance_sigup %>% filter(plot_heading %in% c("Meta panT", "Meta Naive", "Meta Memory",              "Meta CD8+")) %>% pull(entry) %>% unique()
meta_all_but_CD8_list      <-  v31_LUX_data_prot_diff_abundance_sigup %>% filter(plot_heading %in% c("Meta panT", "Meta Naive", "Meta Memory", "Meta CD4+"             )) %>% pull(entry) %>% unique()
#
meta_all_but_nCD4_list   <-  v31_LUX_data_prot_diff_abundance_sigup   %>% filter(plot_heading %in% c("Meta panT", "Meta Naive", "Meta Memory", "Meta CD4+", "Meta CD8+",               "Memory CD4+", "Naive CD8+",  "Memory CD8+")) %>% pull(entry) %>% unique()
meta_all_but_nnCD4_list  <-  v31_LUX_data_prot_diff_abundance_sigup   %>% filter(plot_heading %in% c("Meta panT", "Meta Naive", "Meta Memory", "Meta CD4+", "Meta CD8+", "Naive CD4+",                "Naive CD8+",  "Memory CD8+")) %>% pull(entry) %>% unique()
meta_all_but_nCD8_list   <-  v31_LUX_data_prot_diff_abundance_sigup   %>% filter(plot_heading %in% c("Meta panT", "Meta Naive", "Meta Memory", "Meta CD4+", "Meta CD8+", "Naive CD4+", "Memory CD4+"              ,  "Memory CD8+")) %>% pull(entry) %>% unique()
meta_all_but_nnCD8_list  <-  v31_LUX_data_prot_diff_abundance_sigup   %>% filter(plot_heading %in% c("Meta panT", "Meta Naive", "Meta Memory", "Meta CD4+", "Meta CD8+", "Naive CD4+", "Memory CD4+", "Naive CD8+"                )) %>% pull(entry) %>% unique()

# network will be restricted to selected "nice to explain" categroies (weird overlaps ignored or introduced later lets see)
cytoscape_bait_network  <- v31_LUX_data_prot_diff_abundance_sigup %>%
  mutate(bait = case_when(
    # panT master core overlap list
    plot_heading == "Meta panT" ~ "panT",
    # meta subsets
    plot_heading == "Meta Naive"  & !(entry %in% meta_all_but_naive_list)  ~ "Naive",   # 
    plot_heading == "Meta Memory" & !(entry %in% meta_all_but_memory_list) ~ "Memory",
    plot_heading == "Meta CD4+"   & !(entry %in% meta_all_but_CD4_list)    ~ "CD4+",
    plot_heading == "Meta CD8+"   & !(entry %in% meta_all_but_CD8_list)    ~ "CD8+",
    # subset unique
    plot_heading == "Naive CD4+"  & !(entry %in% meta_all_but_nCD4_list )  ~ "Naive CD4+",  # exclude any proteins from unique-to-subset assignment that are also contained in any of meta categories
    plot_heading == "Memory CD4+" & !(entry %in% meta_all_but_nnCD4_list)  ~ "Memory CD4+",
    plot_heading == "Naive CD8+"  & !(entry %in% meta_all_but_nCD8_list )  ~ "Naive CD8+",
    plot_heading == "Memory CD8+" & !(entry %in% meta_all_but_nnCD8_list)  ~ "Memory CD8+",
    TRUE ~ NA_character_
  )) %>% 
  filter(!is.na(bait)) %>%   # eliminate any unassigned proteins
  dplyr::select(entry, entry_name, bait) %>%
  distinct() %>%
  left_join(proteome_upsp_202501 %>% 
              dplyr::select(entry, gene_names_primary, protein_names), 
            by = "entry") %>%
  dplyr::select(bait, gene_names_primary, entry, entry_name, protein_names) %>%
  dplyr::rename("prey" = "gene_names_primary")
  
## mapping resources -------------------------
cytoscape_info_annotation <- v31_LUX_data_prot_diff_abundance_sigup %>% 
  dplyr::select(entry, entry_name, overlap, proximity_id_cout) %>%
  dplyr::filter(entry_name != "LYSC_LYSEN")  %>%
  distinct() %>%
  left_join(cytoscape_bait_network %>% dplyr::select(entry, bait)               , by = "entry") %>% # annotate subset so when removing subsed edges there is still way to select subset nodes
  # account for high confidence (4x sig-up) proteins that don´t fit in "nice to explain categories" (4 subsets + 4 meta and 1 panT subset excludes overlaps e.g. naiveCD8-naiveCD4-metaCD8)
  mutate(bait = ifelse(is.na(bait) & proximity_id_cout >= 4 |            # either ... supported by at least 4 enrichment analysis 
                       is.na(bait) & str_detect(overlap, "Meta|panT") & proximity_id_cout >= 2,  # ... or supported at least by 1 meta subset an 1 other subset
                       "panT_extended", bait))   %>%
  left_join(protein_family_df      %>% dplyr::select(entry, protein_families)   , by = "entry") %>% # annotate subset so when removing subsed edges there is still way to select subset nodes
  left_join(proteome_upsp_202501   %>% dplyr::select(entry, protein_names, gene_names_primary, function_cc, gene_ontology_molecular_function, gene_ontology_i_ds), by = "entry")  %>%
  left_join(string_uniprot_mapping, by = "entry") %>% # allows to add nodes by string id into network 
  left_join(proteome_immune_disease %>% dplyr::select(entry, immune_disease), by = "entry") %>% # allows to add nodes by string id into network 
  left_join(proteome_drug_20250804  %>% dplyr::select(entry, drug_target)   , by = "entry") # allows to add nodes by string id into network 

  
## mapping cytoscape values ----------------
cytoscape_info_annotation_go <- cytoscape_info_annotation  %>%
              mutate(
                abTCR_complex            = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`alpha-beta T cell receptor complex` %>% na.omit())), "Yes", "No"),
                tcr_signalosome          = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`TCR signalosome` %>% na.omit())), "Yes", "No"),
                tcr_signaloseome_assembly= if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`TCR signalosome assembly` %>% na.omit())), "Yes", "No"),
                
                tcr_signaling_pathway    = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`T cell receptor signaling pathway` %>% na.omit())), "Yes", "No"),
                tcr_signaling_plus       = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`positive regulation of T cell receptor signaling pathway` %>% na.omit())), "Yes", "No"),
                tcr_signaling_minus      = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`negative regulation of T cell receptor signaling pathway` %>% na.omit())), "Yes", "No"),
                t_act                    = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`T cell activation` %>% na.omit())), "Yes", "No"),
                t_act_plus               = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`positive regulation of T cell activation` %>% na.omit())), "Yes", "No"),
                t_act_minus              = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`negative regulation of T cell activation` %>% na.omit())), "Yes", "No"),
                t_co_stim                = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`T cell costimulation` %>% na.omit())), "Yes", "No"),
                ISP_reg_plus             = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`positive regulation of immune system process` %>% na.omit())), "Yes", "No"),
                ISP_reg_minus            = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`negative regulation of immune system process` %>% na.omit())), "Yes", "No"),
                
                surfaceome_receptor_sign = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`cell surface receptor signaling pathway` %>% na.omit())), "Yes", "No"),
                receptor_complex         = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`receptor complex` %>% na.omit())), "Yes", "No"),
                receptor_clustering      = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`receptor clustering` %>% na.omit())), "Yes", "No"),
                receptor_internalization = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`receptor internalization` %>% na.omit())), "Yes", "No"),
                receptor_recycling       = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`receptor recycling` %>% na.omit())), "Yes", "No"),
                receptor_metabolism      = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`receptor metabolic process` %>% na.omit())), "Yes", "No"),
                signaling_receptor       = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`signaling receptor activity` %>% na.omit())), "Yes", "No"),
                signal_transduction      = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`signal transduction` %>% na.omit())), "Yes", "No"),
                signaling                = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`signaling` %>% na.omit())), "Yes", "No"),
                
                GEF                      = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`guanyl-nucleotide exchange factor activity` %>% na.omit())), "Yes", "No"),
                GAP                      = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`GTPase activator activity` %>% na.omit())), "Yes", "No"),
                GTPase                   = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`GTPase activity` %>% na.omit())), "Yes", "No"),
                kinase                   = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`protein kinase activity` %>% na.omit())), "Yes", "No"),
                phosphatase              = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`phosphoprotein phosphatase activity` %>% na.omit())), "Yes", "No"),
                tetraspaning_enr_microdo = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`tetraspanin-enriched microdomain` %>% na.omit())), "Yes", "No"),

                gpcr                     = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`G protein-coupled receptor activity` %>% na.omit())), "Yes", "No"),
                enzyme                   = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`catalytic activity` %>% na.omit())), "Yes", "No"),
                adapter                  = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`protein-macromolecule adaptor activity` %>% na.omit())), "Yes", "No"),
                adapter_in_signaling     = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`signaling adaptor activity` %>% na.omit())), "Yes", "No"),
                ubiquitin_ligase_complex = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`ubiquitin ligase complex` %>% na.omit())), "Yes", "No"),
                ubiquitin_ligase         = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`ubiquitin protein ligase activity` %>% na.omit())), "Yes", "No"),
                transporter              = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`transporter activity` %>% na.omit())), "Yes", "No"),
                transporter_activity     = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`transporter activity` %>% na.omit())), "Yes", "No"),
                ion_channel              = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`monoatomic ion channel activity` %>% na.omit())), "Yes", "No"),
                abc_transporter          = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`ABC-type transporter activity` %>% na.omit())), "Yes", "No"),
                peptidase                = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`peptidase activity` %>% na.omit())), "Yes", "No"),
                cytokine_receptor        = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`cytokine receptor activity` %>% na.omit())), "Yes", "No"),           
                
                antigen_binding          = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`antigen binding` %>% na.omit())), "Yes", "No"),
                mhc_complex              = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`MHC protein complex` %>% na.omit())), "Yes", "No"),                
                growth_factor_binding    = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`growth factor binding` %>% na.omit())), "Yes", "No"),
                immune_response          = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`immune response` %>% na.omit())), "Yes", "No"),
                cell_adhesion            = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`cell adhesion` %>% na.omit())), "Yes", "No"),
                integrin_binding         = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`integrin binding` %>% na.omit())), "Yes", "No"),
                hydrolase_activity       = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`hydrolase activity` %>% na.omit())), "Yes", "No"),
                oxidoreductase_activity  = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`oxidoreductase activity` %>% na.omit())), "Yes", "No"),
                apoptotic_process        = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`apoptotic process` %>% na.omit())), "Yes", "No"),
                redox_homeostasis        = if_else(str_split(gene_ontology_i_ds, ";") %>% map_lgl(~ any(str_trim(.)    %in%   GO_descendants_df$`cell redox homeostasis` %>% na.omit())), "Yes", "No")
              ) %>% 
              dplyr::select(-gene_ontology_i_ds) 

cytoscape_info_annotation_full <- cytoscape_info_annotation_go %>%
  mutate(                
    # protein families / protein name based annotations
    ADP_ribosyl_cyclase      = ifelse(grepl("ADP-ribosyl cyclase family"        , protein_families), "Yes", "No"),       # 2x (of total 2)  --> 100% !!!
    recep_of_compl_actRCQ    = ifelse(grepl("Receptors of complement activation", protein_families), "Yes", "No"),       # 3x (of total 4)  --> 75%
    CD99_family              = ifelse(grepl("CD99 family"                       , protein_families), "Yes", "No"),       # 2x (of total 3)  --> 66%
    SLAM_famliy              = ifelse(grepl("SLAM family member"                , protein_names   ), "Yes", "No"),       # 6x (of total 9)  --> 66%
    ICAM_family              = ifelse(grepl("Immunoglobulin superfamily, ICAM family", protein_families), "Yes", "No"),  # 3x (of total 5)  --> 60 %
    ITA_family               = ifelse(grepl("Integrin alpha chain family"       , protein_families), "Yes", "No"),       # 9x (of total 18) --> 50%
    MHC_1_famliy             = ifelse(grepl("MHC class I family"                , protein_families), "Yes", "No"),       # 5x (of total 16) --> 33% 
    tetraspanin              = ifelse(grepl("Tetraspanin|tetraspanin"           , protein_families), "Yes", "No"),       # 4x (of total 37) 
    calomdulin_family        = ifelse(grepl("Calmodulin family"                 , protein_families), "Yes", "No"),       # 1x (of total 6)
    # protein domains
    ITIM         = ifelse(entry %in% c(proteome_domain %>% dplyr::filter(ITIM == "Yes") %>% pull(entry) %>% unique()), "Yes", "No"),
    ITAM         = ifelse(entry %in% c(proteome_domain %>% dplyr::filter(ITAM == "Yes") %>% pull(entry) %>% unique()), "Yes", "No"),
    ITSM         = ifelse(entry %in% c(proteome_domain %>% dplyr::filter(ITSM == "Yes") %>% pull(entry) %>% unique()), "Yes", "No"),
    ITxM         = ifelse(entry %in% c(proteome_domain %>% dplyr::filter(ITxM == "Yes") %>% pull(entry) %>% unique()), "Yes", "No"),
    # APMS TCR ppis
    APMS_CD3z_FC2         = ifelse(tolower(gene_names_primary) %in% tolower(APMS_TCR_ppis$CD3z_NS_Romain2022_FC2)   , "Yes", "No"),
    APMS_CD3e_FC2         = ifelse(tolower(gene_names_primary) %in% tolower(APMS_TCR_ppis$CD3e_NS_Romain2019_FC2)   , "Yes", "No"),
    APMS_ZAP70_FC2        = ifelse(tolower(gene_names_primary) %in% tolower(APMS_TCR_ppis$ZAP70_NS_Romain2022_FC2)  , "Yes", "No"),
    APMS_Lck              = ifelse(tolower(gene_names_primary) %in% tolower(APMS_TCR_ppis$Lck_NS_Romain2019_FC2)    , "Yes", "No"),
    APMS_FC2_overlap_count= rowSums( across( c(APMS_CD3z_FC2, APMS_CD3e_FC2, APMS_ZAP70_FC2, APMS_Lck),~ . == "Yes")),
    APMS_combined_NS_CD3  = ifelse(#APMS_FC2_overlap_count >= 2 & # for now trust 1-time hits
      tolower(gene_names_primary) %in% unique(c(tolower(APMS_TCR_ppis$CD3e_NS_Romain2019_FC2) , tolower(APMS_TCR_ppis$CD3z_NS_Romain2022_FC2)))     , "Yes", "No"),
    APMS_combined_NS      = ifelse(#APMS_FC2_overlap_count >= 2 & # for now trust 1-time hits
                                     tolower(gene_names_primary) %in% unique(c(tolower(APMS_TCR_ppis$Lck_NS_Romain2019_FC2)     , tolower(APMS_TCR_ppis$CD3e_NS_Romain2019_FC2)      , tolower(APMS_TCR_ppis$CD3z_NS_Romain2022_FC2)   , tolower(APMS_TCR_ppis$ZAP70_NS_Romain2022_FC2)))     , "Yes", "No"),
    APMS_combined_meta    = ifelse(#APMS_FC2_overlap_count >= 2 & # for now trust 1-time hits
                                   tolower(gene_names_primary) %in% unique(c(tolower(APMS_TCR_ppis$Lck_anytime_Romain2019_FC2), tolower(APMS_TCR_ppis$CD3e_allStim_Romain2019_FC2) , tolower(APMS_TCR_ppis$CD3z_allStim_Romain2022_FC2), tolower(APMS_TCR_ppis$ZAP70_allStim_Romain2022_FC2))), "Yes", "No"),
    # poi table based
    CPI_candidate     = ifelse(entry_name %in% c(poi_table %>% dplyr::filter(CPIs != "")  %>% pull(CPIs)), "Yes", "No"),
    antibodies_bought = ifelse(entry_name %in% c(poi_table %>% dplyr::filter(antibodies_bought != "")  %>% pull(antibodies_bought)), "Yes", "No"),
    # GO based annotations
    abTCR_and_CD3_chains  = ifelse(entry_name %in% abTCR_chains      , "Yes", "No"),
         meta_surfaceome       = ifelse(entry      %in% surface_annotations$cspa_2015surfy_2018tcsa_2021cd_antigen_veneer_proteome_high, "Yes", "No"),
           # string interactor types ~ outline width
    string_AB_physical_400   = ifelse(entry %in% string_AB_physical_400, 10, 3),   # border width
    string_AB_full_400       = ifelse(entry %in% string_AB_full_400    , 10, 3),   # border width string_AB_full_400
    string_AB_full_700       = ifelse(entry %in% string_AB_full_700    , 10, 3),   # border width string_AB_full_400
    meta_tcr_signaling = case_when(
              tcr_signaling_plus  == "Yes" & tcr_signaling_minus == "Yes" ~ "Dual Role", # DEFINE FIRST - will pick first match so order important- account for dual roles
              tcr_signaling_plus  == "Yes"  ~ "Plus",
              tcr_signaling_minus == "Yes"  ~ "Minus",
              TRUE ~ "Other"),
    meta_Tact = case_when(
              t_act_plus  == "Yes" & t_act_minus == "Yes" ~"Dual Role", # DEFINE FIRST - will pick first match so order important- account for dual roles
              t_act_plus  == "Yes" ~ "Plus",
              t_act_minus == "Yes" ~ "Minus",
              TRUE ~ "Other"),
    meta_ISP = case_when(
              ISP_reg_plus  == "Yes" & ISP_reg_minus == "Yes" ~"Dual Role", # DEFINE FIRST - will pick first match so order important- account for dual roles
              ISP_reg_plus  == "Yes" ~ "Plus",
              ISP_reg_minus == "Yes" ~ "Minus",
              TRUE ~ "Other"),
    cytoscape_shape = case_when(
             kinase == "Yes" | phosphatase == "Yes"                               ~ "DIAMOND",  # signaling modulators
             ubiquitin_ligase == "Yes" | ubiquitin_ligase_complex   == "Yes"      ~ "PARALLELOGRAM",
             GTPase == "Yes" |  GAP == "Yes" | GEF == "Yes"                       ~ "HEXAGON",  # signaling modulators
             abTCR_and_CD3_chains == "Yes" | tcr_signaling_pathway == "Yes" | tcr_signaling_plus == "Yes" | tcr_signaling_minus == "Yes" |  meta_tcr_signaling != "Other"    ~ "ROUND_RECTANGLE", # TCR components & TCR signaling modulators
             t_co_stim  == "Yes" | t_act_plus == "Yes"  | t_act_minus == "Yes"    ~ "RECTANGLE", 
             peptidase        == "Yes"                                            ~ "OCTAGON", 
             transporter      == "Yes" | abc_transporter  == "Yes" | transporter_activity == "Yes"  | ion_channel   == "Yes"  ~ "VEE",
             adapter          == "Yes" | adapter_in_signaling == "Yes"            ~ "TRIANGLE",                     
             TRUE ~ "ELLIPSE"
             ),
    cytoscape_fill = case_when(
             abTCR_and_CD3_chains  == "Yes"         ~ "#1a9850", # top priority to highlight these (uniformly)
             # ADP_ribosyl_cyclase    recep_of_compl_actRCQ  CD99_family SLAM_famliy  ITA_family
             ADP_ribosyl_cyclase    == "Yes"                                                     ~ "#b3b3e6",
             recep_of_compl_actRCQ  == "Yes"                                                     ~ "#cc99ff",
             CD99_family            == "Yes"                                                     ~ "#ffccff",
             # SLAM_famliy          == "Yes"                                                     ~ "#b3b3b3",
             ICAM_family            == "Yes"                                                     ~ "#ff99ce",
             # ITA_family           == "Yes"                                                     ~ "#b3d9ff",
             # tetraspanin          == "Yes"                                                     ~ "#c7bee7",
             
             # meta_tcr_signaling  == "Dual Role" | GTPase == "Yes" | ITSM == "Yes" | meta_Tact == "Dual Role"  ~ "#ffe5b4",    # 
             # #
             # tcr_signaling_plus  == "Yes" | t_co_stim  == "Yes" | t_act_plus == "Yes" | tcr_signaling_pathway == "Yes" | 
             # kinase == "Yes" | GAP == "Yes"  | ITAM == "Yes"                                            ~ "#91cf60",    #  
             # #
             # tcr_signaling_minus == "Yes"| t_act_minus == "Yes" | 
             # phosphatase == "Yes" | GEF == "Yes" |  ITIM == "Yes" |
             # ubiquitin_ligase  == "Yes" | ubiquitin_ligase_complex   == "Yes"                           ~ "#ff6961",  
             
             # meta_tcr_signaling  == "Dual Role"          ~ "#fed976",
             meta_Tact == "Dual Role"                    ~ "#fed976",    # ffe5b4  # | ITSM == "Yes"
             #
             tcr_signaling_plus  == "Yes"                                                         ~ "#1a9850",    #  3d73a9   # | tcr_signaling_pathway == "Yes"  
             tcr_signaling_minus == "Yes"                                                         ~ "#d73027",
             
             t_co_stim  == "Yes" | t_act_plus == "Yes" | ISP_reg_plus == "Yes"                    ~ "#91cf60", # | ITAM == "Yes" 
             t_act_minus == "Yes" | ISP_reg_minus == "Yes"                                        ~ "#fc8d59", # | ITIM == "Yes"  
             
             #
             ITAM == "Yes"               | ITSM == "Yes"   | ITIM == "Yes"                        ~ "#6699ff", # "#e6e600" , #d2a993
             GTPase             == "Yes" | GEF == "Yes"    | GAP == "Yes"    | 
             phosphatase        == "Yes" | kinase == "Yes" |
             ubiquitin_ligase   == "Yes" | ubiquitin_ligase_complex == "Yes" |
             cytokine_receptor  == "Yes" | signaling_receptor    == "Yes"    | surfaceome_receptor_sign == "Yes"  ~ "#99ccff", # "#e6e600" , #d2a993 
             
             cytoscape_shape == "VEE"      ~ "#ffffff",  # d2a993
             cytoscape_shape == "TRIANGLE" ~ "#ffffff",  # d5eff6
             cytoscape_shape == "OCTAGON"  ~ "#ffffff",
             
             TRUE ~ "#ffffff" # "#d9d9d9"
             ),

        cytoscape_label = case_when(cytoscape_fill %in% c("#1a9850", "#d73027", "#6699ff")  ~ "#ffffff",   TRUE ~ "#000000")
        ) 
# colors
easy_colors <- c("#40bf40", "#ffa31a", "#ff471a", 
                 "#000000", "#595959",  "#bfbfbf", "#b2b266", "#d2a479", "#ff6666", "#ff9999", "#ffd480", "#F6F3a9", "#ffff00", "#66ff66", "#ccff66", "#9fdfbf", "#6699ff", "#99ccff", "#66ccff",  "#d279d2", "#cc99ff", "#c7bee7", "#ffbbe1")
# filter non-... values
keep_cols <- !sapply(cytoscape_info_annotation_full, function(x) all(is.na(x) | x == "" | x =="No" | x =="Other"))
cytoscape_info_annotation_full <- cytoscape_info_annotation_full[, keep_cols] 
#
write.csv(cytoscape_info_annotation_full, "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/cytoscape_info_annotation.csv", row.names = FALSE)
#
## qid merging of previous manual annotation with new dataframe
# xxx <- read.csv("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/shs2.27_golden-analysis_v2/shs2.22_goden_analysis_merge_file.csv", header = TRUE)
# xxxx <- cytoscape_info_annotation_full %>% dplyr::select(
#   entry, entry_name, overlap, proximity_id_cout, bait, protein_families,
#   protein_names, gene_names_primary, function_cc, gene_ontology_molecular_function, ENSPid_string,
#   immune_disease, abTCR_complex, abTCR_and_CD3_chains, meta_tcr_signaling, meta_Tact,
#   CPI_candidate, meta_surfaceome, meta_ISP, ITxM, signaling,
#   cell_adhesion, kinase, phosphatase, ubiquitin_ligase_complex, ubiquitin_ligase,
#   APMS_CD3z_FC2, APMS_CD3e_FC2, APMS_ZAP70_FC2, APMS_Lck, APMS_FC2_overlap_count,
#   APMS_combined_NS_CD3, APMS_combined_NS, APMS_combined_meta) %>%
#   left_join(xxx, by = "entry")
# write.csv(xxxx, "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/shs2.27_golden-analysis_v2/golden_analysis.csv", row.names = FALSE)

paste(c("CD3 and TCR chains in datast: ", cytoscape_info_annotation_full %>% dplyr::filter(abTCR_and_CD3_chains  == "Yes") %>% pull(entry_name) %>% unique()))
length(cytoscape_info_annotation_full %>% dplyr::filter(abTCR_and_CD3_chains  == "Yes") %>% pull(entry_name) %>% unique())


library(RCy3)
cytoscapePing()  # Confirms Cytoscape is reachable
# # cytoscape  "" fix +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # nodes <- getTableColumns('node') %>%
# #   mutate(across(c("shared name", entry, name, entry_name, protein_names), # 
# #                 ~gsub('"', '', .x)))cyto  
# # Get the node table as a dataframe & Remove quotes from the column of interest (e.g., 'names')
# nodes <- getTableColumns('node') %>%
#   mutate(across(c("shared name", entry, name, entry_name, protein_names, Matching.Attribute), # 
#                 ~gsub('"', '', .x)))
# # Push updated table back to Cytoscape (assumes names and SUIDs match)
# loadTableData(nodes, data.key.column = "SUID", table.key.column = "SUID")
# _____________________________________________________________________________________________________________________________________________________________________________________________________________


# append info annotation table to cytoscape node table (existing colums are updated)
loadTableData(
  data = cytoscape_info_annotation_full,
  data.key.column = "entry",       # Key column in your data frame (node IDs)
  table = "node",                  # Load into node table
  table.key.column = "entry"       # Corresponding key column in Cytoscape node table (usually "name")
)
# clone the netwokr
cloneNetwork()

cytoscape_info_annotation_full %>% 
  dplyr::filter(gene_names_primary == "PDCD1") %>% 
  dplyr::select(bait, gene_names_primary, meta_surfaceome, CPI_candidate, tcr_signaling_pathway, meta_tcr_signaling, meta_Tact, meta_ISP, ITxM) %>%
  pivot_longer(
    cols = c(meta_surfaceome, CPI_candidate, tcr_signaling_pathway, meta_tcr_signaling, meta_Tact, meta_ISP, ITxM),
    names_to = "annotation_type",
    values_to = "annotation_value"
  )


# manually add nodes: 
# panT_extended_nodes <- cytoscape_info_annotation_full %>% dplyr::filter(bait == "panT_extended") %>% pull(gene_names_primary) %>% unique()
# addCyNodes( c("TRAV26-1",   "TRAV6", "TRAV8-1",    "TRBV7-3")) # panT_extended_nodes)
# manually add edges
# addCyEdges(list(c('Protein_X', 'Protein_Y')), edgeType = 'interacts with')

# overwrite fill with category of interest
c("APMS_combined_zero","APMS_combined_anytime", "immune_disease", "abTCR_complex", "tcr_signalosome", "receptor_clustering",
  "receptor_internalization", "receptor_recycling", "gpcr", "growth_factor_binding",
  "cell_adhesion", "integrin_binding", "hydrolase_activity", "oxidoreductase_activity",
  "apoptotic_process", "ITxM", "CPI_candidate", "antibodies_bought",
  "meta_surfaceome",
  "APMS_combined_NS", "APMS_combined_NS_CD3" ,"APMS_combined_meta", "APMS_CD3z_FC2", "APMS_CD3e_FC2", "APMS_ZAP70_FC2", "APMS_Lck", "APMS_FC2_overlap")

cytoscape_column_yellow <- "APMS_combined_NS"  # alternative
table(cytoscape_info_annotation_full[[cytoscape_column_yellow]] )
cytoscape_info_annotation_full_qid <- cytoscape_info_annotation_full %>%
  mutate(cytoscape_fill = case_when(cytoscape_info_annotation_full[[cytoscape_column_yellow]] == "Yes" ~ "#ffff00", TRUE ~ "#bfbfbf"),
         cytoscape_label = "#000000")
# append info annotation table to cytoscape node table (existing colums are updated)
loadTableData(
  data = cytoscape_info_annotation_full_qid,
  data.key.column = "entry",       # Key column in your data frame (node IDs)
  table = "node",                  # Load into node table
  table.key.column = "entry"       # Corresponding key column in Cytoscape node table (usually "name")
)


## complex_portal ppi mapping ####################################
sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/thesis_figures_functions", pattern = "\\.R$", full.names = TRUE), source)
# complex_portal_ppi <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_human_ppi_MG20250306.csv")
# comp_cp_overall_LUX <- complex_compp(query_list = sort(unique(data_LUX_prot_diff_v31_meta_full                                                      %>% pull(entry)) ), set = "overall")  
comp_cp_nCD4  <- complex_compp(query_list = sort(unique(v31_LUX_data_prot_diff_abundance_sigup %>% dplyr::filter(plot_heading == "Naive CD4+" ) %>% pull(entry)) ), set = "Naive CD4+")  
comp_cp_nnCD4 <- complex_compp(query_list = sort(unique(v31_LUX_data_prot_diff_abundance_sigup %>% dplyr::filter(plot_heading == "Memory CD4+") %>% pull(entry)) ), set = "Memory CD4+")  
comp_cp_nCD8  <- complex_compp(query_list = sort(unique(v31_LUX_data_prot_diff_abundance_sigup %>% dplyr::filter(plot_heading == "Naive CD8+" ) %>% pull(entry)) ), set = "Naive CD8+")  
comp_cp_nnCD8 <- complex_compp(query_list = sort(unique(v31_LUX_data_prot_diff_abundance_sigup %>% dplyr::filter(plot_heading == "Memory CD8+") %>% pull(entry)) ), set = "Memory CD8+")  

## biogrid ppi mapping - could be used to add additional links not present in STRING

complex_portal_ppi <- rbind(comp_cp_nCD4, comp_cp_nnCD4, comp_cp_nCD8, comp_cp_nnCD8) %>% 
  dplyr::filter(complex_recall_cp >= (2/3) ) %>%
  separate(matches, into = c("A_entry", "B_entry"), sep = ", ") %>%
  dplyr::select(A_entry, B_entry) %>%
  left_join(proteome_20250804 %>% dplyr::select(entry, string), by = c("A_entry" = "entry")) %>%
  dplyr::rename("string_A" = "string") %>%
  left_join(proteome_20250804 %>% dplyr::select(entry, string), by = c("B_entry" = "entry")) %>%
  dplyr::rename("string_B" = "string") %>%
  mutate(name = paste(string_A, " (ppp) ", string_B),
         name = gsub(";", "", name),  # remove ";" which is contained in upsp file
         edge_color           = "blue",
         edge_weight          = 3,
         shared_name          = name,
         interaction          = "ppp",
         `shared interaction` = "ppp") %>%
  dplyr::select(-c(A_entry, B_entry, string_A, string_B))


loadTableData(
  data = complex_portal_ppi,
  data.key.column  = "name",      # Key column in your data frame (edge IDs)
  table            = "edge",      # Load into edge table
  table.key.column = "name"       # Corresponding key column in Cytoscape edge table (usually "name")
)


setEdgeColorMapping(
  table.column        = 'relationship',                     # Edge data column
  table.column.values = c('activation', 'inhibition'),
  colors              = c('#A0A0A0', '#FF7755')
)
setEdgeLineWidthMapping(
  table.column        = 'weight',                     # Numeric edge attribute
  table.column.values = c(1, 5),               # Range of edge weights
  widths              = c(1, 10)                            # Corresponding widths
)

# cytoscape_edge_table <- getTableColumns(table = "edge") %>%  # gets cytoscape edge table
#   mutate(edge_width = 1,
#          edge_color = )
# 
# 
# setdiff(names(complex_portal_ppi), names(cytoscape_edge_table))
# 


# # enrichment analysis background - not working currently... - load into cytoscape. manually "stringyfy" network tools>STRING>...
# protein_bg_v31 <- data.frame(id = c(read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v31_LUX_FP20_HoxHoxox_semi_6aa__4ss/_output_1-2_ludo_adjp_0.7string_shs2.27_ss_meta/_data_prot_level.csv") %>% pull(entry), 
#   v31_LUX_data_prot_diff_abundance_sigup_ssEXTENDED %>% pull(entry)) %>% unique(),   stringsAsFactors = FALSE  )
# createNetworkFromDataFrames(nodes = protein_bg_v31, title = "LUX protein ID background", collection = "LUX protein ID background")







# ##  OmicsVisualizer ------------
# # annotation_file <- normalizePath("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/cytoscape_info_annotation_TCRsignaling.csv")
# annotation_file <- normalizePath("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Chapter_3_R/cytoscape_info_annotation.csv")
# # OmicsVisualizer: load table  
# commandsPOST(paste0('ov load file="', annotation_file, '" newTableName="MyOmicsData"'))
# # OmicsVisualizer: connect table to nodes
# commandsPOST('ov connect mappingColNet="entry" mappingColTable="entry"')
# # Define your annotation column and desired colors
# column = "tetraspaning_enr_microdo"
# unique(cytoscape_info_annotation_full[[column]])
# setNodeColorMapping(
#   table.column        = column,
#   table.column.values = c("Yes", "No"),             #   c("Minus" , "Other", "Plus" ),
#   colors              = c("orange"    , "grey" ),   # c("red"   , "grey" , "green"),   # Choose distinct colors (hex or named),
#   mapping.type        = "d",  # 'd' for discrete
#   style.name          = "AP-MS"
# )
# cytoscape_shape_column = "cytoscape_shape"
# unique(cytoscape_info_annotation_full[[cytoscape_shape_column]])
# setNodeShapeMapping(table.column        = cytoscape_shape_column, 
#                     table.column.values = unique(cytoscape_info_annotation_full[[cytoscape_shape_column]]),
#                     # mapping.type        = "p",  # 'd' for discrete
#                     style.name          = "potato")#, identifier_values, shape_values)


 


# cytoscape_bait_network annotation of supporting data

cytoscape_bait_network <- cytoscape_bait_network %>%
  mutate(string_TCR = case_when(entry_name %in% string_targets$TCR ~ "x", TRUE ~ ""),
         string_CD4 = case_when(entry_name %in% string_targets$CD4 ~ "x", TRUE ~ ""),
         string_CD8 = case_when(entry_name %in% string_targets$CD8 ~ "x", TRUE ~ ""))













#





################################################################################################################################################
# __________________________________________________________________________________________________________________________________________________________________________________________________________________
# Figure Fig3.3.1_sup1: Surfaceome contriburtion plots (Volcano sig-up) ..........................................................................................................................................................................................................................................
surfaceome_sigup_table <- v31_LUX_data_prot_diff_abundance_sigup %>%
  mutate(meta_surfaceome = case_when(entry_name %in% abTCR_chains[!abTCR_chains %in% c("CD3D_HUMAN", "CD3G_HUMAN", "CD3Z_HUMAN")] ~ "abTCR", TRUE ~ meta_surfaceome)) %>%
  # mutate(comparison = gsub("0.5", "0,5", comparison)) %>%
  group_by(plot_heading, meta_surfaceome) %>%
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
Fig3.3.1_sup1 <- surfaceome_sigup_table %>%
  pivot_longer(cols = c(Other, Surfaceome, abTCR), names_to = "Annotation", values_to = "count") %>%
  mutate(Annotation   = factor(Annotation, levels = c("Other", "abTCR", "Surfaceome")),
         plot_heading = factor(plot_heading, levels = c("Meta panT", 
                                                        "Meta CD4+", "Meta CD8+", "Meta Naive", "Meta Memory", 
                                                        "Naive CD4+", "Memory CD4+", "Naive CD8+", "Memory CD8+"))
  ) %>%
  ggplot(aes(x = plot_heading, y = count, fill = Annotation)) +
  geom_bar(stat = "identity") +
  labs(title = "LUX Enrichment Annotation", 
       x     = "TCR LUX Condition", 
       y     = "Significantly\nEnriched Proteins [%]",
       fill  = "Category") +
  scale_fill_manual(values = c("Other" = "darkgrey", "Surfaceome" = "#cc0000", "abTCR" = "#0000EE")) +
  plot_theme()

Fig3.3.1_sup1
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.1_sup1.png",
  plot = Fig3.3.1_sup1,
  width  = 44,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

# __________________________________________________________________________________________________________________________________________________________________________________________________________________
# Fig3.3.1_sup2: abTCR-LUX  ..........................................................................................................................................................................................................................................

upset_data <- bind_rows( # ""  ""  "" "" ""   ""   ""   ""  ""
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Naive CD4+")  %>% pull(entry) %>% unique(),     set = "Naive CD4+"),
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Naive CD8+")  %>% pull(entry) %>% unique(),     set = "Naive CD8+"),
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Memory CD4+")  %>% pull(entry) %>% unique(),    set = "Memory CD4+"),
  tibble(entry = v31_LUX_data_prot_diff_abundance_sigup  %>% filter(plot_heading == "Memory CD8+")  %>% pull(entry) %>% unique(),    set = "Memory CD8+"),
) %>% 
  group_by(entry) %>% 
  summarise(sets = list(set)) %>%  # Critical: list column of set memberships
  ungroup() 

Fig3.3.1_sup2 <- upset_data %>% 
  ggplot(aes(x = sets)) +
  geom_bar(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_upset(
    sets = unique(v31_LUX_data_prot_diff_abundance_sigup$plot_heading),
    name = "",
    n_intersections = 30  ) +
  labs(y     = "Intersection",
       title = "Protein Community Intersection") +
  plot_theme() # +
# theme(axis.text.y = element_text(size = 14))

Fig3.3.1_sup2
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig3.3.1_sup2.png",
  plot = Fig3.3.1_sup2,
  width  = 15,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  

sapply(list.files(path = "/Users/mgesell/Desktop/currentR/git/surfaceome_hybrid_script/thesis_figures_functions", pattern = "\\.R$", full.names = TRUE), source)



##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
#     ___ _                 _              _  _   
#   / ___| |__   __ _ _ __ | |_ ___ _ __  | || |  
#  | |   | '_ \ / _` | '_ \| __/ _ \ '__| | || |_ 
#  | |___| | | | (_| | |_) | ||  __/ |    |__   _|
#  \____|_| |_|\__,_| .__/ \__\___|_|       |_|  
#                   |_|                          
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
# AAV poi
# core known TCR community
AAV_poi <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv") %>% 
  filter(AAV2_poi != "") %>% 
  pull(AAV2_poi)
AAV_Pillay2016 <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/POI_lists/POI_lists.csv") %>% 
  filter(AAV2_Pillay2016 != "") %>% 
  pull(AAV2_Pillay2016)
# CSC data HeLa, HepG2, HEK, K562
AAV_CSC_cell_line <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/AAV_CSC_v2_cell_lines/_output_1-2_ludo_adjp_0.7string_shs2.27/_data_prot_level.csv")

# Hela AAV-SOG timecourse
AAV_LUX_Hela        <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/AAV_1Y_AAV2-Hela_timecourse/_output_1-2_no-impute_shs2.27/_data_prot_level.csv")
AAV_LUX_Hela_diff   <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/AAV_1Y_AAV2-Hela_timecourse/_output_1-2_no-impute_shs2.27/_data_prot_diff_abundance.csv")
#________________________________________________________________________________________________________________________________________________________________________________________________________

### Fig4.3.1_a Cell line expression of AAV pois ........................................................................................................
AAV_CSC_cell_line_poi <- AAV_CSC_cell_line %>%
  filter(entry_name %in% AAV_poi) %>%
  dplyr::select(condition, entry_name, condition_median_imp_log2_intensity) %>%
  mutate(entry_name = gsub("_HUMAN", "", entry_name))

Fig4.3.1_a <-   ggplot(AAV_CSC_cell_line_poi, aes(x = entry_name, y = condition, fill = condition_median_imp_log2_intensity)) +
  geom_tile(color = "white", linewidth = 0.3) +  # Add white borders between tiles
  scale_fill_gradient(
    name = "log2(median\nprotein intensity)", 
    low  = "#e6e6ff", 
    high = "#001a4d"   ) + # Dark blue gradient (adjust as needed)
  labs(
    x = "Protein", 
    y = "Cell line", 
    title = "CSC protein signal\nfor known AAV interactors"  ) +
  plot_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.0),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),  # Center title
    panel.grid = element_blank()  # Remove grid lines
  )
Fig4.3.1_a
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig4.3.1_a.png",
  plot = Fig4.3.1_a,
  width  = 17,  # document is 16 cm wide             (before 12cm used)
  height = 8,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
) 

### Fig4.3.1_b Cell line expression of AAV pois ........................................................................................................
# overlap of glyco IDs for the 4 CSC profiled cell lines
upset_data <- bind_rows(
  tibble(entry = AAV_CSC_cell_line %>% filter(condition == "HeLa", imputed_prot_intensity_log2 > 0) %>% pull(entry) %>% unique(), 
         set = "HeLa"),
  tibble(entry = AAV_CSC_cell_line %>% filter(condition == "HepG2", imputed_prot_intensity_log2 > 0) %>% pull(entry) %>% unique(),         
         set = "HepG2"),
  tibble(entry = AAV_CSC_cell_line %>% filter(condition == "K562", imputed_prot_intensity_log2 > 0) %>% pull(entry) %>% unique(),        
         set = "K562"),
  tibble(entry = AAV_CSC_cell_line %>% filter(condition == "HEK", imputed_prot_intensity_log2 > 0) %>% pull(entry) %>% unique(), 
         set = "HEK") 
  # tibble(entry = intersect(intersect(CSPA, CSPA_Jurkat), panT_csc_ids), # to not display full CSPA but show which proteins already annotate in CSPA  
  #                               set = "CSPA subset"),
) %>% 
  group_by(entry) %>% 
  summarise(sets = list(set)) %>%  # Critical: list column of set memberships
  ungroup() 

Fig4.3.1_b <- upset_data %>% 
  ggplot(aes(x = sets)) +
  geom_bar(fill = "black", color = "white", linewidth = 0.3) +
  scale_x_upset(
    sets = c("HeLa", "HepG2", "K562", "HEK"),
    name = "",
    n_intersections = 20
  ) +
  labs(
    y     = "Intersection",
    title = "Protein Overlap"
  ) +
  plot_theme() +
  theme(axis.text.y = element_text(size = 14))

Fig4.3.1_b
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig4.3.1_b.png",
  plot = Fig4.3.1_b,
  width  = 9,  # document is 16 cm wide             (before 12cm used)
  height = 8.00,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
)  


### Fig4.3.1_h LUX AAV poi protein signal  ........................................................................................................
# # median signal of meta_surfaceome proteins (surfaceome abundance control for heatmaps)
# AAV_LUX_Hela_median_surfaceome_intensity_per_condition <- AAV_LUX_Hela %>%
#   filter(meta_surfaceome == "yes") %>%
#   dplyr::select(entry_name, condition, condition_median_imp_log2_intensity) %>%
#   group_by(condition) %>%
#   mutate(condition_median_imp_log2_intensity = median(condition_median_imp_log2_intensity,  na.rm = TRUE),
#          entry_name = "median_surface_ctrl") %>%
#   ungroup() %>%
#   distinct()
# median signal of CSC identified HeLA proteins (surfaceome abundance control for heatmaps)
AAV_LUX_Hela_median_CSCidentified_intensity_per_condition <- AAV_LUX_Hela %>%
  filter(entry_name %in% (AAV_CSC_cell_line %>% filter(condition == "HeLa") %>% pull(entry_name) %>% unique())) %>%
  dplyr::select(entry_name, condition, condition_median_imp_log2_intensity) %>%
  group_by(condition) %>%
  mutate(condition_median_imp_log2_intensity = median(condition_median_imp_log2_intensity,  na.rm = TRUE),
         entry_name = "median_CSC") %>%
  ungroup() %>%
  distinct()

AAV_LUX_Hela <- AAV_LUX_Hela %>%
  dplyr::select(condition, entry_name, condition_median_imp_log2_intensity) %>%
  rbind(AAV_LUX_Hela_median_CSCidentified_intensity_per_condition) %>%
  # rbind(AAV_LUX_Hela_median_surfaceome_intensity_per_condition) %>%
  mutate(condition  = factor(condition , levels = rev(c("10_NL", "10", "20", "30", "40", "60")))) %>%
  # for 0-1 scaling:
  group_by(entry_name) %>%
  mutate(normalized_expression = percent_rank(condition_median_imp_log2_intensity)) %>%
  ungroup()

AAV_LUX_Hela_poi <- AAV_LUX_Hela %>%
  filter(entry_name %in% c(AAV_poi, "CAPSD_AAV2S", "median_surface_ctrl", "median_CSC"), condition != "") %>%
  mutate(entry_name = gsub("_HUMAN"     , ""        , entry_name),
         entry_name = gsub("CAPSD_AAV2S", "AAV2-VP1", entry_name),
         entry_name = factor(entry_name,  levels = c("median_CSC", sort(setdiff(unique(entry_name), "median_CSC")))))

Fig4.3.1_h <-   ggplot(AAV_LUX_Hela_poi, aes(x = entry_name, y = condition, fill = normalized_expression)) +
  geom_tile(color = "white", linewidth = 0.3) +  # Add white borders between tiles
  scale_fill_gradient(
    name = "Normalized log2\n(median protein\nintensity)", 
    low  = "#e6e6ff", 
    high = "#001a4d"   ) + # Dark blue gradient (adjust as needed)
  labs(
    x = "Protein", 
    y = "Cell line", 
    title = "LUX protein signal for known AAV2 interactors"  ) +
  plot_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.0),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),  # Center title
    panel.grid = element_blank()  # Remove grid lines
  )
Fig4.3.1_h
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig4.3.1_h.png",
  plot = Fig4.3.1_h,
  width  = 30,  # document is 16 cm wide             (before 12cm used)
  height = 9,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
) 
# ....
AAV_LUX_Hela_poi_pillay <- AAV_LUX_Hela %>%
  filter(entry_name %in% AAV_Pillay2016, condition != "") %>%
  mutate(entry_name = gsub("_HUMAN", "", entry_name))

Fig4.3.1_i <-   ggplot(AAV_LUX_Hela_poi_pillay, aes(x = entry_name, y = condition, fill = normalized_expression)) +
  geom_tile(color = "white", linewidth = 0.3) +  # Add white borders between tiles
  scale_fill_gradient(
    name = "Normalized log2\n(median protein\nintensity)", 
    low  = "#e6e6ff", 
    high = "#001a4d"   ) + # Dark blue gradient (adjust as needed)
  labs(
    x = "Protein", 
    y = "Cell line", 
    title = "LUX protein signal for top hits in Pillay 2016 KO screen"  ) +
  plot_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.0),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),  # Center title
    panel.grid = element_blank()  # Remove grid lines
  )
Fig4.3.1_i
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig4.3.1_i.png",
  plot = Fig4.3.1_i,
  width  = 30,  # document is 16 cm wide             (before 12cm used)
  height = 8,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
) 
#.....
# GO Term	Biological Role	AAV Relevance
# GO:0019062 virion attachment to host cell
# GO:0019065 receptor-mediated endocytosis of virus by host cell
# GO:0001618 virus receptor activity
# GO:0016032 viral process
AAV_go_poi <- proteome_upsp_202501 %>% 
  filter(grepl("GO:0019062|GO:0019065|GO:0001618", gene_ontology_go)) %>%
  dplyr::select(entry_name) %>% 
  pull(entry_name)

AAV_LUX_Hela_poi_go <- AAV_LUX_Hela %>%
  filter(entry_name %in% AAV_go_poi) %>%
  mutate(entry_name = gsub("_HUMAN", "", entry_name))

Fig4.3.1_j <-   ggplot(AAV_LUX_Hela_poi_go, aes(x = entry_name, y = condition, fill = normalized_expression)) +
  geom_tile(color = "white", linewidth = 0.3) +  # Add white borders between tiles
  scale_fill_gradient(
    name = "Normalized log2\n(median protein\nintensity)", 
    low  = "#e6e6ff", 
    high = "#001a4d"   ) + # Dark blue gradient (adjust as needed)
  labs(
    x = "Protein", 
    y = "Cell line", 
    title = "LUX protein signal for viral GO interactors"  ) +
  plot_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.0),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),  # Center title
    panel.grid = element_blank()  # Remove grid lines
  )
Fig4.3.1_j
ggsave(
  filename = "/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/thesis_figures/Fig4.3.1_j.png",
  plot = Fig4.3.1_j,
  width  = 30,  # document is 16 cm wide             (before 12cm used)
  height = 8,  # 4/3 width/high ratio is common      (before 8 cm used)
  units  = "cm",
  dpi    = 300    # default for good quality
) 


library(gitcreds)
gitcreds_set()

