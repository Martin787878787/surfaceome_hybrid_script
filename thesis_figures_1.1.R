#Load libraries required ----------------------------------------------------------------------------------------------------------
rm(list = ls())           # Purge workspace
set.seed(123)
script_version = "_1.1" # version stamp on output directory
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
v25_LLOQ_CSC <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_CSC_panT_lowinput/_output_1-2_ludo_adjp_0.7string_shs2.26/_data_prot_level.csv") %>%
  mutate(condition = gsub("0_5", "0.5", condition)) 

v25_LLOQ_CSC_pep <- read_protti("/Users/mgesell/Desktop/currentR/2025-01__local_reanalysis_paper_candi_experiements/v25_CSC_panT_lowinput/_output_1-2_ludo_adjp_0.7string_shs2.26/_data_feature_level.csv") %>%
  mutate(condition = gsub("0_5", "0.5", condition)) %>%
  mutate(N_deam_psm_n = str_count(peptide_sequence_mod,   "N\\[0.9840\\]"))

table(v25_LLOQ_CSC_pep %>% dplyr::filter(csc_signature_psm == "yes", !is.na(ng_sites)) %>%  select(entry, csc_signature_psm_n, peptide_sequence_mod, ng_sites) %>% distinct() %>% pull(csc_signature_psm_n))

# QC of deamidation signatures - how many in NxSTC signature - should be majority otherwise likely chem deam possible too
NsiteQC <- v25_LLOQ_CSC_pep %>% select(peptide_sequence_mod, N_deam_psm_n, csc_signature_psm_n) %>% 
  distinct() %>% # keep only unique modified peptides (no charge states)
  select(-peptide_sequence_mod) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Count")
  
ggplot(NsiteQC %>%
         group_by(Type) %>%
         summarise(Sum = sum(Count)), aes(x = Type, y = Sum, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("N_deam_psm_n" = "#ff7f0e", "csc_signature_psm_n" = "black")) +
  plot_theme()

bad_csc_sites   <- v25_LLOQ_CSC_pep %>% filter(csc_signature_psm == "no" , N_deam_psm_n >0) %>% select(entry_name, semi_stripped_peptide_sequence, meta_surfaceome, cspa_2015) %>% distinct()
weird_csc_sites <- v25_LLOQ_CSC_pep %>% filter(csc_signature_psm == "yes", N_deam_psm_n >1) %>% # require at least 2 sites (CSC and non-csc site)
  select(entry_name, semi_stripped_peptide_sequence, meta_surfaceome, cspa_2015) %>% distinct()
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
  select(entry, entry_name, condition, csc_signature_psm, peptide_sequence_mod) %>%
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

paste("median CV 5e6 cells CSC:  ", round(v25_LLOQ_CSC %>% filter(condition == "5") %>% select(entry, cv_protein) %>% distinct() %>% pull(cv_protein) %>% median(na.rm = TRUE),1))
paste("median CV 1e6 cells CSC:  ", round(v25_LLOQ_CSC %>% filter(condition == "1") %>% select(entry, cv_protein) %>% distinct() %>% pull(cv_protein) %>% median(na.rm = TRUE),1))
paste("median CV 5e5 cells CSC:  ", round(v25_LLOQ_CSC %>% filter(condition == "0.5") %>% select(entry, cv_protein) %>% distinct() %>% pull(cv_protein) %>% median(na.rm = TRUE),1))

# Fig2.3.1_f: CSC Precursor Signal ......................................................................................................................
Fig2.3.1_f <- v25_LLOQ_CSC_pep %>%
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
novel_Ng_proteins_surface_annotated  <- evaluate_csc_sites_result[[2]]
cat(novel_Ng_proteins_surface_annotated$entry %>% unique(), sep = "\n")
#
query_protein = "P05556"
query_protein_peptides <- extract_csc_peptides_for_protti(csc_output__data_raw_up_surf_gly_nZ = v25_LLOQ_CSC_pep,
                                                          csc_sites_plus_upsp_annotation      = csc_sites_plus_upsp_annotation, 
                                                          query_protein                       = query_protein)
#
Ng_summary <- csc_sites_plus_upsp_annotation %>%
  select(entry, site, summary) %>%
  distinct() %>%
  mutate(summary = gsub("NOVEL site", "novel", summary)) %>%
  rename(CSC_site = site) %>%
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



# =============================================================================================================================================================
# Figure 2.3.1_k: TCR-LUX LLOQ panT  ..........................................................................................................................................................................................................................................
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
#     ____ _                 _              _____ 
#   / ___| |__   __ _ _ __ | |_ ___ _ __  |___ / 
#  | |   | '_ \ / _` | '_ \| __/ _ \ '__|   |_ \ 
#  | |___| | | | (_| | |_) | ||  __/ |     ___) |
#  \____|_| |_|\__,_| .__/ \__\___|_|    |____/ 
#                   |_|                         
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################






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
  select(condition, entry_name, condition_median_imp_log2_intensity) %>%
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
#   select(entry_name, condition, condition_median_imp_log2_intensity) %>%
#   group_by(condition) %>%
#   mutate(condition_median_imp_log2_intensity = median(condition_median_imp_log2_intensity,  na.rm = TRUE),
#          entry_name = "median_surface_ctrl") %>%
#   ungroup() %>%
#   distinct()
# median signal of CSC identified HeLA proteins (surfaceome abundance control for heatmaps)
AAV_LUX_Hela_median_CSCidentified_intensity_per_condition <- AAV_LUX_Hela %>%
  filter(entry_name %in% (AAV_CSC_cell_line %>% filter(condition == "HeLa") %>% pull(entry_name) %>% unique())) %>%
  select(entry_name, condition, condition_median_imp_log2_intensity) %>%
  group_by(condition) %>%
  mutate(condition_median_imp_log2_intensity = median(condition_median_imp_log2_intensity,  na.rm = TRUE),
         entry_name = "median_CSC") %>%
  ungroup() %>%
  distinct()

AAV_LUX_Hela <- AAV_LUX_Hela %>%
  select(condition, entry_name, condition_median_imp_log2_intensity) %>%
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
  select(entry_name) %>% 
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



