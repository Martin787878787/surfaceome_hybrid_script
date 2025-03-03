# belongs to shs_downstream
load_uniprot_and_annotate <- function(proteome) {
  
  # human proteome
  proteome  <- proteome %>%
    # annotate PCOIs (protein categories of interest)
    mutate(   # annotate surfaceome categories of interest (extend list if needed and communicate suggestion to Martin)
      t_tcr_signaling_plus       = as.integer(str_detect(gene_ontology_i_ds, "GO:0050852")),    # 114  reviewed human proteins in UP 2025-01-17
      t_tcr_signaling_minus      = as.integer(str_detect(gene_ontology_i_ds, "GO:0050860")),    # 26   reviewed human proteins in UP 2025-01-17
      t_act_plus                 = as.integer(str_detect(gene_ontology_i_ds, "GO:0042110")),    # 337  reviewed human proteins in UP 2025-01-17
      t_act_minus                = as.integer(str_detect(gene_ontology_i_ds, "GO:0050868")),    # 130  reviewed human proteins in UP 2025-01-17
      
      r_complex         = as.integer(str_detect(gene_ontology_i_ds, "GO:0043235")),             # 545  reviewed human proteins in UP 2025-01-17
      r_clustering      = as.integer(str_detect(gene_ontology_i_ds, "GO:0043113")),             # 48   reviewed human proteins in UP 2025-01-17
      r_transactivation = as.integer(str_detect(gene_ontology_i_ds, "GO:0035624")),             # 4    reviewed human proteins in UP 2025-01-17
      r_internalization = as.integer(str_detect(gene_ontology_i_ds, "GO:0031623")),             # 75   reviewed human proteins in UP 2025-01-17
      r_recycling       = as.integer(str_detect(gene_ontology_i_ds, "GO:0001881")),             # 23   reviewed human proteins in UP 2025-01-17
      r_metabolism      = as.integer(str_detect(gene_ontology_i_ds, "GO:0043112")),             # 42   reviewed human proteins in UP 2025-01-17
      r_growth_factor_binding    = as.integer(str_detect(gene_ontology_i_ds, "GO:0019838")),    # 129  reviewed human proteins in UP 2025-01-17 
      r_antigen_binding          = as.integer(str_detect(gene_ontology_i_ds, "GO:0003823")),    # 184  reviewed human proteins in UP 2025-01-17
      r_mhc_complex              = as.integer(str_detect(gene_ontology_i_ds, "GO:0042611")),    # 25   reviewed human proteins in UP 2025-01-17
      r_growth_factor_binding    = as.integer(str_detect(gene_ontology_i_ds, "GO:0019838")),    # 129  reviewed human proteins in UP 2025-01-17
      r_cytokine_receptor        = as.integer(str_detect(gene_ontology_i_ds, "GO:0004896")),    # 96   reviewed human proteins in UP 2025-01-17
      r_cell_adhesion            = as.integer(str_detect(gene_ontology_i_ds, "GO:0007155")),    # 968  reviewed human proteins in UP 2025-01-17
      r_integrin_binding         = as.integer(str_detect(gene_ontology_i_ds, "GO:0005178")),    # 159  reviewed human proteins in UP 2025-01-17
      
      s_signal_transduction      = as.integer(str_detect(gene_ontology_i_ds, "GO:0007165")),    # 4821 reviewed human proteins in UP 2025-01-17
      s_kinase                   = as.integer(str_detect(gene_ontology_i_ds, "GO:0004672")),    # 563  reviewed human proteins in UP 2025-01-17
      s_phosphatase              = as.integer(str_detect(gene_ontology_i_ds, "GO:0016791|GO:0004725|GO:0016311")), # 272+98+189 reviewed human proteins in UP 2025-01-17
      
      f_immune_response          = as.integer(str_detect(gene_ontology_i_ds, "GO:0006955")),    # 1667 reviewed human proteins in UP 2025-01-17
      f_transporter_activity     = as.integer(str_detect(gene_ontology_i_ds, "GO:0005215")),    # 1229 reviewed human proteins in UP 2025-01-17
      f_ion_channel              = as.integer(str_detect(gene_ontology_i_ds, "GO:0005216")),    # 427  reviewed human proteins in UP 2025-01-17
      f_hydrolase_activity       = as.integer(str_detect(gene_ontology_i_ds, "GO:0016787")),    # 2440 reviewed human proteins in UP 2025-01-17
      f_oxidoreductase_activity  = as.integer(str_detect(gene_ontology_i_ds, "GO:0016491")),    # 714  reviewed human proteins in UP 2025-01-17
      f_apoptotic_process        = as.integer(str_detect(gene_ontology_i_ds, "GO:0006915")),    # 1043 reviewed human proteins in UP 2025-01-17
      f_redox_homeostasis        = as.integer(str_detect(gene_ontology_i_ds, "GO:0045454")),    # 40   reviewed human proteins in UP 2025-01-17
      f_protein_targeting        = as.integer(str_detect(gene_ontology_i_ds, "GO:0006605")),    # 250  reviewed human proteins in UP 2025-01-17
      f_complex_organization     = as.integer(str_detect(gene_ontology_i_ds, "GO:0043933")),    # 1443 reviewed human proteins in UP 2025-01-17 
      
      pf_tetraspanin_galectin    = as.integer(entry_name %in% protein_families$tetraspanin & entry_name != "" | entry_name %in% protein_families$galectin    & entry_name != ""), #   # ifelse(..., "tetraspanin ", ""),    # 33 proteins 2025-01-17
      pf_integrin         = as.integer(entry_name %in% protein_families$integrin    & entry_name != ""),#            ..., "integrin "   , ""),    # 25 proteins 2025-01-17
      pf_disintegrin      = as.integer(entry_name %in% protein_families$disintegrin & entry_name != ""),#            ..., "disintegrin ", ""),    # 44 proteins 2025-01-17
      pf_selectin         = as.integer(entry_name %in% protein_families$selectin    & entry_name != ""),#            ..., "selectin "   , ""),    # 3  proteins 2025-01-17
      pf_icam             = as.integer(entry_name %in% protein_families$icam        & entry_name != ""),#            ..., "icam "       , ""),    # 6  proteins 2025-01-17
      pf_gpcr             = as.integer(str_detect(gene_ontology_i_ds, "GO:0004930")),        # 876  reviewed human proteins in UP 2025-01-17
      pf_tnfr             = as.integer(str_detect(gene_ontology_i_ds, "GO:0005031")),        # 10   reviewed human proteins in UP 2025-01-17
      
      poi_LUX_targets.CPIs         = ifelse(entry_name %in% poi_reference$LUX_targets.CPIs[poi_reference$LUX_targets.CPIs                    != ""], 1, 0),
      poi_tcr_signaling_GO.0050852 = ifelse(entry_name %in% poi_reference$tcr_signaling_GO.0050852[poi_reference$tcr_signaling_GO.0050852    != ""], 1, 0),
      poi_abTCR_chain              = ifelse(entry_name %in% poi_reference$abTCR_chains_cd3_mgmanual[poi_reference$abTCR_chains_cd3_mgmanual  != ""], 1, 0),
      poi_ydTCR_chain              = ifelse(entry_name %in% poi_reference$ydTCR_chains_cd3_mgmanual[poi_reference$ydTCR_chains_cd3_mgmanual  != ""], 1, 0),
      poi_CPIs                     = ifelse(entry_name %in% poi_reference$CPIs[poi_reference$CPIs                                            != ""], 1, 0),
      poi_POI_Ben                  = ifelse(entry_name %in% poi_reference$POI_Ben[poi_reference$POI_Ben                                      != ""], 1, 0),
      poi_FACS_candidates          = ifelse(entry_name %in% poi_reference$FACS_candidates[poi_reference$FACS_candidates                      != ""], 1, 0),
      poi_Tcell_Act_GO.0042110     = ifelse(entry_name %in% poi_reference$Tcell_Act_GO.0042110[poi_reference$Tcell_Act_GO.0042110            != ""], 1, 0),
      poi_Thermo_Marker_Tact       = ifelse(entry_name %in% poi_reference$Thermo_Marker_Tact[poi_reference$Thermo_Marker_Tact                != ""], 1, 0),
      poi_Thermo_Marker_Tsubset    = ifelse(entry_name %in% poi_reference$Thermo_Marker_Tsubset[poi_reference$Thermo_Marker_Tsubset          != ""], 1, 0),
      poi_Tact_markers             = ifelse(entry_name %in% poi_reference$Tact_markers[poi_reference$Tact_markers                            != ""], 1, 0),
      poi_BindingBlockListPharmacoscopy_shilts2022_fig4d = ifelse(entry_name %in% poi_reference$BindingBlockListPharmacoscopy_shilts2022_fig4d[poi_reference$BindingBlockListPharmacoscopy_shilts2022_fig4d != ""], 1, 0),
      poi_CD_antigen               = ifelse(entry_name %in% poi_reference$CD_antigen[poi_reference$CD_antigen                                != ""], 1, 0),
      
      string_CD4 = ifelse(entry_name %in% string[string$CD4_LUX == 1,]$protein2_entry_name , 1, 0),  # 551 of these 331 in proteome 2025-01-17  
      string_CD8 = ifelse(entry_name %in% string[string$CD8_LUX == 1,]$protein2_entry_name , 1, 0),  # 507 of these 287 in proteome 2025-01-17
      string_TCR = ifelse(entry_name %in% string[string$TCR_LUX == 1,]$protein2_entry_name , 1, 0)   # 821 of these 372 in proteome 2025-01-17
    )
  # define pcoi categories (negative so code still works when further pcoi categroies added) for loops
  pcoi_list <- colnames( proteome %>% select(-c("entry", "reviewed", "entry_name", "protein_names", "gene_names", "organism", "length", 
                                                "gene_ontology_go", "gene_ontology_biological_process", "gene_ontology_cellular_component", 
                                                "gene_ontology_molecular_function", "gene_ontology_i_ds", "subcellular_location_cc", "transmembrane", 
                                                "function_cc", "pathway", "gene_names_primary", "disulfide_bond", "glycosylation", "lipidation", "intramembrane", 
                                                "topological_domain", "ensembl", "sequence", "interacts_with")))
  return(proteome)
}