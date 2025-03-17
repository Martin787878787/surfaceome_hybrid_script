# belongs to shs_downstream
complex_compp <- function(query_list, set) {

  # #### subset and unboxing of complex portal data ##############################################################################################################################
  # complex_portal__ <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/ComplexPortal_human9606_20250218.tsv")   %>%
  #   filter(taxonomy_identifier == "9606") %>%
  #   dplyr::select(number_complex_ac, recommended_name, description, identifiers_and_stoichiometry_of_molecules_in_complex, experimental_evidence, complex_assembly)
  # 
  # ## load complex_portal complexes --------------------------------------------------------------------------------------------------------------------------------------------------------
  # # !!! CAVEAT i assume [] indicates "one or the other" in complex. so far did not account for this and just considered as ppis also withing bracket ... !!!
  # complex_portal_ <- complex_portal__   %>%
  #   mutate(subunits_uniprot_ids_cp = gsub("\\([0-9]+\\)", ""     , identifiers_and_stoichiometry_of_molecules_in_complex),
  #          # subunits_uniprot_ids_cp = gsub("\\[|\\]"         , ""     , subunits_uniprot_ids_cp),  # destroy sub-complex structure
  #          subunits_uniprot_ids_cp = gsub("\\|" , ","    , subunits_uniprot_ids_cp),  # destroy sub-complex structure
  #          subunits_uniprot_ids_cp = gsub("-PRO_"       , ",PRO_",  subunits_uniprot_ids_cp))  # replace - (which i think indicate protein RNA complex) with , to indicate complex members
  # 
  # # some complexes are described as complexes of complexes - in order to resolve this relationship to protein only annotation ...
  # # Function to expand complex identifiers
  # expand_complex <- function(complex_id, df) {
  #   # Locate the row where number_complex_ac matches the complex_id
  #   row <- df %>% filter(number_complex_ac == complex_id)
  #   # If the complex_id is not found, return the original complex_id
  #   if (nrow(row) == 0) {
  #     return(complex_id)
  #   }
  #   # Return the protein identifiers for that complex
  #   return(row$subunits_uniprot_ids_cp)
  # }
  # # Recursive replacement function
  # replace_complexes <- function(data) {
  #   new_data <- data %>%
  #     mutate(subunits_uniprot_ids_cp = sapply(subunits_uniprot_ids_cp, function(ids) {
  #       while (str_detect(ids, "CPX-")) {
  #         # Extract the first complex identifier found
  #         complex_id <- str_extract(ids, "CPX-\\d+")
  #         # Check if a complex_id was actually found
  #         if (is.na(complex_id)) {
  #           break # Exit the while loop if no complex ID is found
  #         }
  #         # Expand the complex identifier
  #         expanded_ids <- expand_complex(complex_id, data)
  #         # Replace the complex identifier with the expanded IDs
  #         ids <- str_replace(ids, complex_id, expanded_ids)
  #       }
  #       return(ids)
  #     }))
  #   return(new_data)
  # }
  # 
  # # Apply the recursive replacement
  # complex_portal_expanded <- replace_complexes(complex_portal_)
  # 
  # # ## remove CHEBI and RNAs from complexes
  # complex_portal_expanded$subunits_uniprot_ids_cp <- gsub("CHEBI:\\d+[,]?\\s*", "", complex_portal_expanded$subunits_uniprot_ids_cp)
  # #    # fyi these are unique terms according to perplexity
  # #    # CHEBI:15414: Acyclic dialdehyde that is glyoxal.
  # #    # CHEBI:15996: A N-acyl-L-amino acid that is N-acetylalanine.
  # #    # CHEBI:16238: A monocarboxylic acid amide resulting from the formal condensation of the carboxy group of nicotinic acid with ammonia.
  # #    # CHEBI:17489: A 2-hydroxy fatty acid that is 2-hydroxypalmitic acid.
  # #    # CHEBI:18420: A C18 steroid that is estratriol.
  # #    # CHEBI:21137: A purine nucleotide containing guanine as the nucleobase and ribose 5'-phosphate as the sugar phosphate moiety.
  # #    # CHEBI:26355: An aminopurine nucleotide that is deoxyadenosine monophosphate.
  # #    # CHEBI:27363: An oligosaccharide composed of two beta-D-galactopyranose units joined by a glycosidic linkage.
  # #    # CHEBI:29034: A monosaccharide derivative resulting from the formal condensation of D-galactose with one molecule of ammonia.
  # #    # CHEBI:29035: A hexosamine that is a component of hyaluronan and keratan sulfate.
  # #    # CHEBI:29105: A hexosamine that is D-glucosamine.
  # #    # CHEBI:29108: An N-acetylhexosamine that is N-acetyl-D-glucosamine.
  # #    # CHEBI:30408: A conjugate base of guanylic acid, arising from deprotonation of the phosphate groups.
  # #    # CHEBI:37565: A tertiary alpha-hydroxy ketone and tertiary alcohol. It has a role as a human metabolite.
  # #    # CHEBI:49601: A purine 5'-monophosphate having adenine as the nucleobase.
  # #    # CHEBI:49883: A ribonucleotide 5'-monophosphate having cytidine as the nucleobase.
  # #    # CHEBI:57692: A pyrimidine 5'-monophosphate having uracil as the nucleobase.
  # #    # CHEBI:58189: A compound having a structure based on that of cholesterol, but with the double bond between positions 5 and 6 saturated.
  # #    # CHEBI:58210: A purine 5'-diphosphate having guanine as the nucleobase.
  # #    # CHEBI:58937: A purine 5'-monophosphate having guanine as the nucleobase.
  # #    # CHEBI:59789: A 20-hydroxyecdysone.
  # #    # CHEBI:60240: A C27 steroid that is cholesta-5,7-dien-3beta-ol.
  # #    # CHEBI:64607: A purine 5'-triphosphate having guanine as the nucleobase.
  # #    # CHEBI:190135: An N-acyl-L-amino acid that is N-lauroylglycine.
  # complex_portal_expanded$subunits_uniprot_ids_cp <- gsub("URS\\d+[A-Z0-9]+_\\d+[,]?\\s*", "", complex_portal_expanded$subunits_uniprot_ids_cp)
  # # # URS000075C8FA_9606
  # # # URS000044DFF6_9606
  # # # URS000047A7F4_9606
  # # # URS000080E36F_9606
  # # # URS0000704D22_9606
  # # # URS00000F9D45_9606
  # # # URS000075D341_9606
  # # # URS0000ABD82A_9606
  # # # URS0000000A8C_9606
  # # # URS000058CD99_9606
  # # # URS00005CCF67_9606
  # # # URS0000A76E3A_9606
  # # # URS0000631BD4_9606
  # # # URS00006767A8_9606
  # # # URS0000715A86_9606
  # # # URS00006A1489_9606
  # # # URS00000478B7_9606
  # # # URS000013F331_9606
  # 
  # # write.csv(### do not overwrite !! ## complex_portal_expanded %>% filter(grepl("\\[", subunits_uniprot_ids_cp)),
  # #           "/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_anoying_combo_complexes.csv", row.names = FALSE)
  # manual_rescue <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_corum_anoying_combo_complexes.csv", head = TRUE) %>% dplyr::select(-integrated, -ignored)
  #       # remove this overly elaborate 100 combo copmlexes
  #       # #   CPX-5223	40S cytosolic small ribosomal subunit	Component of the ribosome, the site of protein biosynthesis resulting from translation of messenger RNA (mRNA). Acts as the decoding centre of the ribosome which brings mRNA and aminoacylated transfer (t)RNAs together, with the 16S ribosomal (r)RNA being required for the selection of the cognate tRNA.	P08708(0)|P08865(0)|P15880(0)|P23396(0)|P25398(0)|P39019(0)|P46781(0)|P46782(0)|P46783(0)|P60866(0)|P61247(0)|P62081(0)|P62241(0)|P62244(0)|P62249(0)|P62263(0)|P62266(0)|P62269(0)|P62273(0)|P62277(0)|P62280(0)|P62753(0)|P62841(0)|P62847(0)|P62851(0)|P62854(0)|P62857(0)|P62861(0)|P62979(0)|P63220(0)|P63244(0)|URS0000704D22_9606(0)|[P42677,Q71UM5](0)|[P62701](0)	-	-	P08708,P08865,P15880,P23396,P25398,P39019,P46781,P46782,P46783,P60866,P61247,P62081,P62241,P62244,P62249,P62263,P62266,P62269,P62273,P62277,P62280,P62753,P62841,P62847,P62851,P62854,P62857,P62861,P62979,P63220,P63244,	P42677,Q71UM5,	P62701
  #       # #   CPX-6904	Vacuolar proton translocating ATPase complex, ATP6V0A2 variant	Translocates protons across a lipid bilayer via an ATP-driven rotary mechanism, thus acidifing the lumen of its resident organelle. Membrane-bound ion transporters/proton exchangers use the pH gradient to sequester metal ions to the vacuole and other cellular organelles. The combined action of the V-ATPase and membrane transporters plays a key role in maintaining cellular homoeostasis. The ATP6V0A2 variant localizes to the Golgi apparatus.	O75787(1)|P27449(9)|P38606(3)|Q15904(1)|Q16864(1)|Q6P5S7(1)|Q99437(1)|Q9UI12(0)|Q9Y487(1)|Q9Y5K8(1)|[O15342,Q8NHE4](1)|[O95670,Q96LB4,O75348](3)|[P15313,P21281](3)|[P21283,Q8NEY4](1)|[Q8N8Y2,P61421](1)|[Q96A05,P36543](3)	-	Hetero 32-mer	O75787,P27449,P38606,Q15904,Q16864,Q6P5S7,Q99437,Q9UI12,Q9Y487,Q9Y5K8,	O15342,Q8NHE4,	O95670,Q96LB4,O75348,	P15313,P21281,	P21283,Q8NEY4,	Q8N8Y2,P61421,	Q96A05,P36543
  #       # #   CPX-6905	Vacuolar proton translocating ATPase complex, ATP6V0A3 variant	Translocates protons across a lipid bilayer via an ATP-driven rotary mechanism, thus acidifing the lumen of its resident organelle. Membrane-bound ion transporters/proton exchangers use the pH gradient to sequester metal ions to the vacuole and other cellular organelles. The combined action of the V-ATPase and membrane transporters plays a key role in maintaining cellular homoeostasis. The ATP6V0A3 variant localizes to endolysosomal compartments and can also be trafficked to the plasma membrane, being found at the ruffled border of osteoclasts.	O75787(1)|P27449(9)|P38606(3)|Q13488(1)|Q15904(1)|Q16864(1)|Q6P5S7(1)|Q99437(1)|Q9UI12(1)|Q9Y5K8(1)|[O15342,Q8NHE4](1)|[O95670,Q96LB4,O75348](3)|[P15313,P21281](3)|[P21283,Q8NEY4](1)|[Q8N8Y2,P61421](1)|[Q96A05,P36543](3)	-	Hetero 32-mer	O75787,P27449,P38606,Q13488,Q15904,Q16864,Q6P5S7,Q99437,Q9UI12,Q9Y5K8,	O15342,Q8NHE4,	O95670,Q96LB4,O75348,	P15313,P21281,	P21283,Q8NEY4,	Q8N8Y2,P61421,	Q96A05,P36543
  #       # #   CPX-6912	Vacuolar proton translocating ATPase complex, ATP6V0A4 variant	Translocates protons across a lipid bilayer via an ATP-driven rotary mechanism, thus acidifing the lumen of its resident organelle. Membrane-bound ion transporters/proton exchangers use the pH gradient to sequester metal ions to the vacuole and other cellular organelles. The combined action of the V-ATPase and membrane transporters plays a key role in maintaining cellular homoeostasis. The ATP6V0A4 variant is trafficked to the plasma membrane, being found at the apical membrane of type A intercalated cells.	O75787(1)|P27449(9)|P38606(3)|Q15904(1)|Q16864(1)|Q6P5S7(1)|Q99437(1)|Q9HBG4(1)|Q9UI12(1)|Q9Y5K8(1)|[O15342,Q8NHE4](1)|[O95670,Q96LB4,O75348](3)|[P15313,P21281](3)|[P21283,Q8NEY4](1)|[Q8N8Y2,P61421](1)|[Q96A05,P36543](3)	-	Hetero 32-mer	O75787,P27449,P38606,Q15904,Q16864,Q6P5S7,Q99437,Q9HBG4,Q9UI12,Q9Y5K8,	O15342,Q8NHE4,	O95670,Q96LB4,O75348,	P15313,P21281,	P21283,Q8NEY4,	Q8N8Y2,P61421,	Q96A05,P36543
  #       # #   CPX-5183	60S cytosolic large ribosomal subunit	Component of the ribosome, the site of protein biosynthesis resulting from translation of messenger RNA (mRNA). Responsible for the catalytic activity of the ribosome, the peptidyltransferase activity required to catalyze peptide bond formation. The nascent polypeptides leave the ribosome through a tunnel in the large subunit and interact with protein factors that function in enzymatic processing, targeting, and the membrane insertion of nascent chains at the exit of the ribosomal tunnel.	P05386(0)|P05387(0)|P05388(0)|P18077(1)|P18124(0)|P18621(1)|P26373(1)|P27635(1)|P30050(1)|P32969(0)|P36578(1)|P39023(1)|P40429(1)|P42766(1)|P46776(1)|P46777(1)|P46778(1)|P46779(1)|P47914(1)|P49207(1)|P50914(1)|P61313(1)|P61353(1)|P61513(1)|P61927(1)|P62424(0)|P62750(1)|P62829(1)|P62888(1)|P62891(1)|P62899(1)|P62906(1)|P62910(1)|P62913(1)|P62917(0)|P62945(1)|P62987(1)|P63173(1)|P83731(1)|P84098(1)|Q02543(1)|Q02878(1)|Q07020(1)|Q9Y3U8(1)|URS00000F9D45_9606(0)|URS000075D341_9606(0)|URS0000ABD82A_9606(0)|[P35268,Q6P5R6](1)|[Q969Q0,P83881](1)|[Q9UNX3,P61254](1)	-	Hetero 50-mer	P05386,P05387,P05388,P18077,P18124,P18621,P26373,P27635,P30050,P32969,P36578,P39023,P40429,P42766,P46776,P46777,P46778,P46779,P47914,P49207,P50914,P61313,P61353,P61513,P61927,P62424,P62750,P62829,P62888,P62891,P62899,P62906,P62910,P62913,P62917,P62945,P62987,P63173,P83731,P84098,Q02543,Q02878,Q07020,Q9Y3U8,	P35268,Q6P5R6,	Q969Q0,P83881,	Q9UNX3,P61254
  #       # #   CPX-7664	60S cytosolic large ribosomal subunit, testis-specific variant	Component of the ribosome, the site of protein biosynthesis resulting from translation of messenger RNA (mRNA). Responsible for the catalytic activity of the ribosome, the peptidyltransferase activity required to catalyze peptide bond formation. The nascent polypeptides leave the ribosome through a tunnel in the large subunit and interact with protein factors that function in enzymatic processing, targeting, and the membrane insertion of nascent chains at the exit of the ribosomal tunnel. This variant is found only in the testis and regulates the folding of a subset of male germ-cell-specific proteins that are essential for the formation of sperm.	P05386(0)|P05387(0)|P05388(0)|P18077(1)|P18124(0)|P18621(1)|P26373(1)|P30050(1)|P32969(0)|P36578(1)|P39023(1)|P40429(1)|P42766(1)|P46776(1)|P46777(1)|P46778(1)|P46779(1)|P47914(1)|P49207(1)|P50914(1)|P61313(1)|P61353(1)|P61513(1)|P61927(1)|P62424(0)|P62750(1)|P62829(1)|P62888(1)|P62899(1)|P62906(1)|P62910(1)|P62913(1)|P62917(0)|P62945(1)|P62987(1)|P63173(1)|P83731(1)|P84098(1)|Q02543(1)|Q02878(1)|Q07020(1)|Q96EH5(1)|Q96L21(1)|Q9Y3U8(1)|URS00000F9D45_9606(0)|URS000075D341_9606(0)|URS0000ABD82A_9606(0)|[P35268,Q6P5R6](1)|[Q969Q0,P83881](1)|[Q9UNX3,P61254](1)	-	Hetero 50-mer	P05386,P05387,P05388,P18077,P18124,P18621,P26373,P30050,P32969,P36578,P39023,P40429,P42766,P46776,P46777,P46778,P46779,P47914,P49207,P50914,P61313,P61353,P61513,P61927,P62424,P62750,P62829,P62888,P62899,P62906,P62910,P62913,P62917,P62945,P62987,P63173,P83731,P84098,Q02543,Q02878,Q07020,Q96EH5,Q96L21,Q9Y3U8,	P35268,Q6P5R6,	Q969Q0,P83881,	Q9UNX3,P61254
  #       # #   CPX-7665	60S cytosolic large ribosomal subunit, striated muscle variant	Component of the ribosome, the site of protein biosynthesis resulting from translation of messenger RNA (mRNA). Responsible for the catalytic activity of the ribosome, the peptidyltransferase activity required to catalyze peptide bond formation. The nascent polypeptides leave the ribosome through a tunnel in the large subunit and interact with protein factors that function in enzymatic processing, targeting, and the membrane insertion of nascent chains at the exit of the ribosomal tunnel. This variant is found only in the heart and skeletal muscle and regulates striated muscle function.	P05386(0)|P05387(0)|P05388(0)|P18077(1)|P18124(0)|P18621(1)|P26373(1)|P27635(1)|P30050(1)|P32969(0)|P36578(1)|P40429(1)|P42766(1)|P46776(1)|P46777(1)|P46778(1)|P46779(1)|P47914(1)|P49207(1)|P50914(1)|P61313(1)|P61353(1)|P61513(1)|P61927(1)|P62424(0)|P62750(1)|P62829(1)|P62888(1)|P62891(1)|P62899(1)|P62906(1)|P62910(1)|P62913(1)|P62917(0)|P62945(1)|P62987(1)|P63173(1)|P83731(1)|P83881(1)|P84098(1)|Q02543(1)|Q02878(1)|Q07020(1)|Q92901(1)|Q9Y3U8(1)|URS00000F9D45_9606(0)|URS000075D341_9606(0)|URS0000ABD82A_9606(0)|[P35268,Q6P5R6](1)|[Q9UNX3,P61254](1)	-	Hetero 50-mer	P05386,P05387,P05388,P18077,P18124,P18621,P26373,P27635,P30050,P32969,P36578,P40429,P42766,P46776,P46777,P46778,P46779,P47914,P49207,P50914,P61313,P61353,P61513,P61927,P62424,P62750,P62829,P62888,P62891,P62899,P62906,P62910,P62913,P62917,P62945,P62987,P63173,P83731,P83881,P84098,Q02543,Q02878,Q07020,Q92901,Q9Y3U8,	P35268,Q6P5R6,	Q9UNX3,P61254
  #       # #   CPX-2470	Vacuolar proton translocating ATPase complex, ATP6V0A1 variant	Translocates protons across a lipid bilayer via an ATP-driven rotary mechanism, thus acidifing the lumen of its resident organelle. Membrane-bound ion transporters/proton exchangers use the pH gradient to sequester metal ions to the vacuole and other cellular organelles. The combined action of the V-ATPase and membrane transporters plays a key role in maintaining cellular homoeostasis. The ATP6V0A1 variant localizes to synaptic vesicles in presynaptic neurons and to late endo/lysosomal compartments in non-neuronal tissues.	O75787(1)|P27449(9)|P38606(3)|Q15904(1)|Q16864(1)|Q6P5S7(1)|Q93050(1)|Q99437(1)|Q9UI12(1)|Q9Y5K8(1)|[O15342,Q8NHE4](1)|[O95670,Q96LB4,O75348](3)|[P15313,P21281](3)|[P21283,Q8NEY4](1)|[Q8N8Y2,P61421](1)|[Q96A05,P36543](3)	-	Hetero 32-mer	O75787,P27449,P38606,Q15904,Q16864,Q6P5S7,Q93050,Q99437,Q9UI12,Q9Y5K8,	O15342,Q8NHE4,	O95670,Q96LB4,O75348,	P15313,P21281,	P21283,Q8NEY4,	Q8N8Y2,P61421,	Q96A05,P36543
  # # combine manual combo complexes with rest of dataframe
  # complex_portal_expanded <- rbind(complex_portal_expanded %>% filter(!grepl("\\[", subunits_uniprot_ids_cp)),
  #       manual_rescue)
  # 
  # # only after removing non-proteins  count memberes
  # complex_portal_expanded <- complex_portal_expanded %>%
  #   mutate(complex_size = str_count(subunits_uniprot_ids_cp, ",") +1,   # each component separated by ","
  #          complex_counter_cp = 0, # Initialize complex_counter_cp column in complex_portal with zeros
  #          matches = "")        # Initialize matches column
  # 
  #  # get pro mapping (proteoforms are saved in complex portal - in order to match them with dataset translate to uniprot reviewed format)
  #  human_unreviewd_isoforms_mapping <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resource_uniprot/human_UNREVIEWED_isoform_mapping_20250306.tsv", header = TRUE) %>%
  #    # Extract all PRO_ values and collapse them into a single string per row
  #    mutate(across(c("chain", "propeptide", "initiator_methionine",
  #                    "post_translational_modification", "lipidation",
  #                    "cross_link", "disulfide_bond", "glycosylation",
  #                    "modified_residue", "peptide", "transit_peptide",
  #                    "signal_peptide"), # specify your columns here
  #                  ~ str_extract_all(., "PRO_\\d+", simplify = TRUE) %>%
  #                    apply(1, paste, collapse = " ") %>%
  #                    str_trim()
  #    )) %>%
  #    mutate(all_pro_values = paste(chain, propeptide, initiator_methionine,
  #                                  post_translational_modification, lipidation,
  #                                  cross_link, disulfide_bond, glycosylation,
  #                                  modified_residue, peptide, transit_peptide,
  #                                  signal_peptide, sep = " ")) %>%
  # 
  #    dplyr::select(entry, entry_name, reviewed, all_pro_values) %>%
  #    mutate(all_pro_values = str_squish(all_pro_values)) %>%  # clean up unnecessary spaces (e.g., trailing spaces or double spaces)
  #    filter(all_pro_values != "") %>%                         # remove empty rows
  #    separate_longer_delim(all_pro_values, delim = " ") # get each PRO_ in to separate colum (sometimes multiple isoforms per protein --> to match them all assign individually (row wise) to entry)
  # 
  #  # Create a named vector for mapping PRO IDs to UniProt IDs
  #  pro_to_uniprot <- setNames(human_unreviewd_isoforms_mapping$entry, human_unreviewd_isoforms_mapping$all_pro_values)
  #  # Function to replace PRO IDs with UniProt IDs
  #  replace_pro_ids <- function(x) {
  #    ids <- unlist(strsplit(x, ","))
  #    replaced_ids <- ifelse(grepl("^PRO_", ids), pro_to_uniprot[ids], ids)
  #    paste(replaced_ids, collapse = ",")
  #  }
  #  # Apply the replacement function to the subunits_uniprot_ids_cp column
  #  complex_portal_expanded$subunits_uniprot_ids_cp <- sapply(complex_portal_expanded$subunits_uniprot_ids_cp, replace_pro_ids)
  # 
  #  # # remove isoform information
  #  complex_portal_expanded <- complex_portal_expanded %>%
  #    mutate(subunits_uniprot_ids_cp = str_replace_all(subunits_uniprot_ids_cp, "-\\d+", "")) %>%
  #    distinct() %>%
  #    rename("interaction_identifiers" = "number_complex_ac") %>%
  #    mutate(interaction_identifiers = gsub("CPX-", "complexportal_", interaction_identifiers))
  # 
  #  write.csv(complex_portal_expanded, "/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_human_MG20250305.csv", row.names = FALSE)
  #  complex_portalMG <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_human_MG20250305.csv")
  # 
  #  # complex portal ppi extraction ---------------------------------------------
  #  complex_portal_MG_ppi <- complex_portalMG %>%
  #    dplyr::select(interaction_identifiers, subunits_uniprot_ids_cp, complex_size) %>%
  #    # get complex format (e.g. Q92828,Q13227,O15379,O75376,O60907,Q9BZK7) to ppi format interactor_a - interactor_b
  #    rename("subunits_uniprot_id" = "subunits_uniprot_ids_cp") %>%
  #    separate_rows(subunits_uniprot_id, sep = ",") %>%
  #    group_by(interaction_identifiers) %>%
  #    mutate(
  #      interactor_a = first(subunits_uniprot_id),
  #      interactor_b = subunits_uniprot_id) %>%
  #    filter(interactor_a != interactor_b) %>%
  #    ungroup() %>%
  #    dplyr::select(interaction_identifiers, interactor_a, interactor_b)
  #  write.csv(complex_portal_MG_ppi, "/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_ppi_MG20250306.csv", row.names = FALSE)
  #  complex_portal_MG_ppi <- read_protti("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_ppi_MG20250306.csv")
  #  ##########################################################################################################################################################################

   # reimport relevant complex portal format
   complex_portalMG <- read.csv("/Users/mgesell/Desktop/currentR/git/shs_resources/resources_complex/_complex_portal_human_MG20250305.csv", head=TRUE)
   # Iterate through each entry in query_list
   for (entry in query_list) {
     # Use grepl to find matches in complex_portalMG$subunits_uniprot_ids_cp
     matches <- grepl(entry, complex_portalMG$subunits_uniprot_ids_cp, fixed = TRUE)
     # Increment complex_counter_cp for matching rows
     complex_portalMG$complex_counter_cp[matches] <- complex_portalMG$complex_counter_cp[matches] + 1
     # Append the matched entry to the 'matches' column
     complex_portalMG$matches[matches] <- ifelse(is.na(complex_portalMG$matches[matches]) | complex_portalMG$matches[matches] == "",
                                                 entry,
                                                 paste(complex_portalMG$matches[matches], entry, sep = ", "))
   }

  
  complex_portal_query_subset <- complex_portalMG %>%
    mutate(subunits_uniprot_ids_cp = gsub("\\([0-9]+\\)", "", identifiers_and_stoichiometry_of_molecules_in_complex),
           subunits_uniprot_ids_cp = gsub("\\|", ", ", subunits_uniprot_ids_cp)) %>%
    dplyr::select(interaction_identifiers, description, recommended_name, subunits_uniprot_ids_cp, matches, complex_size, complex_counter_cp) %>%
    filter(complex_counter_cp >= 2) %>%                          # single protein is no complex ...
    mutate(complex_recall_cp = complex_counter_cp/complex_size) %>% # calculate recall
    arrange(desc(complex_recall_cp))  %>%
    mutate(comparison = paste(set),
           subunits_uniprot_ids_cp = gsub(";", ", ", subunits_uniprot_ids_cp),
           subunits_uniprot_ids_cp = sapply(strsplit(subunits_uniprot_ids_cp, ", "), function(x) paste(sort(x), collapse = ", " ))) %>%  # alphabetical display of complex units (easier to compare with matches when ordered)
    rename("complex_name_cp" = "interaction_identifiers", "description_cp" = "description")
  
  return(complex_portal_query_subset)
  
}
