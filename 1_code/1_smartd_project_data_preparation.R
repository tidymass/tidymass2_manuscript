# library(r4projects)
# setwd(get_project_wd())
# # source("R/enrich_pathways.R")
# 
# load("demo_data/2_smartd_project/urine_metabolomics_data.rda")
# load("demo_data/2_smartd_project/peak_marker")
# 
# load("demo_data/2_smartd_project/peak_marker")
# load("demo_data/2_smartd_project/urine_metabolomics_data.rda")
# load("demo_data/2_smartd_project/kegg_compound_ms1.rda")
# 
# load("demo_data/2_smartd_project/metabolic_network.rda")
# load("demo_data/2_smartd_project/kegg_hsa_pathway.rda")
# 
# library(tidymass)
# library(tidyverse)
# library(igraph)
# library(ggraph)
# 
# up_marker <-
#   peak_marker %>%
#   dplyr::filter(score > 0 & correlation1 > 0.3 & p1_adj < 0.05)  %>%
#   dplyr::select("Gene ID", score, p1_adj)
# 
# down_marker <-
#   peak_marker %>%
#   dplyr::filter(score < 0 &
#                   correlation1 < -0.3 & p1_adj < 0.05)  %>%
#   dplyr::select("Gene ID", score, p1_adj)
# 
# up_marker <-
#   up_marker %>%
#   dplyr::left_join(urine_metabolomics_data@variable_info,
#                    by = c("Gene ID" = "variable_id"))
# 
# down_marker <-
#   down_marker %>%
#   dplyr::left_join(urine_metabolomics_data@variable_info,
#                    by = c("Gene ID" = "variable_id"))
# 
# feature_table_marker_up <-
#   up_marker %>%
#   dplyr::select("Gene ID", mz, rt, score, p1_adj) %>%
#   dplyr::rename(variable_id = "Gene ID",
#                 degree = score,
#                 p_value = p1_adj) %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(variable_id, "POS") ~ "positive",
#     stringr::str_detect(variable_id, "NEG") ~ "negative"
#   ))
# 
# feature_table_marker_down <-
#   down_marker %>%
#   dplyr::select("Gene ID", mz, rt, score, p1_adj) %>%
#   dplyr::rename(variable_id = "Gene ID",
#                 degree = score,
#                 p_value = p1_adj) %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(variable_id, "POS") ~ "positive",
#     stringr::str_detect(variable_id, "NEG") ~ "negative"
#   ))
# 
# feature_table_marker <-
#   rbind(feature_table_marker_up, feature_table_marker_down)
# 
# feature_table_all <-
#   urine_metabolomics_data@variable_info %>%
#   dplyr::mutate(polarity = case_when(
#     stringr::str_detect(variable_id, "POS") ~ "positive",
#     stringr::str_detect(variable_id, "NEG") ~ "negative"
#   )) %>%
#   dplyr::select(variable_id, mz, rt, polarity)
# 
# 
# dir.create("demo_data/1_manuscript/3_data_analysis/1_smartd_project_data_preparation")
# setwd("demo_data/1_manuscript/3_data_analysis/1_smartd_project_data_preparation")
# 
# save(feature_table_marker_up, file = "feature_table_marker_up.rda")
# save(feature_table_marker_down, file = "feature_table_marker_down.rda")
# save(feature_table_marker, file = "feature_table_marker.rda")
# save(feature_table_all, file = "feature_table_all.rda")
# save(urine_metabolomics_data, file = "urine_metabolomics_data.rda")
# save(kegg_hsa_pathway, file = "kegg_hsa_pathway.rda")
# save(metabolic_network, file = "metabolic_network.rda")
# save(kegg_compound_database, file = "kegg_compound_database.rda")
# 
