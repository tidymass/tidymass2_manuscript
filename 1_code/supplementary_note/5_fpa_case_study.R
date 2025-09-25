library(r4projects)
setwd(get_project_wd())
rm(list = ls())

load("2_data/feature_analysis_module/importance_down_metabolite.rda")
load("2_data/feature_analysis_module/importance_up_metabolite.rda")
load("2_data/peaks/variable_info")

library(tidyverse)
library(metpath)


feature_up_clean <- 
  importance_up_metabolite %>% 
  dplyr::select(variable_id, mz, rt, score, p1) %>% 
  dplyr::rename(
    degree = score,
    p_value = p1
  )

feature_down_clean <-
  importance_down_metabolite %>% 
  dplyr::select(variable_id, mz, rt, score, p1) %>% 
  dplyr::rename(
    degree = score,
    p_value = p1
  )

feature_table_marker_all <- 
  rbind(
    feature_up_clean,
    feature_down_clean
  ) %>%
  mutate(
    polarity = ifelse(grepl("_POS$", variable_id), "positive", "negative")
  )

save(feature_table_marker_all, file = "3_data_analysis/6_supplementary_note/fpa_input/feature_table_marker_all.rda")

feature_table_all_clean <- 
  variable_info %>% 
  dplyr::select(variable_id, mz, rt) %>%
  mutate(
    polarity = ifelse(grepl("_POS$", variable_id), "positive", "negative")
  ) 

load("3_data_analysis/6_supplementary_note/fpa_input/kegg_compound_database.rda")
load("3_data_analysis/6_supplementary_note/fpa_input/kegg_hsa_pathway.rda")
load("3_data_analysis/6_supplementary_note/fpa_input/metabolic_network.rda")

fpa_result <- 
  perform_fpa(
    feature_table_marker = feature_table_marker_all,
    feature_table_all = feature_table_all_clean,
    metabolite_database = kegg_compound_database,
    column = "rp",
    adduct.table = NULL,
    ms1.match.ppm = 25,
    rt.match.tol = 5,
    ##unit is second
    mz.ppm.thr = 400,
    threads = 3,
    include_hidden_metabolites = FALSE,
    metabolic_network = metabolic_network,
    pathway_database = kegg_hsa_pathway
  )

save(fpa_result, file = "3_data_analysis/6_supplementary_note/fpa_result.rda")

