library(r4projects)
setwd(get_project_wd())
# source("R/enrich_pathways.R")
rm(list = ls())

load("2_data/urine_pregnancy/urine_metabolomics_data.rda")
load("2_data/urine_pregnancy/peak_marker")

load("2_data/urine_pregnancy/peak_marker")
load("2_data/urine_pregnancy/urine_metabolomics_data.rda")
load("2_data/urine_pregnancy/kegg_compound_ms1.rda")

load("2_data/urine_pregnancy/metabolic_network.rda")
load("2_data/urine_pregnancy/kegg_hsa_pathway.rda")

load("2_data/ms1_database.rda")

library(tidymass)
library(tidyverse)
library(igraph)
library(ggraph)

up_marker <-
  peak_marker %>%
  dplyr::filter(score > 0 & correlation1 > 0.3 & p1_adj < 0.05)  %>%
  dplyr::select("Gene ID", score, p1_adj)

down_marker <-
  peak_marker %>%
  dplyr::filter(score < 0 &
                  correlation1 < -0.3 & p1_adj < 0.05)  %>%
  dplyr::select("Gene ID", score, p1_adj)

up_marker <-
  up_marker %>%
  dplyr::left_join(urine_metabolomics_data@variable_info,
                   by = c("Gene ID" = "variable_id"))

down_marker <-
  down_marker %>%
  dplyr::left_join(urine_metabolomics_data@variable_info,
                   by = c("Gene ID" = "variable_id"))

feature_table_marker_up <-
  up_marker %>%
  dplyr::select("Gene ID", mz, rt, score, p1_adj) %>%
  dplyr::rename(variable_id = "Gene ID",
                degree = score,
                p_value = p1_adj) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  ))

feature_table_marker_down <-
  down_marker %>%
  dplyr::select("Gene ID", mz, rt, score, p1_adj) %>%
  dplyr::rename(variable_id = "Gene ID",
                degree = score,
                p_value = p1_adj) %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  ))

feature_table_marker <-
  rbind(feature_table_marker_up, feature_table_marker_down)

feature_table_all <-
  urine_metabolomics_data@variable_info %>%
  dplyr::mutate(polarity = case_when(
    stringr::str_detect(variable_id, "POS") ~ "positive",
    stringr::str_detect(variable_id, "NEG") ~ "negative"
  )) %>%
  dplyr::select(variable_id, mz, rt, polarity)


dir.create("3_data_analysis/1_smartd_project_data_preparation")
setwd("3_data_analysis/1_smartd_project_data_preparation")

annotation_table <-
  urine_metabolomics_data@variable_info

annotation_table_hmdb <-
  annotation_table %>%
  dplyr::filter(!is.na(HMDB.ID))

annotation_table_kegg <-
  annotation_table %>%
  dplyr::filter(!is.na(KEGG.ID)) %>%
  dplyr::filter(!Lab.ID %in% annotation_table_hmdb$Lab.ID)

annotation_table_cas <-
  annotation_table %>%
  dplyr::filter(!is.na(CAS.ID)) %>%
  dplyr::filter(!Lab.ID %in% annotation_table_hmdb$Lab.ID) %>%
  dplyr::filter(!Lab.ID %in% annotation_table_kegg$Lab.ID)

annotation_table_hmdb <-
  annotation_table_hmdb %>%
  dplyr::distinct(HMDB.ID, .keep_all = TRUE) %>%
  dplyr::left_join(
    ms1_database %>%
      dplyr::filter(!is.na(HMDB_ID)) %>%
      dplyr::distinct(HMDB_ID, .keep_all = TRUE) %>%
      dplyr::select(c(HMDB_ID, from_human:from_which_food)),
    by = c("HMDB.ID" = "HMDB_ID")
  )

annotation_table_kegg <-
  annotation_table_kegg %>%
  dplyr::distinct(KEGG.ID, .keep_all = TRUE) %>%
  dplyr::left_join(
    ms1_database %>%
      dplyr::filter(!is.na(KEGG_ID)) %>%
      dplyr::distinct(KEGG_ID, .keep_all = TRUE) %>%
      dplyr::select(c(KEGG_ID, from_human:from_which_food)),
    by = c("KEGG.ID" = "KEGG_ID")
  )

annotation_table_cas <-
  annotation_table_cas %>%
  dplyr::distinct(CAS.ID, .keep_all = TRUE) %>%
  dplyr::left_join(
    ms1_database %>%
      dplyr::filter(!is.na(CAS_ID)) %>%
      dplyr::distinct(CAS_ID, .keep_all = TRUE) %>%
      dplyr::select(c(CAS_ID, from_human:from_which_food)),
    by = c("CAS.ID" = "CAS_ID")
  )

dim(annotation_table_hmdb)
dim(annotation_table_kegg)
dim(annotation_table_cas)

nrow(annotation_table_hmdb) +
  nrow(annotation_table_kegg) +
  nrow(annotation_table_cas)

sum(!is.na(annotation_table$Compound.name))

annotation_table <-
  rbind(annotation_table_hmdb, annotation_table_kegg, annotation_table_cas)

urine_metabolomics_data@annotation_table <-
  annotation_table

save(feature_table_marker_up, file = "feature_table_marker_up.rda")
save(feature_table_marker_down, file = "feature_table_marker_down.rda")
save(feature_table_marker, file = "feature_table_marker.rda")
save(feature_table_all, file = "feature_table_all.rda")
save(kegg_hsa_pathway, file = "kegg_hsa_pathway.rda")
save(metabolic_network, file = "metabolic_network.rda")
save(kegg_compound_database, file = "kegg_compound_database.rda")
save(urine_metabolomics_data, file = "urine_metabolomics_data.rda")
