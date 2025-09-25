library(r4projects)
setwd(get_project_wd())
rm(list = ls())

load("3_data_analysis/4_smartd_project_fpa/fpa_result.rda")
load("3_data_analysis/6_supplementary_note/object_in_house_annotation.rda")

library(tidyverse)

# unique id for fpa annotation results
length(unique(fpa_result$annotation_table$variable_id))

fpa_annotation_table <- fpa_result$annotation_table

fpa_id <- fpa_annotation_table$variable_id %>% unique()

# matched id for fpa and urine data annotation table
urine_annotation <- object_new@annotation_table %>% dplyr::filter(!is.na(KEGG.ID))

# urine_annotation <- urine_annotation %>% 
#   dplyr::filter(Level == "1")
# urine_annotation2 <- urine_metabolomics_data@annotation_table %>% dplyr::filter(!is.na(SS))
urine_id <- urine_annotation$variable_id %>% unique()

matched_id <- intersect(fpa_id, urine_id)

# fpa annotation table for matched id
fpa_matched <- fpa_annotation_table %>% dplyr::filter(variable_id %in% matched_id)

# urine annotation table for matched id
urine_matched <- urine_annotation %>% dplyr::filter(variable_id %in% matched_id)

urine_long <- urine_matched %>%
  dplyr::select(variable_id, CAS.ID, HMDB.ID, KEGG.ID) %>%
  pivot_longer(
    cols = c(CAS.ID, HMDB.ID, KEGG.ID),
    names_to = "id_type",
    values_to = "id"
  ) %>%
  dplyr::filter(!is.na(id)) %>%
  distinct(variable_id, id_type, id)

fpa_long <- fpa_matched %>%
  mutate(.row_id = row_number()) %>%
  dplyr::select(.row_id, variable_id, CAS.ID, HMDB.ID, KEGG.ID) %>%
  pivot_longer(
    cols = c(CAS.ID, HMDB.ID, KEGG.ID),
    names_to = "id_type",
    values_to = "id"
  ) %>%
  dplyr::filter(!is.na(id))

hit_rows <- fpa_long %>%
  inner_join(urine_long, by = c("variable_id", "id_type", "id")) %>%
  distinct(.row_id)

fpa_validation <- fpa_matched %>%
  mutate(.row_id = row_number()) %>%
  mutate(matched_result = if_else(.row_id %in% hit_rows$.row_id, "yes", "no")) %>%
  dplyr::select(-.row_id)

# length(unique(fpa_validation %>% dplyr::filter(matched_result == "yes") %>% .$variable_id))
hit_id <- unique(fpa_validation %>% dplyr::filter(matched_result == "yes") %>% .$variable_id)

save(fpa_validation, file = "3_data_analysis/6_supplementary_note/fpa_validation.rda")




