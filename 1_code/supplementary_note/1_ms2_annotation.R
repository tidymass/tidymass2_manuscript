library(r4projects)
setwd(get_project_wd())
rm(list = ls())

load("3_data_analysis/1_smartd_project_data_preparation/urine_metabolomics_data.rda")

library(tidymass)

old_object <- urine_metabolomics_data

old_object@annotation_table <- data.frame()
# old_object@annotation_table <- old_object@annotation_table[0, ]

variable_info <- old_object@variable_info

variable_info <- variable_info %>% 
  dplyr::select(variable_id, mz, rt)

variable_info_pos <- variable_info[grepl("_POS$", variable_info$variable_id), ]
variable_info_neg <- variable_info[grepl("_NEG$", variable_info$variable_id), ]

object_pos <- old_object
object_neg <- old_object

object_pos@variable_info <- variable_info_pos
object_neg@variable_info <- variable_info_neg

# add ms2 to the objects
object_pos2 =
  mutate_ms2(
    object = object_pos,
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/MS2/POS/mgf"
  )

object_neg2 =
  mutate_ms2(
    object = object_neg,
    column = "rp",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/MS2/NEG/mgf"
  )

# load the in house database
load("2_data/ms2_original/mpsnyder_rplc_ms2.rda")

# positive mode
object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    rt.match.tol = 30,
    polarity = "positive",
    column = "rp",
    database = gnps_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

# negative mode
object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    rt.match.tol = 30,
    polarity = "negative",
    column = "rp",
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )


object_new <- object_pos2
annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table
annotation_table <- rbind(annotation_table_pos, annotation_table_neg)

object_new@annotation_table <- annotation_table

variable_info_pos <- object_pos2@variable_info
variable_info_neg <- object_neg2@variable_info
variable_info <- rbind(variable_info_pos, variable_info_neg)
object_new@variable_info <- variable_info

save(object_new, file = "3_data_analysis/6_supplementary_note/object_in_house_annotation.rda")




