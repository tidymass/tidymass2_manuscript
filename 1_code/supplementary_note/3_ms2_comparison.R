library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(ggplot2)
library(tidyverse)

load("3_data_analysis/6_supplementary_note/object_in_house_annotation.rda")
load("2_data/ms2_original/mpsnyder_rplc_ms2.rda")
load("3_data_analysis/6_supplementary_note/fpa_validation.rda")

urine_annotation <- object_new@annotation_table

ms2_data <- object_new@ms2_data[[1]]

# for M363T580_POS
which(ms2_data@variable_id == "M363T580_POS")

spec_1 <- ms2_data@ms2_spectra[[3160]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_679"]][["NCE25"]]

p <- masstools::ms2_plot(spec_1, 
                              spec_1_standard,
                              spectrum1_name = "M363T580_POS",
                              spectrum2_name = "Cortisol")
                              

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M363T580_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)




# for M311T602_POS
which(ms2_data@variable_id == "M311T602_POS")

spec_1 <- ms2_data@ms2_spectra[[2603]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_403"]][["NCE50"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M311T602_POS",
                         spectrum2_name = "L-Octanoylcarnitine")

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE50/M311T602_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)


# for M274T390_POS
which(ms2_data@variable_id == "M274T390_POS")

spec_1 <- ms2_data@ms2_spectra[[2081]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_259"]][["NCE25"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M274T390_POS",
                         spectrum2_name = "N-acetylneuraminate")

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M274T390_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)


# for M292T390_POS
which(ms2_data@variable_id == "M292T390_POS")

spec_1 <- ms2_data@ms2_spectra[[2325]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_259"]][["NCE25"]]


p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M292T390_POS",
                         spectrum2_name = "N-acetylneuraminate")

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M292T390_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)


# for M232T328_2_POS
which(ms2_data@variable_id == "M232T328_2_POS")

spec_1 <- ms2_data@ms2_spectra[[1526]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_399"]][["NCE25"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M232T328_2_POS",
                         spectrum2_name = "O-Butanoylcarnitine")

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M232T328_2_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)

