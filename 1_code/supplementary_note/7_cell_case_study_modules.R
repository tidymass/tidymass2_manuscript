library(r4projects)
setwd(get_project_wd())

library(metpath)

load("3_data_analysis/6_supplementary_note/fpa_input/feature_table_marker_all.rda")
load("3_data_analysis/6_supplementary_note/fpa_result.rda")
load("3_data_analysis/6_supplementary_note/fpa_input/metabolic_network.rda")

fpa_result <- result
feature_table_marker <- feature_table_marker_all

for (i in 1:17) {
  cat(i, " ")
  
  plot <- plot_metabolic_module_fpa(
    fpa_result = fpa_result,
    feature_table_marker = feature_table_marker,
    include_feature = TRUE,
    include_hidden_metabolites = FALSE,
    add_compound_name = TRUE,
    metabolic_module_index = i,
    layout = "fr",
    add_pathways = TRUE
  )
  
  ggsave(
    plot,
    filename = paste0("4_manuscript/Supplementary_note/cell_case_study/pathways/module_", i, ".pdf"),
    width = 14,
    height = 8
  )
  
  module_name <- paste0("metabolic_module_", i)
  
  enriched_data <- fpa_result$enriched_pathways_list[[module_name]][, c("pathway_id", "pathway_name", "p_value_adjust", "mapped_id")] %>%
    dplyr::filter(p_value_adjust < 0.05)
  
  write.csv(
    enriched_data,
    file = paste0("4_manuscript/Supplementary_note/cell_case_study/pathways/metabolic_module_", i, ".csv"),
    row.names = FALSE
  )
}