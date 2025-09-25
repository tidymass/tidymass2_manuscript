library(r4projects)
setwd(get_project_wd())

library(metpath)

load("3_data_analysis/6_supplementary_note/fpa_input/feature_table_marker_all.rda")
load("3_data_analysis/6_supplementary_note/fpa_result.rda")
load("3_data_analysis/6_supplementary_note/fpa_input/metabolic_network.rda")

fpa_result <- result
feature_table_marker <- feature_table_marker_all

###all enriched pathways
plot <-
  enrich_scatter_plot(fpa_result$enriched_pathways)


write.csv(fpa_result$enriched_pathways@result[,c("pathway_id", "pathway_name", "p_value_adjust", "mapped_id")] %>% 
            dplyr::filter(p_value_adjust < 0.05), 
          file = "4_manuscript/Supplementary_note/cell_case_study/pathways/pathway_overall.csv", row.names = FALSE)

plot
ggsave(plot,
       filename = "4_manuscript/Supplementary_note/cell_case_study/module_plots/all_pathways_enrich_scatter_plot.pdf",
       width = 8,
       height = 5)


plot_without_feature <-
  plot_metabolic_network_fpa(
    fpa_result = fpa_result,
    feature_table_marker = feature_table_marker,
    include_feature = FALSE,
    node_color_by_module = TRUE,
    include_hidden_metabolites = FALSE,
    add_compound_name = TRUE,
    layout = "fr"
  )

plot_without_feature


ggsave(
  plot_without_feature,
  filename = "4_manuscript/Supplementary_note/cell_case_study/module_plots/plot_without_feature.pdf",
  width = 14,
  height = 10,
  device = cairo_pdf
)

