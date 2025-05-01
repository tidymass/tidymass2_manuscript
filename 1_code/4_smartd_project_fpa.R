library(r4projects)
setwd(get_project_wd())
# source("R/1_metpath_logo.R")
# source("R/7_utils.R")
# source("R/8_zzz.R")
# source("R/9_hmdb_database.R")
# source("R/10_kegg_database.R")
# source("R/11_feature_based_pathway_enrichment.R")
# source("R/12_enrich_hmdb.R")
# source("R/14_enrich_pathways.R")
# source("R/15_enrich_result_class.R")
# source("R/16_filter_pathways.R")
# source("R/17_pathway_database-class.R")
# source("R/18_dplyr-enrich_result.R")
# source("R/19_dplyr-pathway_database.R")
# source("R/20_pathway_database-methods.R")
# source("R/21_graphics.R")

load("3_data_analysis/1_smartd_project_data_preparation/feature_table_all.rda")
load("3_data_analysis/1_smartd_project_data_preparation/feature_table_marker.rda")
load("3_data_analysis/1_smartd_project_data_preparation/metabolic_network.rda")



load(
  "3_data_analysis/1_smartd_project_data_preparation/urine_metabolomics_data.rda"
)
load("3_data_analysis/1_smartd_project_data_preparation/kegg_hsa_pathway.rda")
load("3_data_analysis/1_smartd_project_data_preparation/kegg_compound_database.rda")

dir.create("3_data_analysis/4_smartd_project_fpa", showWarnings = FALSE)
setwd("3_data_analysis/4_smartd_project_fpa")

library(tidygraph)
library(ggraph)

temp <-
  metabolic_network %>%
  activate(nodes) %>%
  as_tibble() %>%
  pull(degree)

sum(temp <= 10) / length(temp)



library(metid)

# fpa_result <-
#   perform_fpa(
#     feature_table_marker = feature_table_marker,
#     feature_table_all = feature_table_all,
#     metabolite_database = metabolite_database,
#     column = column,
#     adduct.table = NULL,
#     threads = 8,
#     include_hidden_metabolites = FALSE,
#     metabolic_network = metabolic_network,
#     pathway_database = kegg_hsa_pathway
#   )
#
# save(fpa_result, file = "fpa_result.rda")
load("fpa_result.rda")

write.csv(fpa_result$dysregulated_metabolic_module[,c("Name", "p_value", "Total_metabolite_number", "Total_metabolite_id")] %>% 
            dplyr::filter(p_value < 0.05), 
          file = "fpa_result_dysregulated_metabolic_module.csv")

length(unique(fpa_result$annotation_table$variable_id)) / 1044

library(metpath)

fpa_result$enriched_pathways@result %>% dplyr::filter(p_value_adjust < 0.05) %>%
  pull(pathway_name)

plot_with_feature <-
  plot_metabolic_network_fpa(
    fpa_result = fpa_result,
    feature_table_marker = feature_table_marker,
    include_feature = TRUE,
    node_color_by_module = FALSE,
    include_hidden_metabolites = FALSE,
    add_compound_name = TRUE
  )

plot_with_feature

library(extrafont)
loadfonts()

# ggsave(
#   plot_with_feature,
#   filename = "plot_with_feature.pdf",
#   width = 10,
#   height = 10
# )

dim(feature_table_marker)

sum(feature_table_marker$variable_id %in% fpa_result$annotation_table$variable_id)

marker_with_ms2_annotation <-
  feature_table_marker %>%
  dplyr::left_join(urine_metabolomics_data@variable_info[, c("variable_id",
                                                             "KEGG.ID",
                                                             "Compound.name",
                                                             "Level",
                                                             "HMDB.ID",
                                                             "Database")], by = "variable_id") %>%
  dplyr::filter(!is.na(Level))


sum(
  marker_with_ms2_annotation$variable_id %in% fpa_result$annotation_table$variable_id
)

marker_with_ms2_annotation[which(
  !marker_with_ms2_annotation$variable_id %in% fpa_result$annotation_table$variable_id
), ]

marker_with_ms2_annotation %>%
  dplyr::left_join(fpa_result$annotation_table, by = "variable_id") %>%
  dplyr::arrange(variable_id) %>%
  dplyr::select(variable_id,
                Compound.name.x,
                Compound.name.y,
                KEGG.ID.x,
                KEGG.ID.y)

###all enriched pathways
plot <-
  enrich_scatter_plot(fpa_result$enriched_pathways)


write.csv(fpa_result$enriched_pathways@result[,c("pathway_id", "pathway_name", "p_value_adjust", "mapped_id")] %>% 
            dplyr::filter(p_value_adjust < 0.05), 
          file = "fpa_result_enriched_pathways.csv")

plot
# ggsave(plot,
#        filename = "all_pathways_enrich_scatter_plot.pdf",
#        width = 8,
#        height = 5)


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

library(extrafont)
loadfonts()

# ggsave(
#   plot_without_feature,
#   filename = "plot_without_feature.pdf",
#   width = 14,
#   height = 10
# )

plot_metabolic_module_fpa(
  fpa_result = fpa_result,
  feature_table_marker = feature_table_marker,
  include_feature = TRUE,
  include_hidden_metabolites = FALSE,
  add_compound_name = TRUE,
  metabolic_module_index = 4,
  layout = "fr",
  add_pathways = TRUE
)

####quantitative of modules

annotation_table <-
  fpa_result$annotation_table
# dplyr::filter(isotope == "[M]")

library(plyr)
library(tidymass)

expression_data <-
  urine_metabolomics_data %>%
  scale_data(center = FALSE) %>%
  extract_expression_data()

module_quantitative_data <-
  fpa_result$dysregulated_metabolic_module$Total_metabolite_id[fpa_result$dysregulated_metabolic_module$p_value < 0.05] %>%
  purrr::map(function(x) {
    KEGG_ID <-
      stringr::str_split(x, pattern = "\\{\\}") %>%
      unlist()
    
    temp_data <-
      data.frame(KEGG.ID = KEGG_ID) %>%
      dplyr::left_join(annotation_table[, c("variable_id",
                                            "KEGG.ID",
                                            "compound_class",
                                            "Adduct",
                                            "isotope",
                                            "score")], by = "KEGG.ID") %>%
      dplyr::arrange(KEGG.ID) %>%
      dplyr::distinct(KEGG.ID, variable_id, compound_class, .keep_all = TRUE) %>%
      dplyr::filter(!is.na(variable_id)) %>%
      plyr::dlply(.variables = .(KEGG.ID)) %>%
      purrr::map(function(y) {
        if (max(y$score) >= 80) {
          y <-
            y %>%
            dplyr::filter(score >= 80)
        }
        
        remain_compound_class <-
          y %>%
          dplyr::count(compound_class) %>%
          dplyr::filter(n > 1) %>%
          pull(compound_class)
        
        
        if (length(remain_compound_class) > 0) {
          y <-
            y %>% dplyr::filter(compound_class %in% remain_compound_class)
        }
        y %>%
          dplyr::filter(isotope == "[M]")
      })
    
    remain_idx <-
      temp_data %>%
      lapply(nrow) %>%
      unlist() %>%
      `>`(0) %>%
      which()
    
    temp_data <-
      temp_data[remain_idx]
    
    sam_score <-
      temp_data %>%
      lapply(function(x) {
        feature_table_marker %>%
          dplyr::filter(variable_id %in% x$variable_id) %>%
          pull(degree) %>%
          mean()
      }) %>%
      unlist
    
    p_value <-
      temp_data %>%
      lapply(function(x) {
        feature_table_marker %>%
          dplyr::filter(variable_id %in% x$variable_id) %>%
          pull(p_value) %>%
          min()
      }) %>%
      unlist
    
    data.frame(sam_score = sam_score, p_value = p_value)  %>%
      tibble::rownames_to_column(var = "KEGG.ID")
    
  })

module_quantitative_data <-
  1:length(module_quantitative_data) %>%
  purrr::map(function(i) {
    cat(i, " ")
    data.frame(module_name = fpa_result$dysregulated_metabolic_module$Name[i],
               module_quantitative_data[[i]])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

up_down_color <-
  c(
    increase = RColorBrewer::brewer.pal(name = "RdYlBu", n = 5)[1],
    decrease = RColorBrewer::brewer.pal(name = "RdYlBu", n = 5)[5],
    non = "grey"
  )

plot1 <-
  module_quantitative_data %>%
  dplyr::mutate(module_name = factor(module_name, levels = stringr::str_sort(unique(module_name), numeric = TRUE))) %>%
  dplyr::mutate(change = case_when(sam_score > 0 ~ "increase", sam_score < 0 ~ "decrease")) %>%
  ggplot(aes(module_name, sam_score)) +
  geom_point(aes(size = -log(p_value, 10), fill = change), shape = 21) +
  scale_fill_manual(values = up_down_color) +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = "")

plot1

temp_data <-
  module_quantitative_data %>%
  plyr::dlply(.variables = .(module_name)) %>%
  purrr::map(function(x) {
    data.frame(increase = sum(x$sam_score > 0),
               decrease = sum(x$sam_score < 0))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data$non <-
  temp_data %>%
  apply(1, function(x) {
    max(temp_data$increase + temp_data$decrease) - sum(x)
  })


library(tibble)

plot2 <-
  temp_data %>%
  rownames_to_column(var = "module_name") %>%
  pivot_longer(
    cols = c(increase, decrease, non),
    names_to = "change",
    values_to = "number"
  ) %>%
  dplyr::mutate(change = factor(change, rev(c(
    "increase", "decrease", "non"
  )))) %>%
  ggplot(aes(module_name, number, fill = change)) +
  geom_bar(stat = "identity",
           position = "stack",
           color = "black") +
  scale_fill_manual(values = up_down_color) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    panel.grid = element_blank()
  ) +
  labs(x = "")

plot2

library(tibble)

plot3 <-
  module_quantitative_data %>%
  plyr::dlply(.variables = .(module_name)) %>%
  purrr::map(function(x) {
    sum(x$sam_score) / nrow(x)
  }) %>%
  unlist() %>%
  data.frame(quantitative_score = .) %>%
  tibble::rownames_to_column(var = "module_name") %>%
  dplyr::mutate(module_name = factor(module_name, levels = stringr::str_sort(unique(module_name), numeric = TRUE))) %>%
  dplyr::mutate(change = case_when(
    quantitative_score > 0 ~ "increase",
    quantitative_score < 0 ~ "decrease"
  )) %>%
  ggplot(aes(module_name, quantitative_score)) +
  geom_line(group = 1) +
  geom_point(size = 5, aes(fill = change), shape = 21) +
  scale_fill_manual(values = up_down_color) +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  ) +
  labs(x = "")

plot3

library(patchwork)

plot <-
  plot2 + plot3 + plot1 + plot_layout(ncol = 1, heights = c(1, 1, 2))

# ggsave(plot,
#        filename = "module_quantative_score.pdf",
#        width = 8,
#        height = 8)

for (i in 1:27) {
  cat(i, " ")
  plot <-
    plot_metabolic_module_fpa(
      fpa_result = fpa_result,
      feature_table_marker = feature_table_marker,
      include_feature = TRUE,
      include_hidden_metabolites = FALSE,
      add_compound_name = TRUE,
      metabolic_module_index = i,
      layout = "fr",
      add_pathways = TRUE
    )
  plot
  ggsave(
    plot,
    filename = paste0("module_", i, ".pdf"),
    width = 14,
    height = 8
  )
}

idx <- 4

kegg_id <-
  stringr::str_split(fpa_result$dysregulated_metabolic_module$Total_metabolite_id[idx],
                     "\\{\\}") %>%
  `[[`(1) %>%
  unique()

length(kegg_id)

fpa_result$annotation_table %>%
  dplyr::filter(KEGG.ID %in% kegg_id) %>%
  pull(variable_id) %>%
  unique() %>%
  length()


fpa_result$enriched_pathways_list[[idx]]@result[, c("pathway_name", "p_value_adjust")] %>% head(10)


temp_data <-
  fpa_result$dysregulated_metabolic_module$Total_metabolite_id[fpa_result$dysregulated_metabolic_module$p_value < 0.05] %>%
  purrr::map(function(x) {
    KEGG_ID <-
      stringr::str_split(x, pattern = "\\{\\}") %>%
      unlist()
    
    temp_data <-
      data.frame(KEGG.ID = KEGG_ID) %>%
      dplyr::left_join(annotation_table[, c("variable_id",
                                            "KEGG.ID",
                                            "compound_class",
                                            "Adduct",
                                            "isotope",
                                            "score")], by = "KEGG.ID") %>%
      dplyr::arrange(KEGG.ID) %>%
      dplyr::distinct(KEGG.ID, variable_id, compound_class, .keep_all = TRUE) %>%
      dplyr::filter(!is.na(variable_id)) %>%
      plyr::dlply(.variables = .(KEGG.ID)) %>%
      purrr::map(function(y) {
        if (max(y$score) >= 80) {
          y <-
            y %>%
            dplyr::filter(score >= 80)
        }
        
        remain_compound_class <-
          y %>%
          dplyr::count(compound_class) %>%
          dplyr::filter(n > 1) %>%
          pull(compound_class)
        
        
        if (length(remain_compound_class) > 0) {
          y <-
            y %>% dplyr::filter(compound_class %in% remain_compound_class)
        }
        y %>%
          dplyr::filter(isotope == "[M]")
      })
    
    remain_idx <-
      temp_data %>%
      lapply(nrow) %>%
      unlist() %>%
      `>`(0) %>%
      which()
    
    temp_data <-
      temp_data[remain_idx]
    
    variable_id <-
      temp_data %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::distinct(variable_id, .keep_all = TRUE) %>%
      pull(variable_id)
    
    expression_data[variable_id, ] %>%
      colSums()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(temp_data) <-
  fpa_result$dysregulated_metabolic_module$Name[fpa_result$dysregulated_metabolic_module$p_value < 0.05]

temp_data <-
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "module_name") %>%
  tidyr::pivot_longer(cols = -module_name,
                      names_to = "sample_id",
                      values_to = "intensity") %>%
  dplyr::left_join(urine_metabolomics_data@sample_info[, c("sample_id", "ga_range")], by = "sample_id") %>%
  dplyr::filter(!is.na(ga_range)) %>%
  dplyr::group_by(module_name, ga_range) %>%
  dplyr::summarise(intensity = mean(intensity)) %>%
  dplyr::mutate(ga_range = factor(ga_range, levels = factor(
    ga_range, stringr::str_sort(unique(ga_range), numeric = TRUE)
  ))) %>%
  dplyr::mutate(module_name = factor(module_name, levels = stringr::str_sort(unique(module_name), numeric = TRUE)))

plot <-
  temp_data %>%
  ggplot(aes(ga_range, intensity)) +
  geom_line(aes(group = module_name, color = module_name)) +
  geom_point(aes(color = module_name), size = 2) +
  theme_bw() +
  scale_color_manual(values = colorRampPalette(ggsci::pal_aaas("default")(10))(27)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "GA range", y = "Z-score")

plot

# ggsave(plot,
#        filename = "module_quantitative_plot.pdf",
#        width = 10,
#        height = 9)



####heatmap
temp_data2 <-
  temp_data %>%
  tidyr::pivot_wider(names_from = ga_range, values_from = intensity) %>%
  tibble::column_to_rownames(var = "module_name")

library(ComplexHeatmap)

Heatmap(temp_data2, cluster_columns = FALSE, row_km = 2)


#####module network
temp_data2 <-
  temp_data %>%
  tidyr::pivot_wider(names_from = ga_range, values_from = intensity) %>%
  tibble::column_to_rownames(var = "module_name")


###correlations between modules, and p values
correlation_data <-
  temp_data2 %>%
  t() %>%
  cor(method = "spearman") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "module_name") %>%
  tidyr::pivot_longer(cols = -module_name,
                      names_to = "module_name2",
                      values_to = "correlation") %>%
  dplyr::filter(module_name != module_name2)

p_value <-
  seq_len(nrow(correlation_data)) %>%
  purrr::map(function(i) {
    cor.test(as.numeric(temp_data2[correlation_data$module_name[i], ]), as.numeric(temp_data2[correlation_data$module_name2[i], ]), method = "spearman")$p.value
  }) %>%
  unlist()


correlation_data$p_value <- p_value
correlation_data$p_adjust <- p.adjust(p_value, method = "fdr")

edge_data <-
  correlation_data %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::mutate(abs(correlation) > 0.5)

node_data <-
  data.frame(id = unique(c(
    edge_data$module_name, edge_data$module_name2
  )), class = "module")


library(tidygraph)
library(ggraph)

temp_graph <-
  tbl_graph(nodes = node_data, edges = edge_data) %>%
  dplyr::mutate(degree = centrality_degree(mode = "all"))

plot <-
  temp_graph %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(color = correlation, size = -log(p_adjust, 10)), show.legend = TRUE) +
  geom_node_point(aes(size = degree)) +
  theme_graph() +
  coord_fixed()
plot
ggsave(plot,
       filename = "module_network.pdf",
       width = 10,
       height = 10)
