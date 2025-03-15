library(r4projects)
setwd(get_project_wd())

load("3_data_analysis/1_smartd_project_data_preparation/feature_table_all.rda")
load("3_data_analysis/1_smartd_project_data_preparation/feature_table_marker.rda")
load("3_data_analysis/1_smartd_project_data_preparation/metabolic_network.rda")
load(
  "3_data_analysis/1_smartd_project_data_preparation/urine_metabolomics_data.rda"
)
load("3_data_analysis/1_smartd_project_data_preparation/kegg_hsa_pathway.rda")
load("3_data_analysis/1_smartd_project_data_preparation/kegg_compound_database.rda")


dir.create("3_data_analysis/2_smartd_project_heatmap", showWarnings = FALSE)
setwd("3_data_analysis/2_smartd_project_heatmap")

library(ComplexHeatmap)

temp_data <-
  urine_metabolomics_data[c(feature_table_marker$variable_id), ] %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(ga_range)) %>%
  scale_data(center = TRUE) %>%
  pivot_longer() %>%
  dplyr::left_join(urine_metabolomics_data@sample_info[, c("sample_id", "subject_id", "GA")], by = c("sample_id")) %>%
  dplyr::select(variable_id, GA, value, sample_id) %>%
  dplyr::mutate(GA = case_when(is.na(GA) ~ 42, TRUE ~ GA)) %>%
  dplyr::arrange(GA) %>%
  dplyr::select(-GA) %>%
  pivot_wider(names_from = sample_id, values_from = value) %>%
  tibble::column_to_rownames(var = "variable_id")

library(circlize)

range(temp_data)

###use more beautifule colors

library(viridis)
col_fun <- circlize::colorRamp2(seq(-3, 3, length.out = 11), rev(RColorBrewer::brewer.pal(name = "RdYlBu", n = 11)))

###Split columns by ga

column_split <-
  data.frame(sample_id = colnames(temp_data)) %>%
  dplyr::left_join(urine_metabolomics_data@sample_info[, c("sample_id", "ga_range")], by = c("sample_id")) %>%
  dplyr::mutate(ga_range = case_when(ga_range == "PP" ~ "PP", ga_range != "PP" ~ "During pregnancy")) %>%
  pull(ga_range)

plot <-
  ComplexHeatmap::Heatmap(
    temp_data,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    border = TRUE,
    column_names_rot = 45,
    row_km = 2,
    name = "Z-score",
    col = col_fun,
    show_column_names = FALSE,
    column_split = column_split
  )

plot <- ggplotify::as.ggplot(plot)
plot

# ggsave(plot,
#        filename = "heatmap.pdf",
#        width = 10,
#        height = 8)

####venn diagram
marker_up <-
  feature_table_marker %>%
  dplyr::filter(degree > 0)

marker_down <-
  feature_table_marker %>%
  dplyr::filter(degree < 0)


dim(marker_up)
dim(marker_down)

library(VennDiagram)

####up marker
identified_up_marker <-
  marker_up %>%
  dplyr::left_join(urine_metabolomics_data@variable_info,
                   by = c("variable_id" = "variable_id")) %>%
  dplyr::filter(!is.na(Level)) %>%
  pull(variable_id)

up_marker <-
  marker_up %>%
  pull(variable_id)

venn_up <-
  VennDiagram::venn.diagram(
    x = list(identified_up_marker, up_marker),
    category.names = c("Identified", "All Up metabolic features"),
    filename = NULL,
    output = FALSE,
    fill = c("#d9042b", "#5e606c")
  )

# pdf("up_marker_venn.pdf", width = 5.2, height = 5)
# grid.draw(venn_up)
# dev.off()




####down marker
identified_down_marker <-
  marker_down %>%
  dplyr::left_join(urine_metabolomics_data@variable_info,
                   by = c("variable_id" = "variable_id")) %>%
  dplyr::filter(!is.na(Level)) %>%
  pull(variable_id)

down_marker <-
  marker_down %>%
  pull(variable_id)

venn_down <-
  VennDiagram::venn.diagram(
    x = list(identified_down_marker, down_marker),
    category.names = c("Identified", "All down metabolic features"),
    filename = NULL,
    output = FALSE,
    fill = c("#d9042b", "#5e606c")
  )

pdf("down_marker_venn.pdf", width = 5.2, height = 5)
grid.draw(venn_down)
dev.off()

####recovery_score for each feature
pp_idx <-
  which(urine_metabolomics_data@sample_info$ga_range == "PP")
last_one <-
  which(urine_metabolomics_data@sample_info$ga_range == "(38,42]")

expression_data <-
  urine_metabolomics_data@expression_data

recovery_score <-
  feature_table_marker$variable_id %>%
  purrr::map(function(x) {
    mean(as.numeric(urine_metabolomics_data@expression_data[x, pp_idx])) / mean(as.numeric(urine_metabolomics_data@expression_data[x, last_one]))
  }) %>%
  unlist()

feature_table_marker$recovery_score <- recovery_score

plot <-
  feature_table_marker %>%
  dplyr::mutate(direction = case_when(degree > 0 ~ "up", degree < 0 ~ "down")) %>%
  ggplot(aes(x = degree, y = recovery_score)) +
  geom_point(
    aes(fill = direction, size = -log(p_value)),
    shape = 21,
    color = "black",
    alpha = 0.7
  ) +
  theme_bw() +
  labs(x = "SAM score", y = "Recovery score") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(
    "up" = RColorBrewer::brewer.pal(name = "RdYlBu", n = 11)[1],
    "down" = RColorBrewer::brewer.pal(name = "RdYlBu", n = 11)[11]
  ))

plot

ggsave(plot,
       filename = "sam_vs_recovery_score.pdf",
       width = 7,
       height = 5)

temp_data <-
  feature_table_marker %>%
  dplyr::mutate(direction = case_when(degree > 0 ~ "up", degree < 0 ~ "down"))

sum(temp_data$direction == "up" & temp_data$recovery_score > 1)/sum(temp_data$direction == "up")
sum(temp_data$direction == "up" & temp_data$recovery_score <= 1)/sum(temp_data$direction == "up")


sum(temp_data$direction == "down" & temp_data$recovery_score > 1)/sum(temp_data$direction == "down")
sum(temp_data$direction == "down" & temp_data$recovery_score <= 1)/sum(temp_data$direction == "down")
