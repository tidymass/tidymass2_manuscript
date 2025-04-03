library(r4projects)
setwd(get_project_wd())

load("3_data_analysis/1_smartd_project_data_preparation/feature_table_marker.rda")
load("3_data_analysis/1_smartd_project_data_preparation/feature_table_all.rda")

load(
  "3_data_analysis/1_smartd_project_data_preparation/urine_metabolomics_data.rda"
)

dir.create("3_data_analysis/5_metabolite_origin_analysis",
           showWarnings = FALSE)
setwd("3_data_analysis/5_metabolite_origin_analysis")

library(metid)

object <-
  analyze_metabolite_origins(object = urine_metabolomics_data)

###Overlap
plot <-
  metabolite_origin_upsetplot(object = object)

plot

# ggsave(plot,
#        filename = "metabolite_origin_upsetplot_1.pdf",
#        width = 14,
#        height = 7)

plot <-
  metabolite_origin_upsetplot(object = object, min_size = 5)

plot

# ggsave(plot,
#        filename = "metabolite_origin_upsetplot_2.pdf",
#        width = 7,
#        height = 5)

###metabolite_origin_network
plot <-
  source_metabolite_network(object = object, circle = TRUE)
plot <-
  plot +
  coord_fixed()
plot
extrafont::loadfonts()
# ggsave(plot,
#        filename = "source_metabolite_network.pdf",
#        width = 10,
#        height = 10)

#####human and bacteria
annotation_table <-
  urine_metabolomics_data@annotation_table

human_metabolite <-
  annotation_table %>%
  dplyr::filter(from_human == "Yes") %>%
  dplyr::pull(Lab.ID)

bacteria_metabolite <-
  annotation_table %>%
  dplyr::filter(from_bacteria == "Yes") %>%
  dplyr::pull(Lab.ID)

environment_metabolite <-
  annotation_table %>%
  dplyr::filter(from_environment == "Yes") %>%
  dplyr::pull(Lab.ID)

metabolite_source_color <- c(
  "Human" = "#2c61a1",
  "Plant" = "#78938a",
  "Food" = "#f5eddc",
  "Bacteria" = "#0f1c5c",
  "Animal" = "#d2b48c",
  "Enviornment" = "#8f354b",
  "Drug" = "#000000"
)

##Venn diagram
venn <-
  VennDiagram::venn.diagram(
    x = list(
      human_metabolite,
      bacteria_metabolite,
      environment_metabolite
    ),
    category.names = c("Human", "Bacteria", "Environment"),
    filename = NULL,
    output = FALSE,
    fill = c(
      metabolite_source_color["Human"],
      metabolite_source_color["Bacteria"],
      metabolite_source_color["Enviornment"]
    )
  )

grid::grid.draw(venn)

# pdf("metabolite_origin_venn.pdf",
#     width = 5,
#     height = 5)
# grid::grid.draw(venn)
# dev.off()

####specific metabolites in bacteria and environment
plot <-
  source_network(
    object = object,
    source_id = c("Bacteria"),
    top_specific_source = 3
  )
plot
# ggsave(plot,
#        filename = "bacteria_specific_metabolite_network.pdf",
#        width = 10,
#        height = 10)

plot <-
  source_network(
    object = object,
    source_id = c("Environment"),
    top_specific_source = 3
  )

plot
ggsave(plot,
       filename = "environment_specific_metabolite_network.pdf",
       width = 10,
       height = 10)


plot <-
  source_network(
    object = object,
    source_id = c("Human"),
    top_specific_source = 3
  )

plot
ggsave(plot,
       filename = "human_specific_metabolite_network.pdf",
       width = 10,
       height = 10)



#####Specific metabolites from bacteria
index <-
  which(annotation_table$from_bacteria == "Yes" &
          annotation_table$from_human != "Yes")

specific_metabolite_bacteria <-
  annotation_table[index, c(
    "variable_id",
    "Lab.ID",
    "from_bacteria",
    "from_human",
    "Compound.name",
    "HMDB.ID",
    "KEGG.ID"
  )] %>%
  dplyr::left_join(feature_table_marker[, c("variable_id", "degree", "p_value")], by = "variable_id")
# write.csv(specific_metabolite_bacteria, "specific_metabolite_bacteria.csv", row.names = FALSE)

##Specific metabolites from environment
index <-
  which(annotation_table$from_environment == "Yes" &
          annotation_table$from_human != "Yes")

specific_metabolite_environment <-
  annotation_table[index, c(
    "variable_id",
    "Lab.ID",
    "from_environment",
    "from_human",
    "Compound.name",
    "HMDB.ID",
    "KEGG.ID"
  )] %>%
  dplyr::left_join(feature_table_marker[, c("variable_id", "degree", "p_value")], by = "variable_id")
# write.csv(specific_metabolite_environment, "specific_metabolite_environment.csv", row.names = FALSE)

###Box plot for all these metabolites
###trimesters 1 - 12 weeks,
###trimesters 2 - 13-27 weeks,
###trimesters 3 - 28-40 weeks
sample_info <-
  urine_metabolomics_data@sample_info %>%
  select(sample_id, GA) %>%
  dplyr::filter(!is.na(GA)) %>%
  dplyr::mutate(trimester =
                  case_when(GA < 13 ~ "t1", GA >= 13 &
                              GA < 28 ~ "t2", GA >= 28 ~ "t3"))

temp_data <-
  urine_metabolomics_data@expression_data

###bacteria
temp_data <-
  temp_data[, sample_info$sample_id] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

####Box plot for bacteria metabolites


metabolite_origin_network(
  object = object,
  metabolite_id = object@annotation_table$Lab.ID,
  top_specific_source = 3
)



specific_source_network(
  object = object,
  specific_source_id = c("Urine"),
  top_specific_source = 3
)
