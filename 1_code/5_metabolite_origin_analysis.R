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

length(human_metabolite)

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
# ggsave(plot,
#        filename = "environment_specific_metabolite_network.pdf",
#        width = 10,
#        height = 10)


plot <-
  source_network(
    object = object,
    source_id = c("Human"),
    top_specific_source = 3
  )

plot
# ggsave(plot,
#        filename = "human_specific_metabolite_network.pdf",
#        width = 10,
#        height = 10)


###Box plot for all these metabolites
###trimesters 1 - 12 weeks,
###trimesters 2 - 13-27 weeks,
###trimesters 3 - 28-40 weeks
sample_info <-
  urine_metabolomics_data@sample_info %>%
  select(sample_id, GA) %>%
  dplyr::filter(!is.na(GA)) %>%
  dplyr::mutate(trimester =
                  case_when(
                    GA < 13 ~ "Trimester 1",
                    GA >= 13 &
                      GA < 28 ~ "Trimester 2",
                    GA >= 28 ~ "Trimester 3"
                  ))

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
  )]
# dplyr::left_join(feature_table_marker[, c("variable_id", "degree", "p_value")], by = "variable_id")

####calculate the different between trimester 2 and trimester 3
fc_p <-
  1:length(specific_metabolite_bacteria$variable_id) %>%
  purrr::map(function(i) {
    cat(i, " ")
    temp <-
      urine_metabolomics_data@expression_data[specific_metabolite_bacteria$variable_id[i], sample_info$sample_id] %>%
      t() %>%
      as.data.frame() %>%
      cbind(sample_info) %>%
      as.data.frame()
    
    colnames(temp)[1] <- "value"
    
    p_value <-
      wilcox.test(temp$value[temp$trimester == "Trimester 2"], temp$value[temp$trimester == "Trimester 3"])$p.value
    
    fc <-
      mean(temp$value[temp$trimester == "Trimester 3"], na.rm = TRUE) /
      mean(temp$value[temp$trimester == "Trimester 2"], na.rm = TRUE)
    
    data.frame(
      variable_id = specific_metabolite_bacteria$variable_id[i],
      p_value = p_value,
      fc = fc
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

specific_metabolite_bacteria <-
  specific_metabolite_bacteria %>%
  dplyr::left_join(fc_p, by = "variable_id") %>%
  dplyr::arrange(p_value)
library(openxlsx)
# openxlsx::write.xlsx(
#   specific_metabolite_bacteria,
#   file = "specific_metabolite_bacteria.xlsx",
#   rowNames = FALSE
# )

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
  )]

####calculate the different between trimester 2 and trimester 3
fc_p <-
  1:length(specific_metabolite_environment$variable_id) %>%
  purrr::map(function(i) {
    cat(i, " ")
    temp <-
      urine_metabolomics_data@expression_data[specific_metabolite_environment$variable_id[i], sample_info$sample_id] %>%
      t() %>%
      as.data.frame() %>%
      cbind(sample_info) %>%
      as.data.frame()
    
    colnames(temp)[1] <- "value"
    
    p_value <-
      wilcox.test(temp$value[temp$trimester == "Trimester 2"], temp$value[temp$trimester == "Trimester 3"])$p.value
    
    fc <-
      mean(temp$value[temp$trimester == "Trimester 3"], na.rm = TRUE) /
      mean(temp$value[temp$trimester == "Trimester 2"], na.rm = TRUE)
    
    data.frame(
      variable_id = specific_metabolite_environment$variable_id[i],
      p_value = p_value,
      fc = fc
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


specific_metabolite_environment <-
  specific_metabolite_environment %>%
  dplyr::left_join(fc_p, by = "variable_id") %>%
  dplyr::arrange(p_value)


# openxlsx::write.xlsx(
#   specific_metabolite_environment,
#   file = "specific_metabolite_environment.xlsx",
#   rowNames = FALSE
# )


####Box plot for bacteria metabolites
dir.create("bacteria_metabolite_boxplot", showWarnings = FALSE)

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



for (i in 1:length(specific_metabolite_bacteria$variable_id)) {
  cat(i, " ")
  if (specific_metabolite_bacteria$p_value[i] > 0.05) {
    next()
  }
  temp <-
    temp_data[specific_metabolite_bacteria$variable_id[i], sample_info$sample_id] %>%
    t() %>%
    as.data.frame() %>%
    cbind(sample_info) %>%
    as.data.frame()
  colnames(temp)[1] <- "value"
  
  library(gghalves)
  library(ggsignif)
  
  plot <-
    temp %>%
    ggplot(aes(trimester, value)) +
    geom_half_boxplot(
      side = "l",
      alpha = 0.5,
      aes(fill = trimester),
      show.legend = FALSE
    ) +
    geom_half_dotplot(
      aes(fill = trimester),
      binwidth = 0.1,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_half_violin(
      aes(fill = trimester),
      side = "r",
      alpha = 0.5,
      show.legend = FALSE
    ) +
    theme_bw() +
    ggsignif::geom_signif(
      comparisons = list(
        c("Trimester 1", "Trimester 2"),
        c("Trimester 2", "Trimester 3")
      ),
      test = "wilcox.test",
      map_signif_level = TRUE,
      textsize = 5,
      size = 0.5
    ) +
    labs(x = "",
         y =  "Scale intensity",
         title = specific_metabolite_bacteria$Compound.name[i])
  
  ggsave(
    plot,
    filename = paste0(
      "bacteria_metabolite_boxplot/",
      specific_metabolite_bacteria$variable_id[i],
      ".pdf"
    ),
    width = 5,
    height = 5
  )
}

####Box plot for environment metabolites
dir.create("environment_metabolite_boxplot", showWarnings = FALSE)
temp_data <-
  urine_metabolomics_data@expression_data
###environment
temp_data <-
  temp_data[, sample_info$sample_id] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()
for (i in 1:length(specific_metabolite_environment$variable_id)) {
  cat(i, " ")
  if (specific_metabolite_environment$p_value[i] > 0.05) {
    next()
  }
  temp <-
    temp_data[specific_metabolite_environment$variable_id[i], sample_info$sample_id] %>%
    t() %>%
    as.data.frame() %>%
    cbind(sample_info) %>%
    as.data.frame()
  colnames(temp)[1] <- "value"
  
  library(gghalves)
  library(ggsignif)
  
  plot <-
    temp %>%
    ggplot(aes(trimester, value)) +
    geom_half_boxplot(
      side = "l",
      alpha = 0.5,
      aes(fill = trimester),
      show.legend = FALSE
    ) +
    geom_half_dotplot(
      aes(fill = trimester),
      binwidth = 0.1,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_half_violin(
      aes(fill = trimester),
      side = "r",
      alpha = 0.5,
      show.legend = FALSE
    ) +
    theme_bw() +
    ggsignif::geom_signif(
      comparisons = list(
        c("Trimester 1", "Trimester 2"),
        c("Trimester 2", "Trimester 3")
      ),
      test = "wilcox.test",
      map_signif_level = TRUE,
      textsize = 5,
      size = 0.5
    ) +
    labs(x = "",
         y =  "Scale intensity",
         title = specific_metabolite_environment$Compound.name[i])
  
  ggsave(
    plot,
    filename = paste0(
      "environment_metabolite_boxplot/",
      specific_metabolite_environment$variable_id[i],
      ".pdf"
    ),
    width = 5,
    height = 5
  )
}

plot <-
  data.frame(part = sort(unlist(
    stringr::str_split(object@annotation_table$from_which_part, "\\{\\}")
  ))) %>%
  dplyr::count(part) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(part != "Unknown") %>%
  dplyr::mutate(part = factor(part, levels = part)) %>%
  dplyr::filter(n > 4) %>%
  ggplot(aes(part, n)) +
  geom_bar(stat = "identity",
           fill = metabolite_source_color["Human"],
           color = "black") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Metabolite number")

plot

ggsave(plot,
       filename = "metabolite_human_origin_part.pdf",
       width = 6,
       height = 5)

plot <-
  specific_source_network(
    object = object,
    specific_source_id = c("Urine", "Blood"),
    top_specific_source = 3
  )
plot
ggsave(plot,
       filename = "urine_blood_specific_metabolite_network.pdf",
       width = 7,
       height = 7)

temp <-
  annotation_table %>%
  dplyr::filter(!is.na(from_which_part))

blood_idx <-
  stringr::str_detect(as.character(temp$from_which_part), "Blood") %>%
  which()
urine_idx <-
  stringr::str_detect(temp$from_which_part, "Urine") %>%
  which()

###venn diagram
venn <-
  VennDiagram::venn.diagram(
    x = list(temp$Lab.ID[blood_idx], temp$Lab.ID[urine_idx]),
    category.names = c("Blood", "Urine"),
    filename = NULL,
    output = FALSE,
    fill = c("#f50538", "#e5c185")
  )

grid::grid.draw(venn)

pdf("metabolite_origin_venn_blood_urine.pdf",
    width = 5,
    height = 5)
grid::grid.draw(venn)
dev.off()

specific_source_network(
  object = object,
  specific_source_id = c("Urine"),
  top_specific_source = 3
)



metabolite_origin_network(
  object = object,
  metabolite_id = object@annotation_table$Lab.ID,
  top_specific_source = 3
)
