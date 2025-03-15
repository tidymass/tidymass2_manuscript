library(r4projects)
setwd(get_project_wd())
# source("R/enrich_pathways.R")
# source("R/7_utils.R")
# source("R/11_feature_based_pathway_enrichment.R")
# load("data/hmdb_compound_database.rda")

load("3_data_analysis/1_smartd_project_data_preparation/feature_table_all.rda")
load("3_data_analysis/1_smartd_project_data_preparation/feature_table_marker.rda")
load("3_data_analysis/1_smartd_project_data_preparation/metabolic_network.rda")
load(
  "3_data_analysis/1_smartd_project_data_preparation/urine_metabolomics_data.rda"
)
load("3_data_analysis/1_smartd_project_data_preparation/kegg_hsa_pathway.rda")
load("3_data_analysis/1_smartd_project_data_preparation/kegg_compound_database.rda")

dir.create("3_data_analysis/3_smartd_project_pathway_enrichment",
           showWarnings = FALSE)
setwd("3_data_analysis/3_smartd_project_pathway_enrichment")

library(tidymass)

hmdb_kegg_id <-
  feature_table_marker %>%
  dplyr::left_join(urine_metabolomics_data@variable_info, by = "variable_id") %>%
  dplyr::select(KEGG.ID, HMDB.ID) %>%
  dplyr::filter(!is.na(KEGG.ID) | !is.na(HMDB.ID))

data("kegg_hsa_pathway", package = "metpath")
kegg_hsa_pathway
data("hmdb_pathway", package = "metpath")

# result_kegg =
#   enrich_pathways(
#     query_id = hmdb_kegg_id$KEGG.ID[!is.na(hmdb_kegg_id$KEGG.ID)],
#     query_type = "compound",
#     id_type = "KEGG",
#     pathway_database = kegg_hsa_pathway,
#     only_primary_pathway = FALSE,
#     p_cutoff = 0.05,
#     p_adjust_method = "fdr",
#     threads = 6
#   )
#
# save(result_kegg, file = "result_kegg.rda")
load("result_kegg.rda")
enrich_bar_plot(result_kegg)

pathway_class =
  metpath::pathway_class(hmdb_pathway)

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

remain_idx

hmdb_pathway =
  hmdb_pathway[remain_idx]

# result_hmdb =
#   enrich_pathways(
#     query_id = hmdb_kegg_id$HMDB.ID[!is.na(hmdb_kegg_id$HMDB.ID)],
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = hmdb_pathway,
#     only_primary_pathway = TRUE,
#     p_cutoff = 0.05,
#     p_adjust_method = "fdr",
#     threads = 6
#   )
# save(result_hmdb, file = "result_hmdb.rda")
load("result_hmdb.rda")

enrich_bar_plot(result_kegg, top = 20)
enrich_bar_plot(result_hmdb, top = 20)

plot1 <-
  enrich_scatter_plot(result_kegg)
plot1

ggsave(plot1,
       filename = "kegg_enrichment.pdf",
       width = 5,
       height = 4)

plot2 <-
  enrich_scatter_plot(result_hmdb)
plot
ggsave(plot2,
       filename = "hmdb_enrichment.pdf",
       width = 5,
       height = 4)

result_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05)

result_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05)

data("wikipathway_hsa_pathway", package = "metpath")
wikipathway_hsa_pathway

# result_wikipathway =
#   enrich_pathways(
#     query_id = hmdb_kegg_id$HMDB.ID[!is.na(hmdb_kegg_id$HMDB.ID)],
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = wikipathway_hsa_pathway,
#     only_primary_pathway = FALSE,
#     p_cutoff = 0.05,
#     p_adjust_method = "fdr",
#     threads = 6
#   )
#
# save(result_wikipathway, file = "result_wikipathway.rda")
load("result_wikipathway.rda")
enrich_bar_plot(result_wikipathway)

data("reactome_hsa_pathway", package = "metpath")
reactome_hsa_pathway

# result_reactome =
#   enrich_pathways(
#     query_id = hmdb_kegg_id$HMDB.ID[!is.na(hmdb_kegg_id$HMDB.ID)],
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = reactome_hsa_pathway,
#     only_primary_pathway = FALSE,
#     p_cutoff = 0.05,
#     p_adjust_method = "fdr",
#     threads = 6
#   )
#
# save(result_reactome, file = "result_reactome.rda")
load("result_reactome.rda")
enrich_bar_plot(result_reactome)

id1 <-
  result_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  pull(all_id) %>%
  stringr::str_split(";") %>%
  `[[`(1) %>%
  data.frame(KEGG.ID = .) %>%
  dplyr::left_join(hmdb_compound_database@spectra.info[, c("KEGG.ID", "HMDB.ID")], by = "KEGG.ID") %>%
  pull(HMDB.ID)

id1 <- id1[!is.na(id1)]

id2 <-
  result_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  pull(all_id) %>%
  stringr::str_split(";") %>%
  `[[`(1)

length(intersect(id1, id2)) / length(union(id1, id2))

intersect(id1, id2)

edge1 <-
  data.frame(from = "Steroid hormone biosynthesis", to = id1, class = "KEGG")

node1 <-
  data.frame(
    id = c("Steroid hormone biosynthesis", id1),
    class = c("Pathway", rep("metabolite", length(id1)))
  )


edge2 <-
  data.frame(from = "Steroidogenesis", to = id2, class = 'HMDB')

node2 <-
  data.frame(id = c("Steroidogenesis", id2),
             class = c("Pathway", rep("metabolite", length(id2))))


edge <- rbind(edge1, edge2)
node <- rbind(node1, node2) %>%
  dplyr::distinct(id, .keep_all = TRUE)

temp_graph <-
  tidygraph::tbl_graph(nodes = node, edges = edge)

library(ggraph)

plot <-
  temp_graph %>%
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_arc(strength = 0.1, edge_colour = "black") +
  ggraph::geom_node_point(aes(colour = class, size = class)) +
  geom_node_text(aes(label = ifelse(class == "Pathway", id, NA)), repel = TRUE) +
  scale_color_manual(values = c(
    "Pathway" = "#eeca40",
    "metabolite" = "#23b9c7"
  )) +
  ggraph::theme_graph()
plot
library(extrafont)
loadfonts()
ggsave(plot, filename = "pathway_network.pdf", width = 7, height = 6)
