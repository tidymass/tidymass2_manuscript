library(r4projects)
setwd(get_project_wd())

load("2_data/metabolic_pathways/hmdb_compound_ms1.rda")
load("2_data/metabolic_pathways/hmdb_pathway.rda")
load("2_data/metabolic_pathways/kegg_metabolite_ms1.rda")
load("2_data/metabolic_pathways/kegg_hsa_pathway.rda")
load("2_data/metabolic_pathways/reactome_hsa_pathway.rda")
load("2_data/metabolic_pathways/wikipathway_hsa_pathway.rda")

dir.create("3_data_analysis/5_metabolic_pathway_overlap",
           showWarnings = FALSE)
setwd("3_data_analysis/5_metabolic_pathway_overlap")

####kegg_hsa_pathway
remain_idx <-
  kegg_hsa_pathway@compound_list %>%
  lapply(nrow) %>%
  unlist() %>%
  `>`(0) %>%
  which()

kegg_hsa_pathway <-
  filter_pathway(object = kegg_hsa_pathway, remain_idx = remain_idx)

####hmdb_pathway
remain_idx <-
  hmdb_pathway@compound_list %>%
  lapply(nrow) %>%
  unlist() %>%
  `>`(0) %>%
  which()

hmdb_pathway <-
  filter_pathway(object = hmdb_pathway, remain_idx = remain_idx)


remain_idx <-
  hmdb_pathway@pathway_class %>%
  lapply(function(x) {
    x == "Metabolic;primary_pathway"
  }) %>%
  unlist() %>%
  which()

hmdb_pathway <-
  filter_pathway(object = hmdb_pathway, remain_idx = remain_idx)


length(wikipathway_hsa_pathway)
length(reactome_hsa_pathway)


####wikipathways
remain_idx <-
  wikipathway_hsa_pathway@compound_list %>%
  lapply(nrow) %>%
  unlist() %>%
  `>`(0) %>%
  which()

wikipathway_hsa_pathway <-
  filter_pathway(object = wikipathway_hsa_pathway, remain_idx = remain_idx)

###reactome
remain_idx <-
  reactome_hsa_pathway@compound_list %>%
  lapply(nrow) %>%
  unlist() %>%
  `>`(0) %>%
  which()

reactome_hsa_pathway <-
  filter_pathway(object = reactome_hsa_pathway, remain_idx = remain_idx)

# ####jaccard index inside kegg_hsa_pathway
# kegg_jaccard <-
#   1:(length(kegg_hsa_pathway@compound_list) - 1) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = kegg_hsa_pathway@compound_list[[i]]
#     (i + 1):length(kegg_hsa_pathway@compound_list) %>%
#       purrr::map(function(j) {
#         y = kegg_hsa_pathway@compound_list[[j]]
#         jaccard_index <-
#           length(intersect(x$KEGG.ID, y$KEGG.ID)) / length(unique(c(x$KEGG.ID, y$KEGG.ID)))
#         data.frame(
#           from = kegg_hsa_pathway@pathway_id[i],
#           to = kegg_hsa_pathway@pathway_id[j],
#           from_name = kegg_hsa_pathway@pathway_name[i],
#           to_name = kegg_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(kegg_jaccard, file = "kegg_jaccard.rda")
load("kegg_jaccard.rda")

# ###jaccard index inside hmdb_pathway
# hmdb_jaccard <-
#   1:(length(hmdb_pathway@compound_list) - 1) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = hmdb_pathway@compound_list[[i]]
#     (i + 1):length(hmdb_pathway@compound_list) %>%
#       purrr::map(function(j) {
#         y = hmdb_pathway@compound_list[[j]]
#         jaccard_index <-
#           length(intersect(x$HMDB.ID, y$HMDB.ID)) / length(unique(c(x$HMDB.ID, y$HMDB.ID)))
#         data.frame(
#           from = hmdb_pathway@pathway_id[i],
#           to = hmdb_pathway@pathway_id[j],
#           from_name = hmdb_pathway@pathway_name[i],
#           to_name = hmdb_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(hmdb_jaccard, file = "hmdb_jaccard.rda")
load("hmdb_jaccard.rda")


# ###jaccard inside wikipathway_hsa_pathway
# wikipathway_jaccard <-
#   1:(length(wikipathway_hsa_pathway) - 1) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = wikipathway_hsa_pathway@compound_list[[i]]
#     (i + 1):length(wikipathway_hsa_pathway) %>%
#       purrr::map(function(j) {
#         y = wikipathway_hsa_pathway@compound_list[[j]]
#         jaccard_index <-
#           length(intersect(x$HMDB.ID, y$HMDB.ID)) / length(unique(c(x$HMDB.ID, y$HMDB.ID)))
#         data.frame(
#           from = wikipathway_hsa_pathway@pathway_id[i],
#           to = wikipathway_hsa_pathway@pathway_id[j],
#           from_name = wikipathway_hsa_pathway@pathway_name[i],
#           to_name = wikipathway_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(wikipathway_jaccard, file = "wikipathway_jaccard.rda")
load("wikipathway_jaccard.rda")


# ####reactome
# reactome_jaccard <-
#   1:(length(reactome_hsa_pathway) - 1) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = reactome_hsa_pathway@compound_list[[i]]
#     (i + 1):length(reactome_hsa_pathway) %>%
#       purrr::map(function(j) {
#         y = reactome_hsa_pathway@compound_list[[j]]
#         jaccard_index <-
#           length(intersect(x$CHEBI.ID, y$CHEBI.ID)) / length(unique(c(x$CHEBI.ID, y$CHEBI.ID)))
#         data.frame(
#           from = reactome_hsa_pathway@pathway_id[i],
#           to = reactome_hsa_pathway@pathway_id[j],
#           from_name = reactome_hsa_pathway@pathway_name[i],
#           to_name = reactome_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(reactome_jaccard, file = "reactome_jaccard.rda")
load("reactome_jaccard.rda")


#####kegg vs HMDB
# kegg_hmdb_jaccard <-
#   1:length(kegg_hsa_pathway) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = kegg_hsa_pathway@compound_list[[i]]
#     1:length(hmdb_pathway) %>%
#       purrr::map(function(j) {
#         y = hmdb_pathway@compound_list[[j]]
#         jaccard_index <-
#           length(intersect(x$KEGG.ID, y$KEGG.ID)) / length(unique(c(x$KEGG.ID, y$KEGG.ID)))
#         data.frame(
#           from = kegg_hsa_pathway@pathway_id[i],
#           to = hmdb_pathway@pathway_id[j],
#           from_name = kegg_hsa_pathway@pathway_name[i],
#           to_name = hmdb_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# save(kegg_hmdb_jaccard, file = "kegg_hmdb_jaccard.rda")
load("kegg_hmdb_jaccard.rda")

# ###kegg vs wikipathway
# kegg_wikipathway_jaccard <-
#   1:length(kegg_hsa_pathway) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = kegg_hsa_pathway@compound_list[[i]]
#     1:length(wikipathway_hsa_pathway) %>%
#       purrr::map(function(j) {
#         y = wikipathway_hsa_pathway@compound_list[[j]]
#         intersect_id <-
#           intersect(x$KEGG.ID, y$KEGG.ID)
#         intersect_id <-
#           intersect_id[!is.na(intersect_id)]
#         total_id <-
#           unique(c(x$KEGG.ID, y$KEGG.ID))
#         total_id <-
#           total_id[!is.na(total_id)]
#         jaccard_index <-
#           length(intersect_id) / length(total_id)
#         data.frame(
#           from = kegg_hsa_pathway@pathway_id[i],
#           to = wikipathway_hsa_pathway@pathway_id[j],
#           from_name = kegg_hsa_pathway@pathway_name[i],
#           to_name = wikipathway_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(kegg_wikipathway_jaccard, file = "kegg_wikipathway_jaccard.rda")
load("kegg_wikipathway_jaccard.rda")


# ####kegg vs reactome
# kegg_reactome_jaccard <-
#   1:length(kegg_hsa_pathway) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = kegg_hsa_pathway@compound_list[[i]]
#     1:length(reactome_hsa_pathway) %>%
#       purrr::map(function(j) {
#         y = reactome_hsa_pathway@compound_list[[j]]
#         intersect_id <-
#           intersect(x$KEGG.ID, y$KEGG.ID)
#         intersect_id <-
#           intersect_id[!is.na(intersect_id)]
#         total_id <-
#           unique(c(x$KEGG.ID, y$KEGG.ID))
#         total_id <-
#           total_id[!is.na(total_id)]
#         jaccard_index <-
#           length(intersect_id) / length(total_id)
#         data.frame(
#           from = kegg_hsa_pathway@pathway_id[i],
#           to = reactome_hsa_pathway@pathway_id[j],
#           from_name = kegg_hsa_pathway@pathway_name[i],
#           to_name = reactome_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(kegg_reactome_jaccard, file = "kegg_reactome_jaccard.rda")
load("kegg_reactome_jaccard.rda")


# #####hmdb vs wikipathway
# hmdb_wikipathway_jaccard <-
#   1:length(hmdb_pathway) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = hmdb_pathway@compound_list[[i]]
#     1:length(wikipathway_hsa_pathway) %>%
#       purrr::map(function(j) {
#         y = wikipathway_hsa_pathway@compound_list[[j]]
#         intersect_id <-
#           intersect(x$HMDB.ID, y$HMDB.ID)
#         intersect_id <-
#           intersect_id[!is.na(intersect_id)]
#         total_id <-
#           unique(c(x$HMDB.ID, y$HMDB.ID))
#         total_id <-
#           total_id[!is.na(total_id)]
#         jaccard_index <-
#           length(intersect_id) / length(total_id)
#         data.frame(
#           from = hmdb_pathway@pathway_id[i],
#           to = wikipathway_hsa_pathway@pathway_id[j],
#           from_name = hmdb_pathway@pathway_name[i],
#           to_name = wikipathway_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(hmdb_wikipathway_jaccard, file = "hmdb_wikipathway_jaccard.rda")
load("hmdb_wikipathway_jaccard.rda")


# # ####hmdb vs reactome
# hmdb_reactome_jaccard <-
#   1:length(hmdb_pathway) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = hmdb_pathway@compound_list[[i]]
#     1:length(reactome_hsa_pathway) %>%
#       purrr::map(function(j) {
#         y = reactome_hsa_pathway@compound_list[[j]]
#         intersect_id <-
#           intersect(x$HMDB.ID, y$CHEBI.ID)
#         intersect_id <-
#           intersect_id[!is.na(intersect_id)]
#         total_id <-
#           unique(c(x$HMDB.ID, y$CHEBI.ID))
#         total_id <-
#           total_id[!is.na(total_id)]
#         jaccard_index <-
#           length(intersect_id) / length(total_id)
#         data.frame(
#           from = hmdb_pathway@pathway_id[i],
#           to = reactome_hsa_pathway@pathway_id[j],
#           from_name = hmdb_pathway@pathway_name[i],
#           to_name = reactome_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# save(hmdb_reactome_jaccard, file = "hmdb_reactome_jaccard.rda")
load("hmdb_reactome_jaccard.rda")

# ###wiki vs reactome
# wikipathway_reactome_jaccard <-
#   1:length(wikipathway_hsa_pathway) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = wikipathway_hsa_pathway@compound_list[[i]]
#     1:length(reactome_hsa_pathway) %>%
#       purrr::map(function(j) {
#         y = reactome_hsa_pathway@compound_list[[j]]
#         intersect_id <-
#           intersect(x$ChEBI_ID, y$CHEBI.ID)
#         intersect_id <-
#           intersect_id[!is.na(intersect_id)]
#         total_id <-
#           unique(c(x$ChEBI_ID, y$CHEBI.ID))
#         total_id <-
#           total_id[!is.na(total_id)]
#         jaccard_index <-
#           length(intersect_id) / length(total_id)
#         data.frame(
#           from = wikipathway_hsa_pathway@pathway_id[i],
#           to = reactome_hsa_pathway@pathway_id[j],
#           from_name = wikipathway_hsa_pathway@pathway_name[i],
#           to_name = reactome_hsa_pathway@pathway_name[j],
#           jaccard_index = jaccard_index
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
#
# save(wikipathway_reactome_jaccard, file = "wikipathway_reactome_jaccard.rda")
load("wikipathway_reactome_jaccard.rda")

edge_data <-
  rbind(
    # kegg_jaccard %>% dplyr::mutate(from_dataabse = "KEGG", to_database = "KEGG"),
    # hmdb_jaccard %>% dplyr::mutate(from_dataabse = "HMDB", to_database = "HMDB"),
    # wikipathway_jaccard %>% dplyr::mutate(from_dataabse = "Wikipathway", to_database = "Wikipathway"),
    # reactome_jaccard %>% dplyr::mutate(from_dataabse = "Reactome", to_database = "Reactome"),
    kegg_hmdb_jaccard %>% dplyr::mutate(from_dataabse = "KEGG", to_database = "HMDB"),
    kegg_wikipathway_jaccard %>% dplyr::mutate(from_dataabse = "KEGG", to_database = "Wikipathway"),
    kegg_reactome_jaccard %>% dplyr::mutate(from_dataabse = "KEGG", to_database = "Reactome"),
    hmdb_wikipathway_jaccard %>% dplyr::mutate(from_dataabse = "HMDB", to_database = "Wikipathway"),
    hmdb_reactome_jaccard %>% dplyr::mutate(from_dataabse = "HMDB", to_database = "Reactome"),
    wikipathway_reactome_jaccard %>% dplyr::mutate(from_dataabse = "Wikipathway", to_database = "Reactome")
  ) %>%
  dplyr::filter(jaccard_index > 0.5)

node_data <-
  rbind(
    data.frame(
      node = kegg_hsa_pathway@pathway_id,
      node_name = kegg_hsa_pathway@pathway_name,
      size = lapply(kegg_hsa_pathway@compound_list, nrow) %>%
        unlist(),
      type = "KEGG"
    ),
    data.frame(
      node = hmdb_pathway@pathway_id,
      node_name = hmdb_pathway@pathway_name,
      size = lapply(hmdb_pathway@compound_list, nrow) %>%
        unlist(),
      type = "HMDB"
    ),
    data.frame(
      node = wikipathway_hsa_pathway@pathway_id,
      node_name = wikipathway_hsa_pathway@pathway_name,
      size = lapply(wikipathway_hsa_pathway@compound_list, nrow) %>%
        unlist(),
      type = "Wikipathway"
    ),
    data.frame(
      node = reactome_hsa_pathway@pathway_id,
      node_name = reactome_hsa_pathway@pathway_name,
      size = lapply(reactome_hsa_pathway@compound_list, nrow) %>%
        unlist(),
      type = "Reactome"
    )
    
  )

library(ggraph)
library(tidygraph)

temp_graph <-
  tbl_graph(nodes = node_data, edges = edge_data)

plot_all <-
  ggraph(temp_graph, layout = "kk") +
  geom_edge_link(aes(edge_alpha = jaccard_index), edge_width = 1) +
  geom_node_point(aes(size = size, color = type)) +
  # geom_node_text(aes(label = node_name), repel = TRUE) +
  theme_void() +
  theme()

plot_all
ggsave(plot_all,
       filename = "metabolic_pathway_network_all.pdf",
       width = 12,
       height = 10)


node_data <-
  node_data %>%
  dplyr::filter(node %in% c(unique(edge_data$from), unique(edge_data$to))) %>%
  dplyr::filter(size > 5)

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% node_data$node & to %in% node_data$node)

temp_graph <-
  tbl_graph(nodes = node_data, edges = edge_data)

plot2 <-
  ggraph(temp_graph, layout = "kk") +
  geom_edge_link(aes(edge_alpha = jaccard_index), edge_width = 1) +
  geom_node_point(aes(size = size, color = type)) +
  geom_node_text(aes(label = node_name), repel = TRUE) +
  theme_void() +
  theme()

plot2
ggsave(plot2,
       filename = "metabolic_pathway_network_2.pdf",
       width = 12,
       height = 10)
