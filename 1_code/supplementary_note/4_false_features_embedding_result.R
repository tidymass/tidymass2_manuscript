library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(ggplot2)
library(ggforce)
library(dplyr)

load("3_data_analysis/6_supplementary_note/embedding_score.rda")


# ggplot(embedding_result, aes(x = Feature, y = Similarity, fill = Feature)) +
#   geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
#   geom_jitter(width = 0.1, size = 3, alpha = 0.8) +  
#   theme_bw(base_size = 14) +
#   theme(legend.position = "none") +
#   ylab("Biotext Similarity") +
#   xlab("Features") +
#   ggtitle("Biotext Similarity Distribution between FPA Results and Standards") +
#   scale_y_continuous(limits = c(0, 1), expand = c(0, 0))

p <- 
ggplot(embedding_result, aes(x = Feature, y = Similarity, fill = Feature)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +  
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  ylab("Biotext Similarity") +
  xlab("Features") +
  ggtitle("Biotext Similarity Distribution between FPA Results and Standards") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  facet_zoom(ylim = c(0.74, 0.80)) 


ggsave("4_manuscript/Supplementary_note/fpa_validation/biotext_similarity_distribution.pdf", plot = p, width = 12, height = 6, dpi = 300)
