library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(ggplot2)
library(dplyr)

load("3_data_analysis/6_supplementary_note/embedding_score.rda")


ggplot(embedding_result, aes(x = Feature, y = Similarity, fill = Feature)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +  
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  ylab("Biotext Similarity") +
  xlab("Features") +
  ggtitle("Biotext Similarity Distribution between FPA Results and Standards")
