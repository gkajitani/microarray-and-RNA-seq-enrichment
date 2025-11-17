setwd("~/Desktop/CS_Paper/rev1_figures")

library(ggplot2)
library(tidyverse)
library(grid)
library(pheatmap)


df <- read.csv("brain_GSEA_heatmap.csv")

result <- df %>%
  dplyr::select(Cluster, Description, NES) %>%
  pivot_wider(names_from = Cluster, values_from = NES) %>%
  as.data.frame()

rownames(result) <- result$Description
result$Description <- NULL

df <- result

write.csv(df,"brain_GSEA_heatmap2.csv")

df <- read.csv("brain_GSEA_heatmap2.csv", row.names = 1)


pheatmap(df, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         na_col = "#7F7F7F",
         color = colorRampPalette(c("#3B9AB2", "white", "#F21A00"))(20),
         #border_color = "white",
         fontsize_number = 20)

p <- pheatmap(df, 
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              na_col = "#7F7F7F",
              color = colorRampPalette(c("#3B9AB2", "white", "#F21A00"))(20),
             # border_color = "white",
              fontsize_number = 20)
p

ggsave("GSEA_brain.png", plot = p, width = 7.13, height = 7.53, dpi = 600, limitsize = FALSE)



#########
# GO:CC #
#########


df <- read.csv("brain_GSEA_GOCC_heatmap.csv")

result <- df %>%
  dplyr::select(Cluster, Description, NES) %>%
  pivot_wider(names_from = Cluster, values_from = NES) %>%
  as.data.frame()

rownames(result) <- result$Description
result$Description <- NULL

df <- result

write.csv(df,"brain_GSEA_GOCC_heatmap2.csv")

df <- read.csv("brain_GSEA_GOCC_heatmap2.csv", row.names = 1)


pheatmap(df, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         na_col = "#7F7F7F",
         color = colorRampPalette(c("#3B9AB2", "white", "#F21A00"))(20),
         #border_color = "white",
         fontsize_number = 20)

p <- pheatmap(df, 
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              na_col = "#7F7F7F",
              color = colorRampPalette(c("#3B9AB2", "white", "#F21A00"))(20),
              # border_color = "white",
              fontsize_number = 20)
p

ggsave("GSEA_brain_GOCC.png", plot = p, width = 8.23, height = 4.95, dpi = 600, limitsize = FALSE)
