# Note: Used release 3.1 of enrichplot. Some of the argument names have been changed.

setwd ("~/Desktop/CS_Paper/") 

library(pheatmap)
library(ggplot2)
library(tidyverse)
library (grid)
library (clusterProfiler)
library (enrichplot)
library(org.Hs.eg.db)
library (tidydr)
# library (ggtreeExtra) # talvez não precise
# library(GOSemSim) # talvez não precise
# library (multienrichjam) # talvez não precise


df <- read.csv("df_brain_4.csv", row.names = 1)

mat <- as.matrix(df)


# Generate the heatmap

p <- pheatmap(df, 
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              na_col = "#7F7F7F",
              color = colorRampPalette(c("#3B9AB2", "white", "#F21A00"))(20),
              border_color = "white",
              fontsize_number = 20)
p

#width=750&height=875
#width=767&height=875

ggsave("fig5c_heatmap3.png", plot = p, width = 7.8, height = 8.75, dpi = 600, limitsize = FALSE)






####################
# Figura vias emap #
###################

setwd ("~/Desktop/CS_Paper/GSEA/mouse") 

full_GOBP <- readRDS("GSEA_GOBP_mouse.rds")

df_GOBP <- as.data.frame(full_GOBP)

unique_terms <- as.data.frame(unique(df_GOBP$Description))

result_df_GOBP <- df_GOBP %>%
  group_by(ID) %>%
  filter(n() >= 3)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_GOBP[!duplicated(result_df_GOBP[[column_name]]), ]

write.csv(df_no_duplicates,"GOBP_common terms.csv")


setwd ("~/Desktop/CS_Paper/") 

selectGOBP <- data.frame("ID" = full_GOBP@compareClusterResult$ID, "name" = full_GOBP@compareClusterResult$Description)
select_GOBP_terms <- df_no_duplicates$ID
selectGOBP <- selectGOBP$ID[!selectGOBP$ID %in% select_GOBP_terms]
granules_GOBP <- dropGO(full_GOBP, term = selectGOBP)
g <- as.data.frame(granules_GOBP)

#xxGOBP <- pairwise_termsim(full_GOBP)  

xxGOBP2 <- pairwise_termsim(granules_GOBP)  

x <- as.data.frame(xxGOBP2)

#write.csv (x, "teste_emap.csv")

selected_pathways_BP <- c(
  "dendrite morphogenesis",
  "lymphocyte chemotaxis",
  "monocyte chemotaxis",
  "cellular response to interleukin-1",
  "neutrophil migration",
  "chemokine-mediated signaling pathway",
  "cell killing",
  "neutrophil chemotaxis",
  "response to chemokine",
  "cellular response to chemokine",
  "granulocyte chemotaxis",
  "positive regulation of cytokine production",
  "type II interferon production",
  "granulocyte migration",
  "vesicle-mediated transport in synapse"
)

emapplot(xxGOBP2, showCategory=selected_pathways_BP, 
         pie.params = list(legend_n = 3),
         cex_label_category = 1.3,layout.params = list(layout = "fr"))

emapplot(xxGOBP2, showCategory=selected_pathways_BP, hclust_method = "average", 
         group_category = TRUE,
         group_legend = TRUE,
         pie.params = list(legend_n = 3),
         cex_label_category = 1.3,layout.params = list(layout = "fr"))

fig5a1 <- emapplot(xxGOBP2, showCategory=selected_pathways_BP, 
                   hclust_method = "average", cex_label_category = 1.3, 
                   pie.params = list(legend_n = 3),
                   layout.params = list(layout = "fr"))

fig5a1 <- emapplot(xxGOBP2, showCategory=selected_pathways_BP, 
                   hclust_method = "average", cex_label_category = 1.3, 
                   pie.params = list(legend_n = 3),
                   layout.params = list(layout = "fr"))

fig5a1 <- emapplot(xxGOBP2, showCategory=selected_pathways_BP
                  )

# determina a cor das bolas

#custom_colors <- c("#F21A00", "#F21A00", "#F21A00" )
#custom_colors <- c("#E03F44", "#E03F44", "#E03F44" )
custom_colors <- c("#F88C7F", "#F88C7F", "#F88C7F" )



fig5a1 +
  scale_fill_manual(values = custom_colors)

fig5A1 <- fig5a1 +
  scale_fill_manual(values = custom_colors)

#width=812&height=905

ggsave("fig5A1_update.png", plot = fig5A1, width = 8.12, height = 9.05, dpi = 600, limitsize = FALSE)

