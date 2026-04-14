library (dplyr)
library (ggplot2)
library (grid)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library (tidydr)
library (ggtreeExtra)
library(GOSemSim)
#library (multienrichjam)
library (clusterProfiler)
library (enrichplot)
library(tidyverse)
library(grid)
library(pheatmap)


setwd("~/Desktop/CS_Paper/rev1_figures")

full_GOBP <- readRDS("updown_GOBP.rds")

a <- full_GOBP@compareClusterResult

unique_cluster <- as.data.frame(unique(a$Cluster))

compareCluster_data <- as.data.frame(full_GOBP)

desired_order <- c("up_CSB_CS1AN_Wang_GOBP", "up_CSB_CS1AN_Okur_GOBP", "up_CSB_iPSC_Liu_GOBP",
                   "up_CSB_MSC_Liu_GOBP", "up_CSB_NeuroB_Wang_GOBP", "up_CSA_NeuroB_Wang_GOBP",
                   "up_CS_Patients_GOBP_Wang_2014","up_CSA_Mouse_1m_GOBP","up_CSA_Mouse_12m_GOBP","up_CSA_Mouse_24m_GOBP",
                   "down_CSB_CS1AN_Egly_GOBP",
                   "down_CSB_CS1AN_Wang_GOBP", "down_CSB_CS1AN_Okur_GOBP", "down_CSB_iPSC_Liu_GOBP",
                   "down_CSB_MSC_Liu_GOBP", "down_CSB_NeuroB_Wang_GOBP", "down_CSA_NeuroB_Wang_GOBP",
                   "down_CS_Patients_GOBP_Wang_2014","down_CSA_Mouse_1m_GOBP","down_CSA_Mouse_12m_GOBP","down_CSA_Mouse_24m_GOBP"
                   )

# Reorder the data frame according to the desired cluster names
compareCluster_data$Cluster <- factor(compareCluster_data$Cluster, levels = desired_order)
compareCluster_data <- compareCluster_data[order(compareCluster_data$Cluster), ]

# Update the compareClusterResult object with the reordered data
full_GOBP@compareClusterResult <- compareCluster_data

df_GOBP <- as.data.frame(full_GOBP)

result_df_GOBP <- df_GOBP %>%
  group_by(ID) %>%
  filter(n() >= 9)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_GOBP[!duplicated(result_df_GOBP[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueGOBP_background.csv")

selectGOBP <- data.frame("ID" = full_GOBP@compareClusterResult$ID, "name" = full_GOBP@compareClusterResult$Description)
select_GOBP_terms <- df_no_duplicates$ID
selectGOBP <- selectGOBP$ID[!selectGOBP$ID %in% select_GOBP_terms]
granules_GOBP <- dropGO(full_GOBP, term = selectGOBP)
g <- as.data.frame(granules_GOBP)

xxGOBP <- pairwise_termsim(full_GOBP)  

xxGOBP2 <- pairwise_termsim(granules_GOBP)  

x <- as.data.frame(xxGOBP2)
write.csv(x,"updown_all_GOBP.csv")


dotplot(xxGOBP2)

options(enrichplot.colours = c("#4F1E4D","#A11B58", "#E03F44","#F7C6A6"))


#zissou - find 2 more colors
#custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
#custom_colors <- c("#501E4D", "#B31657" ,"#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
#custom_colors <- c("#501E4D",  "#E1AF00", "#EBCC2A", "#78B7C5","#F21A00", "#3B9AB2", "#B31657")

#custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")


custom_colors <- c("#F21A00","#3B9AB2","#E1AF00","#78B7C5")

treeplot(xxGOBP2, showCategory = 20, offset_tiplab = 4.4, 
         group_color = custom_colors,
         hclust_method = "average", cluster.params = list(n = 4))


fig_updown_gobp <- treeplot(xxGOBP2, showCategory = 20, offset_tiplab = 4.4, 
                            group_color = custom_colors,
                            hclust_method = "average", cluster.params = list(n = 4))
ggsave("fig_updown_gobp.png", plot = fig_updown_gobp, width = 18.5, height = 7.25, dpi = 600, limitsize = FALSE)
ggsave("fig_updown_gobp2.png", plot = fig_updown_gobp, width = 7.5, height = 7.25, dpi = 600, limitsize = FALSE)
ggsave("fig_updown_gobp3.png", plot = fig_updown_gobp, width = 18.3, height = 4.25, dpi = 600, limitsize = FALSE)


x <- as.data.frame(xxGOBP)
write.csv(x,"updown_all_GOBP2.csv")


selected_pathways <- c("extracellular matrix organization",  
"extracellular structure organization",  
"external encapsulating structure organization",  
"axonogenesis",  
"axon guidance",  
"neuron projection guidance",  
"axon development",  
"cell-cell adhesion via plasma-membrane adhesion molecules",  
"synapse organization",  
"chemotaxis",  
"regulation of vasculature development",  
"ameboidal-type cell migration",  
"kidney development",  
"mesenchyme development",  
"ERK1 and ERK2 cascade",  
"angiogenesis",  
"positive regulation of cell migration")  

custom_colors <- c("#E1AF00","#B31657","#78B7C5", "#3B9AB2","#501E4D", "#F21A00")




##########
# Extras #
#########

treeplot(xxGOBP, showCategory = selected_pathways, offset_tiplab = 7.5, 
                            group_color = custom_colors,
                            hclust_method = "average", cluster.params = list(n = 6))

treeplot(xxGOBP, showCategory = selected_pathways, offset_tiplab = 7.5, 
         group_color = custom_colors,
         hclust_method = "average", cluster.params = list(n = 6),clusterPanel.params = list(colnames_angle = 45))

fig_updown_gobp <- treeplot(xxGOBP, showCategory = selected_pathways, offset_tiplab = 7.5, 
                            group_color = custom_colors,
                            hclust_method = "average", cluster.params = list(n = 6),clusterPanel.params = list(colnames_angle = 45))

ggsave("fig_updown_gobp_all.png", plot = fig_updown_gobp, width = 18.3, height = 5.77, dpi = 600, limitsize = FALSE)
ggsave("fig_updown_gobp_all2.png", plot = fig_updown_gobp, width = 5.3, height = 5.77, dpi = 600, limitsize = FALSE)



df <- read.csv("updown_all_GOBP2.csv", row.names = 1)

result <- df %>%
  dplyr::select(Cluster, Description, p.adjust) %>%
  pivot_wider(names_from = Cluster, values_from = p.adjust) %>%
  as.data.frame()

rownames(result) <- result$Description
result$Description <- NULL

df <- result

write.csv(df,"brain_GSEA_heatmap2.csv")

pheatmap(df, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         na_col = "#7F7F7F",
         color = colorRampPalette(c("#4F1E4D","#A11B58","#E03F44","#F7C6A6"))(20),
         #border_color = "white",
         fontsize_number = 20)



fig4a <- treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 9, 
                  group_color = custom_colors,
                  hclust_method = "average", cluster.params = list(n = 6))

ggsave("UPDATE_fig4a.png", plot = fig4a, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)





library (dplyr)
library (ggplot2)
library (grid)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library (tidydr)
library (ggtreeExtra)
library(GOSemSim)
#library (multienrichjam)
library (clusterProfiler)
library (enrichplot)
library(tidyverse)
library(grid)
library(pheatmap)


setwd("~/Desktop/CS_Paper/rev1_figures")

full_KEGG <- readRDS("updown_KEGG_human.rds")

a <- full_KEGG@compareClusterResult

unique_cluster <- as.data.frame(unique(a$Cluster))

compareCluster_data <- as.data.frame(full_KEGG)

desired_order <- c("up_CSB_CS1AN_Egly_KEGG","up_CSB_CS1AN_Wang_KEGG", "up_CSB_CS1AN_Okur_KEGG", "up_CSB_iPSC_Liu_KEGG",
                   "up_CSB_MSC_Liu_KEGG", "up_CSB_NeuroB_Wang_KEGG", "up_CSA_NeuroB_Wang_KEGG",
                   "up_CS_Patients_KEGG_Wang_2014", 
                   "down_CSB_CS1AN_2014_Egly_Fibroblast_CS1AN",
                   "down_CSB_CS1AN_Wang_KEGG", "down_CSB_CS1AN_Okur_KEGG", "down_CSB_iPSC_Liu_KEGG",
                   "down_CSB_MSC_Liu_KEGG", "down_CSB_NeuroB_Wang_KEGG", "down_CSA_NeuroB_Wang_KEGG",
                   "down_CS_Patients_KEGG_Wang_2014"
)

# Reorder the data frame according to the desired cluster names
compareCluster_data$Cluster <- factor(compareCluster_data$Cluster, levels = desired_order)
compareCluster_data <- compareCluster_data[order(compareCluster_data$Cluster), ]

# Update the compareClusterResult object with the reordered data
full_KEGG@compareClusterResult <- compareCluster_data

df_KEGG <- as.data.frame(full_KEGG)

result_df_KEGG <- df_KEGG %>%
  group_by(ID) %>%
  filter(n() >= 6)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_KEGG[!duplicated(result_df_KEGG[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueKEGG_background.csv")

selectKEGG <- data.frame("ID" = full_KEGG@compareClusterResult$ID, "name" = full_KEGG@compareClusterResult$Description)
select_KEGG_terms <- df_no_duplicates$ID
selectKEGG <- selectKEGG$ID[!selectKEGG$ID %in% select_KEGG_terms]
granules_KEGG <- dropGO(full_KEGG, term = selectKEGG)
g <- as.data.frame(granules_KEGG)

#xxKEGG <- pairwise_termsim(full_KEGG)  

xxKEGG2 <- pairwise_termsim(granules_KEGG)  

x <- as.data.frame(xxKEGG2)
write.csv(x,"updown_all_KEGG.csv")


dotplot(xxKEGG2)

options(enrichplot.colours = c("#4F1E4D","#A11B58", "#E03F44","#F7C6A6"))


#zissou - find 2 more colors
#custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
#custom_colors <- c("#501E4D", "#B31657" ,"#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
#custom_colors <- c("#501E4D",  "#E1AF00", "#EBCC2A", "#78B7C5","#F21A00", "#3B9AB2", "#B31657")

#custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")


custom_colors <- c("#F21A00","#3B9AB2","#E1AF00","#78B7C5")


custom_colors <- c("#3B9AB2","#501E4D","#B31657","#E1AF00","#78B7C5" , "#F21A00")


treeplot(xxKEGG2, offset_tiplab = 5.4, 
         group_color = custom_colors,
         hclust_method = "average", cluster.params = list(n = 6))

fig_updown_KEGG <- treeplot(xxKEGG2, offset_tiplab = 5.4, 
                            group_color = custom_colors,
                            hclust_method = "average", cluster.params = list(n = 6))

ggsave("fig_updown_KEGG_all.png", plot = fig_updown_KEGG, width = 15.47, height = 4, dpi = 600, limitsize = FALSE)
ggsave("fig_updown_KEGG_all3.png", plot = fig_updown_KEGG, width = 18.47, height = 4, dpi = 600, limitsize = FALSE)
ggsave("fig_updown_KEGG_all2.png", plot = fig_updown_KEGG, width = 6.3, height = 4, dpi = 600, limitsize = FALSE)

treeplot(xxKEGG2, offset_tiplab = 7.5, 
         group_color = custom_colors,
         hclust_method = "average", cluster.params = list(n = 6),clusterPanel.params = list(colnames_angle = 45))

fig_updown_KEGG <- treeplot(xxKEGG, showCategory = selected_pathways, offset_tiplab = 7.5, 
                            group_color = custom_colors,
                            hclust_method = "average", cluster.params = list(n = 6),clusterPanel.params = list(colnames_angle = 45))

ggsave("fig_updown_KEGG_all.png", plot = fig_updown_KEGG, width = 18.3, height = 5.77, dpi = 600, limitsize = FALSE)
ggsave("fig_updown_KEGG_all2.png", plot = fig_updown_KEGG, width = 5.3, height = 5.77, dpi = 600, limitsize = FALSE)



df <- read.csv("updown_all_KEGG2.csv", row.names = 1)

result <- df %>%
  dplyr::select(Cluster, Description, p.adjust) %>%
  pivot_wider(names_from = Cluster, values_from = p.adjust) %>%
  as.data.frame()

rownames(result) <- result$Description
result$Description <- NULL

df <- result

write.csv(df,"brain_GSEA_heatmap2.csv")

pheatmap(df, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         na_col = "#7F7F7F",
         color = colorRampPalette(c("#4F1E4D","#A11B58","#E03F44","#F7C6A6"))(20),
         #border_color = "white",
         fontsize_number = 20)



fig4a <- treeplot(xxKEGG2, showCategory = selected_pathways, offset_tiplab = 9, 
                  group_color = custom_colors,
                  hclust_method = "average", cluster.params = list(n = 6))

ggsave("UPDATE_fig4a.png", plot = fig4a, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)

