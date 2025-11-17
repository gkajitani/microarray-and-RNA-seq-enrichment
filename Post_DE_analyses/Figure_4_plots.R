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
library(ggtangle)


# Colocar cores mais claras
options(enrichplot.colours = c("#501E4D","#791F59", "#9F1A5B", "#B31657","#DE2E43"))

###############################################
###############################################
###############################################
###############################################
###############################################

# Note: These were the first versions of the plots in the Figure 4. Further adjustments using image editing software were made to combine them.


##############
# Human GOBP #
##############


setwd ("~/Desktop/CS_Paper/") 


full_GOBP <- readRDS("Humans_GOBP_background.rds")

setwd ("~/Desktop/CS_Paper")


a <- full_GOBP@compareClusterResult

unique_cluster <- as.data.frame(unique(a$Cluster))

compareCluster_data <- as.data.frame(full_GOBP)

desired_order <- c("CSB_CS1AN_Wang_GOBP", "CSB_CS1AN_Okur_GOBP", "CSB_iPSC_Liu_GOBP",
                   "CSB_MSC_Liu_GOBP", "CSB_NeuroB_Wang_GOBP", "CSA_NeuroB_Wang_GOBP",
                   "CS_Patients_GOBP_Wang_2014")

# Reorder the data frame according to the desired cluster names
compareCluster_data$Cluster <- factor(compareCluster_data$Cluster, levels = desired_order)
compareCluster_data <- compareCluster_data[order(compareCluster_data$Cluster), ]

# Update the compareClusterResult object with the reordered data
full_GOBP@compareClusterResult <- compareCluster_data


df_GOBP <- as.data.frame(full_GOBP)

result_df_GOBP <- df_GOBP %>%
  group_by(ID) %>%
  filter(n() >= 6)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_GOBP[!duplicated(result_df_GOBP[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueGOBP_background.csv")

selectGOBP <- data.frame("ID" = full_GOBP@compareClusterResult$ID, "name" = full_GOBP@compareClusterResult$Description)
select_GOBP_terms <- df_no_duplicates$ID
selectGOBP <- selectGOBP$ID[!selectGOBP$ID %in% select_GOBP_terms]
granules_GOBP <- dropGO(full_GOBP, term = selectGOBP)
g <- as.data.frame(granules_GOBP)

#xxGOBP <- pairwise_termsim(full_GOBP)  


xxGOBP2 <- pairwise_termsim(granules_GOBP)  

x <- as.data.frame(xxGOBP2)

dotplot(xxGOBP2)

treeplot(xxGOBP2)

options(enrichplot.colours = c("#4F1E4D","#A11B58", "#E03F44","#F7C6A6"))


#zissou - find 2 more colors
#custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
#custom_colors <- c("#501E4D", "#B31657" ,"#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
#custom_colors <- c("#501E4D",  "#E1AF00", "#EBCC2A", "#78B7C5","#F21A00", "#3B9AB2", "#B31657")

#custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")


selected_pathways <- c("kidney development",
                       "renal system development",
                       "mesenchymal cell differentiation",
                       "mesenchyme development",
                       "muscle cell migration",
                       "endothelial cell migration",
                       "smooth muscle cell migration",
                       "nephron development",
                       "regulation of angiogenesis",
                       "regulation of vasculature development",
                       "muscle cell differentiation",
                       "kidney epithelium development",
                       "vascular process in circulatory system",
                       "blood circulation",
                       "mesenchymal cell migration",
                       "muscle contraction"
)

custom_colors <- c("#B31657","#501E4D","#E1AF00","#F21A00" , "#3B9AB2","#78B7C5" )

treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 9, 
         group_color = custom_colors,
         hclust_method = "average", cluster.params = list(n = 6))

fig4a <- treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 9, 
                  group_color = custom_colors,
                  hclust_method = "average", cluster.params = list(n = 6))

ggsave("UPDATE_fig4a.png", plot = fig4a, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)


selected_pathways <- c("positive regulation of MAPK cascade",
                       "ERK1 and ERK2 cascade",
                       "positive regulation of cell development",
                       "response to xenobiotic stimulus",
                       "cell-cell adhesion via plasma-membrane adhesion molecules",
                       "negative regulation of cell migration",
                       "ameboidal-type cell migration",
                       "tissue migration",
                       "external encapsulating structure organization",
                       "epithelial cell proliferation",
                       "negative regulation of response to external stimulus",
                       "extracellular matrix organization",
                       "developmental growth involved in morphogenesis",
                       "negative regulation of cell development",
                       "response to mechanical stimulus",
                       "potassium ion transport"
)

custom_colors <- c("#F21A00","#501E4D","#E1AF00","#B31657" , "#78B7C5", "#3B9AB2")

treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 9, 
         group_color = custom_colors,
         hclust_method = "average", cluster.params = list(n = 6))

fig4b <- treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 9, 
                  group_color = custom_colors,
                  hclust_method = "average", cluster.params = list(n = 6))

ggsave("UPDATE_fig4b.png", plot = fig4b, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)


selected_pathways <- c("axon guidance",
                       "regulation of nervous system development",
                       "axonogenesis",
                       "axon development",
                       "regulation of neuron projection development",
                       "neuron projection guidance",
                       "glial cell differentiation",
                       "regulation of trans-synaptic signaling",
                       "neuron recognition",
                       "gliogenesis",
                       "negative regulation of neuron projection development",
                       "positive regulation of neurogenesis",
                       "synapse organization",
                       "nerve development",
                       "modulation of chemical synaptic transmission",
                       "positive regulation of nervous system development",
                       "regulation of neurogenesis"
)



custom_colors <- c("#B31657","#501E4D","#F21A00","#78B7C5", "#3B9AB2","#E1AF00")

treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 9, 
         group_color = custom_colors,
         hclust_method = "average", cluster.params = list(n = 6))

fig4c <- treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 9, 
                  group_color = custom_colors,
                  hclust_method = "average", cluster.params = list(n = 6))

ggsave("UPDATE_fig4c.png", plot = fig4c, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)



############################
############################
############################
############################

# Testing new figure #

############################
############################
############################
############################

selected_pathways <- c("kidney development",
                       "renal system development",
                       "mesenchymal cell differentiation",
                       "mesenchyme development",
                       "endothelial cell migration",
                       "regulation of angiogenesis",
                       "muscle cell differentiation",
                       "axon guidance",
                       "regulation of neuron projection development",
                       "glial cell differentiation",
                       "neuron recognition",
                       "gliogenesis",
                       "synapse organization",
                       "modulation of chemical synaptic transmission",
                       "regulation of neurogenesis",
                       "extracellular matrix organization",
                       "positive regulation of MAPK cascade",
                       "ERK1 and ERK2 cascade"
                        )

options(enrichplot.colours = c("#4F1E4D","#A11B58", "#E03F44","#F7C6A6"))



custom_colors <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#4F1E4D")



custom_colors <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#4F1E4D")


custom_colors <- c("#fee090","#d73027","#f46d43","#4575b4","#a50026","#fdae61","#74add1","#abd9e9","#e0f3f8","#313695")

custom_colors <- c("#fdae61","#fee090","#d73027","#4575b4","#b8401e","#A11B60","#74add1","#abd9e9","#e0f3f8","#a50030")



treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 8, 
         group_color = custom_colors, fontsize = 4,
         hclust_method = "average", cluster.params = list(n = 10),
         clusterPanel.params = list(clusterPanel = "heatMap", pie = "equal", legend_n = 3,
                                    colnames_angle = 45),
         offset = 10,
         hilight = TRUE,
         align = "left",
         hexpand = 0.1)


treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 8, 
         group_color = custom_colors,          
         hclust_method = "average", cluster.params = list(n = 10),
         clusterPanel.params = list(clusterPanel = "heatMap", pie = "equal", legend_n = 3,
                                    colnames_angle = 45),
         offset = 10,
         hilight = TRUE,
         align = "left",
         cex_label_category = 20,
         label_format_cladelab = 30,
         cex_category = 20,
         hexpand = 0.1) + theme(text=element_text(size=20))

p1 <- treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 8, 
               group_color = custom_colors, fontsize = 4,
               hclust_method = "average", cluster.params = list(n = 10),
               clusterPanel.params = list(clusterPanel = "heatMap", pie = "equal", legend_n = 3,
                                          colnames_angle = 45),
               offset = 20,
               hilight = TRUE,
               align = "left",
               hexpand = 0.1)


p1 <- treeplot(xxGOBP2, showCategory = selected_pathways, offset_tiplab = 8, 
               group_color = custom_colors, fontsize = 4,
               hclust_method = "average", cluster.params = list(n = 10),
               clusterPanel.params = list(clusterPanel = "heatMap"),
               offset = 20,
               hilight = TRUE,
               align = "left",
               hexpand = 0.1)
p1

p1$layers[[7]]$aes_params$size <- 6
p1

#width=1850&height=1016
ggsave("UPDATE_fig4_final.png", plot = p1, width = 18.5, height = 10.16, dpi = 600, limitsize = FALSE)



###############################################
###############################################
###############################################
###############################################
###############################################


##############
# Human KEGG #
##############


setwd ("~/Desktop/CS_Paper/Todos_datasets_originais") 

full_KEGG <- readRDS("full_KEGG_background.rds")

setwd ("~/Desktop/CS_Paper")

a <- as.data.frame(unique(full_KEGG@compareClusterResult$Cluster))

compareCluster_data <- as.data.frame(full_KEGG)

desired_order <- c("CSB_CS1AN_2014_Egly_Fibroblast_CS1AN", "CSB_CS1AN_Wang_KEGG",
                   "CSB_CS1AN_Okur_KEGG",
                   "CSB_iPSC_Liu_KEGG",
                   "CSB_MSC_Liu_KEGG", "CSB_NeuroB_Wang_KEGG", "CSA_NeuroB_Wang_KEGG",
                   "CS_Patients_KEGG_Wang_2014", "CSA_Mouse_1m_KEGG",
                   "CSA_Mouse_12m_KEGG", "CSA_Mouse_24m_KEGG")

# Reorder the data frame according to the desired cluster names
compareCluster_data$Cluster <- factor(compareCluster_data$Cluster, levels = desired_order)
compareCluster_data <- compareCluster_data[order(compareCluster_data$Cluster), ]

# Update the compareClusterResult object with the reordered data
full_KEGG@compareClusterResult <- compareCluster_data

compareCluster_data <- as.data.frame(full_KEGG)

# Filter out clusters that start with "CSA_Mouse"
filtered_data <- compareCluster_data[!grepl("^CSA_Mouse", compareCluster_data$Cluster), ]

# Update the compareClusterResult object with the filtered data
full_KEGG@compareClusterResult <- filtered_data


#full_KEGG <- readRDS("full_human_KEGG_background.rds")

df_KEGG <- as.data.frame(full_KEGG)

df <- as.data.frame(unique(df_KEGG$Description))

subset_mouse <- df_KEGG[grep("^CSA_M", df_KEGG$Cluster), ]
unique_subset_mouse <- as.data.frame(unique(subset_mouse$Description))

subset_human <- subset(df_KEGG, !grepl("CSA_M", Cluster))
unique_subset_human <- as.data.frame(unique(subset_human$Description))

#write.csv(df_KEGG, "mouse_human_KEGG_background.csv")


unique_values <- as.data.frame(unique(df_KEGG$Cluster))


Human_df_KEGG <- df_KEGG %>%
  filter(!grepl("Mouse", Cluster))

unique_values <- as.data.frame(unique(Human_df_KEGG$Cluster))


column_name <- "ID"

#write.csv(result_human_df_KEGG,"human_KEGG_background.csv")


result_df_KEGG <- Human_df_KEGG %>%
  group_by(ID) %>%
  filter(n() >= 6)

column_name <- "ID"


df_no_duplicates <- result_df_KEGG[!duplicated(result_df_KEGG[[column_name]]), ]


#column_name <- "ID"

#write.csv(result_df_KEGG,"all_KEGG_background.csv")

selectKEGG <- data.frame("ID" = full_KEGG@compareClusterResult$ID, "name" = full_KEGG@compareClusterResult$Description)
select_KEGG_terms <- df_no_duplicates$ID
selectKEGG <- selectKEGG$ID[!selectKEGG$ID %in% select_KEGG_terms]
granules_KEGG <- dropGO(full_KEGG, term = selectKEGG)
granules_KEGG <- dropGO(granules_KEGG)
g <- as.data.frame(granules_KEGG)

#xxKEGG <- pairwise_termsim(full_KEGG)  


xxKEGG2 <- pairwise_termsim(granules_KEGG)  


x <- as.data.frame(xxKEGG2)



custom_colors <- c("#501E4D", "#B31657" ,"#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")


# 1850 : 1016
custom_colors <- c("#501E4D", "#B31657" ,"#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")


fig4b <- cnetplot(xxKEGG2, showCategory = 2, layout = "kk",
                  pie.params = list(legend_n = 2, pie = "Count"),
                  node_label="category", cex_label_category=1.9, cex_gene=1.2)

fig4b +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)

fig4B <- fig4b +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)


ggsave("Fig4B.png", plot = fig4B, width = 18.50, height = 10.16, units = "in", dpi = 600)






cnetplot(xxKEGG2)

cnetplot(xxKEGG2, showCategory = 2, layout = "kk",node_label="category", 
         pie.params = list(legend_n = 2, pie = "Count"),
         cex_label_category=1.9, cex_gene=1.2)

cnetplot(xxKEGG2, showCategory = 2, layout = "kk",node_label="none", color_category = "blue",color_item = custom_colors,
         pie.params = list(legend_n = 2, pie = "Count"),
         cex_label_category=1.9, cex_gene=1.2)
