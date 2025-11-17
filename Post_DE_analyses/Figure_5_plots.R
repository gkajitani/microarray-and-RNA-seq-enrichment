library (dplyr)
library (ggplot2)
library (grid)
library (clusterProfiler)
library (enrichplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#library(DOSE)
library (tidydr)
library (ggtreeExtra)
library(GOSemSim)
library (multienrichjam)
library (wesanderson)

#options(enrichplot.colours = c("#571845", "#900C3E", "#C70039", "#FF5733","#FFC300"))

options(enrichplot.colours = c("#501E4D","#791F59", "#9F1A5B", "#B31657","#DE2E43"))

#options(enrichplot.colours = c("#501E4D","#791F59", "#B11758", "#EB493E","#F7C6A6"))

#options(enrichplot.colours = wes_palette("Zissou1"))

#options(enrichplot.colours = c("#F21A00", "#E1AF00", "#EBCC2A"))


######################################## ALL true brain GOBP 

setwd ("~/Desktop/CS_Paper/GSEA/true_brain")

full_GOBP <-  readRDS("GSEA_GOBP_brain.rds")

df_GOBP <- as.data.frame(full_GOBP)

unique_clusters <- unique(df_GOBP$Cluster)

# Create a list to store subsets of the data frame
subsets <- list()

# Subset the data frame by each unique element in the "Cluster" column
for (cluster in unique_clusters) {
  subset_df <- df_GOBP %>% filter(Cluster == cluster)
  assign(paste0("df_", cluster), subset_df)
}

write.csv(df_CS_Patients_GOBP_Wang_2014,"df_CS_Patients_GOBP_Wang_2014.csv")
write.csv(df_CSA_Mouse_12m_GOBP, "df_CSA_Mouse_12m_GOBP.csv")
write.csv(df_CSA_Mouse_1m_GOBP, "df_CSA_Mouse_1m_GOBP.csv")
write.csv(df_CSA_Mouse_24m_GOBP, "df_CSA_Mouse_24m_GOBP.csv")



compareCluster_data <- as.data.frame(full_GOBP)


unique_datasets <- as.data.frame(unique(df_GOBP$Cluster))

desired_order <- c("CS_Patients_GOBP_Wang_2014",
                   "CSA_Mouse_1m_GOBP",
                   "CSA_Mouse_12m_GOBP", 
                   "CSA_Mouse_24m_GOBP")

# Reorder the data frame according to the desired cluster names
compareCluster_data$Cluster <- factor(compareCluster_data$Cluster, levels = desired_order)
compareCluster_data <- compareCluster_data[order(compareCluster_data$Cluster), ]

# Update the compareClusterResult object with the reordered data
full_GOBP@compareClusterResult <- compareCluster_data


df_GOBP <- as.data.frame(full_GOBP)

write.csv(df_GOBP, "full_GOBP_brain.csv")

unique_terms <- as.data.frame(unique(df_GOBP$Description))


unique_datasets <- as.data.frame(unique(df_GOBP$Cluster))

result_df_GOBP <- df_GOBP %>%
  group_by(ID) %>%
  filter(n() >= 4)

column_name <- "ID"



# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_GOBP[!duplicated(result_df_GOBP[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueGOBP_4.csv")

selectGOBP <- data.frame("ID" = full_GOBP@compareClusterResult$ID, "name" = full_GOBP@compareClusterResult$Description)
select_GOBP_terms <- df_no_duplicates$ID
selectGOBP <- selectGOBP$ID[!selectGOBP$ID %in% select_GOBP_terms]
granules_GOBP <- dropGO(full_GOBP, term = selectGOBP)
g <- as.data.frame(granules_GOBP)

write.csv (g,"GOBP_brain_intersection.csv")

#xxGOBP <- pairwise_termsim(full_GOBP)  

xxGOBP2 <- pairwise_termsim(granules_GOBP)  

x <- as.data.frame(xxGOBP2)

dotplot(xxGOBP2)

treeplot(xxGOBP2, clusterPanel.params = list(colnames_angle = 45) )

treeplot(xxGOBP2, clusterPanel.params = list(colnames_angle = 45), color="NES" )

selected_pathways <- c( "vesicle-mediated transport in synapse",
                        "dendrite morphogenesis",
                        "granulocyte chemotaxis",
                        "positive regulation of cytokine production",
                        "type II interferon production",
                        "chemokine-mediated signaling pathway",
                        "cellular response to chemokine",
                        "adaptive immune response",
                        "cell killing",
                        "cellular response to interleukin-1"
)

all_pathways <- c("morphogenesis of an epithelium",
                  "epithelial cell proliferation",
                  "morphogenesis of a branching epithelium",
                  "positive regulation of cell motility",
                  "morphogenesis of a branching structure",
                  "vesicle-mediated transport in synapse",
                  "phosphatidylinositol 3-kinase/protein kinase B signal transduction",
                  "positive regulation of cytokine production",
                  "type II interferon production",
                  "chemokine-mediated signaling pathway",
                  "response to chemokine",
                  "cellular response to chemokine",
                  "adaptive immune response",
                  "granulocyte migration",
                  "neutrophil migration",
                  "granulocyte chemotaxis",
                  "neutrophil chemotaxis",
)



#zissou
custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
custom_colors <- c("#F21A00", "#EBCC2A", "#78B7C5", "#3B9AB2")


custom_colors <- c("#E1AF00","#E03F44","#4F1E4D","#A11B58","#997700")


custom_colors <- c("#4F1E4D","#A11B58","#E03F44")


treeplot(xxGOBP2,showCategory=selected_pathways,cluster.params = list(n = 3), 
         group_color = custom_colors,
         offset_tiplab = 6.1, hclust_method = "average") 


t2 <- treeplot(xxGOBP2,showCategory=selected_pathways,cluster.params = list(n = 3), 
               group_color = custom_colors,
               offset_tiplab = 6.1, hclust_method = "average") 


#width=1850&height=1016

ggsave("fig5c2.png", plot = t2, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)


########### GOMF


full_GOMF <- readRDS("GSEA_GOMF_brain.rds")

df_GOMF <- as.data.frame(full_GOMF)

unique_clusters <- unique(df_GOMF$Cluster)

# Create a list to store subsets of the data frame
subsets <- list()

# Subset the data frame by each unique element in the "Cluster" column
for (cluster in unique_clusters) {
  subset_df <- df_GOMF %>% filter(Cluster == cluster)
  assign(paste0("df_", cluster), subset_df)
}

write.csv(df_CS_Patients_GOMF_Wang_2014,"df_CS_Patients_GOMF_Wang_2014.csv")
write.csv(df_CSA_Mouse_12m_GOMF, "df_CSA_Mouse_12m_GOMF.csv")
write.csv(df_CSA_Mouse_1m_GOMF, "df_CSA_Mouse_1m_GOMF.csv")
write.csv(df_CSA_Mouse_24m_GOMF, "df_CSA_Mouse_24m_GOMF.csv")

unique_terms <- as.data.frame(unique(df_GOMF$Description))

write.csv(df_GOMF, "full_GOMF_brain.csv")

result_df_GOMF <- df_GOMF %>%
  group_by(ID) %>%
  filter(n() >= 4)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_GOMF[!duplicated(result_df_GOMF[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueGOMF_4.csv")

selectGOMF <- data.frame("ID" = full_GOMF@compareClusterResult$ID, "name" = full_GOMF@compareClusterResult$Description)
select_GOMF_terms <- df_no_duplicates$ID
selectGOMF <- selectGOMF$ID[!selectGOMF$ID %in% select_GOMF_terms]
granules_GOMF <- dropGO(full_GOMF, term = selectGOMF)
g <- as.data.frame(granules_GOMF)

write.csv (g,"GOMF_brain_intersection.csv")

#xxGOMF <- pairwise_termsim(full_GOMF)  

xxGOMF2 <- pairwise_termsim(granules_GOMF)  

x <- as.data.frame(xxGOMF2)

dotplot(xxGOMF2)

treeplot(xxGOMF2)

treeplot(xxGOMF2, showCategory = 20, cluster.params = list(n = 7))


########### GOCC



full_GOCC <- readRDS("GSEA_GOCC_brain.rds")

df_GOCC <- as.data.frame(full_GOCC)

unique_clusters <- unique(df_GOCC$Cluster)

# Create a list to store subsets of the data frame
subsets <- list()

# Subset the data frame by each unique element in the "Cluster" column
for (cluster in unique_clusters) {
  subset_df <- df_GOCC %>% filter(Cluster == cluster)
  assign(paste0("df_", cluster), subset_df)
}

write.csv(df_CS_Patients_GOCC_Wang_2014,"df_CS_Patients_GOCC_Wang_2014.csv")
write.csv(df_CSA_Mouse_12m_GOCC, "df_CSA_Mouse_12m_GOCC.csv")
write.csv(df_CSA_Mouse_1m_GOCC, "df_CSA_Mouse_1m_GOCC.csv")
write.csv(df_CSA_Mouse_24m_GOCC, "df_CSA_Mouse_24m_GOCC.csv")

write.csv(df_GOCC, "full_GOCC_brain.csv")


result_df_GOCC <- df_GOCC %>%
  group_by(ID) %>%
  filter(n() >= 4)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_GOCC[!duplicated(result_df_GOCC[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueGOCC_5.csv")

selectGOCC <- data.frame("ID" = full_GOCC@compareClusterResult$ID, "name" = full_GOCC@compareClusterResult$Description)
select_GOCC_terms <- df_no_duplicates$ID
selectGOCC <- selectGOCC$ID[!selectGOCC$ID %in% select_GOCC_terms]
granules_GOCC <- dropGO(full_GOCC, term = selectGOCC)
g <- as.data.frame(granules_GOCC)

write.csv (g,"GOCC_brain_intersection.csv")



#xxGOCC <- pairwise_termsim(full_GOCC)  

xxGOCC2 <- pairwise_termsim(granules_GOCC)  

x <- as.data.frame(xxGOCC2)

write.csv(x,"GOCC_all.csv")

dotplot(xxGOCC2)

emapplot(xxGOCC2)

treeplot(xxGOCC2)



all_pathways <- c("morphogenesis of an epithelium",
                  "epithelial cell proliferation",
                  "morphogenesis of a branching epithelium",
                  "positive regulation of cell motility",
                  "morphogenesis of a branching structure",
                  "vesicle-mediated transport in synapse",
                  "phosphatidylinositol 3-kinase/protein kinase B signal transduction",
                  "positive regulation of cytokine production",
                  "type II interferon production",
                  "chemokine-mediated signaling pathway",
                  "response to chemokine",
                  "cellular response to chemokine",
                  "adaptive immune response",
                  "granulocyte migration",
                  "neutrophil migration",
                  "granulocyte chemotaxis",
                  "neutrophil chemotaxis",
)



selected_pathways <- c("postsynaptic specialization membrane",
                       "postsynaptic density membrane",
                       "GABA-ergic synapse",
                       "postsynaptic membrane",
                       "neuron to neuron synapse",
                       "postsynaptic density",
                       "postsynaptic specialization",
                       "asymmetric synapse"
)

selected_pathways <- c("GABA-ergic synapse",
                       "postsynaptic membrane",
                       "neuron to neuron synapse",
                       "postsynaptic density",
                       "postsynaptic specialization",
                       "asymmetric synapse"
)

custom_colors <- c("#E1AF00","#E03F44","#4F1E4D")


custom_colors <- c("#997700","#E1AF00")


treeplot(xxGOCC2,showCategory=selected_pathways,cluster.params = list(n = 2), 
         group_color = custom_colors,
         offset_tiplab = 6.1, hclust_method = "average") 


t3 <- treeplot(xxGOCC2,showCategory=selected_pathways,cluster.params = list(n = 2), 
               group_color = custom_colors,
               offset_tiplab = 6.1, hclust_method = "average") 


#width=1850&height=1016

ggsave("fig5c3.png", plot = t3, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)


treeplot(xxGOCC2,cluster.params = list(n = 3), 
         offset_tiplab = 6.6, hclust_method = "average") 


##### GO ALL


full_GOALL <- readRDS("GSEA_GOALL_brain.rds")

df_GOALL <- as.data.frame(full_GOALL)

result_df_GOALL <- df_GOALL %>%
  group_by(ID) %>%
  filter(n() >= 4)

column_name <- "ID"

# Keep only the first oALLurrence of each value in the specified column
df_no_duplicates <- result_df_GOALL[!duplicated(result_df_GOALL[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueGOALL.csv")

selectGOALL <- data.frame("ID" = full_GOALL@compareClusterResult$ID, "name" = full_GOALL@compareClusterResult$Description)
select_GOALL_terms <- df_no_duplicates$ID
selectGOALL <- selectGOALL$ID[!selectGOALL$ID %in% select_GOALL_terms]
granules_GOALL <- dropGO(full_GOALL, term = selectGOALL)
g <- as.data.frame(granules_GOALL)

write.csv (g,"GOALL_brain_intersection.csv")


#xxGOALL <- pairwise_termsim(full_GOALL)  

xxGOALL2 <- pairwise_termsim(granules_GOALL)  






#### KEGG

full_KEGG <- readRDS("GSEA_KEGG_brain.rds")



df_KEGG <- as.data.frame(full_KEGG)


write.csv(df_KEGG, "full_KEGG_brain.csv")


############ Ver se da pra editar os Descriptions e IDs aqui

full_KEGG@compareClusterResult$ID <- gsub("rno", "hsa", full_KEGG@compareClusterResult$ID)
full_KEGG@compareClusterResult$ID <- gsub("mmu", "hsa", full_KEGG@compareClusterResult$ID)
full_KEGG@compareClusterResult$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", full_KEGG@compareClusterResult$Description)



df_KEGG <- as.data.frame(full_KEGG)



unique_clusters <- unique(df_KEGG$Cluster)

# Create a list to store subsets of the data frame
subsets <- list()

# Subset the data frame by each unique element in the "Cluster" column
for (cluster in unique_clusters) {
  subset_df <- df_KEGG %>% filter(Cluster == cluster)
  assign(paste0("df_", cluster), subset_df)
}

write.csv(df_CS_Patients_KEGG_Wang_2014,"df_CS_Patients_KEGG_Wang_2014.csv")
write.csv(df_CSA_Mouse_12m_KEGG, "df_CSA_Mouse_12m_KEGG.csv")
write.csv(df_CSA_Mouse_1m_KEGG, "df_CSA_Mouse_1m_KEGG.csv")
write.csv(df_CSA_Mouse_24m_KEGG, "df_CSA_Mouse_24m_KEGG.csv")

unique_terms <- as.data.frame(unique(df_KEGG$Description))


result_df_KEGG <- df_KEGG %>%
  group_by(ID) %>%
  filter(n() >= 4)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_KEGG[!duplicated(result_df_KEGG[[column_name]]), ]

write.csv(df_no_duplicates, "uniqueKEGG_all.csv")

selectKEGG <- data.frame("ID" = full_KEGG@compareClusterResult$ID, "name" = full_KEGG@compareClusterResult$Description)
select_KEGG_terms <- df_no_duplicates$ID
selectKEGG <- selectKEGG$ID[!selectKEGG$ID %in% select_KEGG_terms]
granules_KEGG <- dropGO(full_KEGG, term = selectKEGG)
g <- as.data.frame(granules_KEGG)

write.csv (g,"KEGG_brain_intersection.csv")


#xxKEGG <- pairwise_termsim(full_KEGG)  

xxKEGG2 <- pairwise_termsim(granules_KEGG)  

x <- as.data.frame(xxKEGG2)



dotplot(xxKEGG2)

emapplot(xxKEGG2)

cnetplot(xxKEGG2)

treeplot(xxKEGG2)

treeplot(xxKEGG2, showCategory = 20, cluster.params = list(n = 7))


########### Reactome

full_Pathway <- readRDS("GSEA_Pathway_brain.rds")

df_Pathway <- as.data.frame(full_Pathway)


write.csv(df_Pathway, "full_Reactome_brain.csv")

full_Pathway@compareClusterResult$ID <- gsub("RNA", "HSA", full_Pathway@compareClusterResult$ID)
full_Pathway@compareClusterResult$ID <- gsub("MMU", "HSA", full_Pathway@compareClusterResult$ID)

df_Pathway <- as.data.frame(full_Pathway)

unique_terms <- as.data.frame(unique(df_Pathway$Description))


unique_clusters <- unique(df_Pathway$Cluster)

# Create a list to store subsets of the data frame
subsets <- list()

# Subset the data frame by each unique element in the "Cluster" column
for (cluster in unique_clusters) {
  subset_df <- df_Pathway %>% filter(Cluster == cluster)
  assign(paste0("df_", cluster), subset_df)
}

write.csv(df_CS_Patients_Pathway_Wang_2014,"df_CS_Patients_Pathway_Wang_2014.csv")
write.csv(df_CSA_Mouse_12m_Pathway, "df_CSA_Mouse_12m_Pathway.csv")
write.csv(df_CSA_Mouse_1m_Pathway, "df_CSA_Mouse_1m_Pathway.csv")
write.csv(df_CSA_Mouse_24m_Pathway, "df_CSA_Mouse_24m_Pathway.csv")



result_df_Pathway <- df_Pathway %>%
  group_by(ID) %>%
  filter(n() >= 4)

column_name <- "ID"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- result_df_Pathway[!duplicated(result_df_Pathway[[column_name]]), ]

write.csv(df_no_duplicates, "uniquePathway_4.csv")

selectPathway <- data.frame("ID" = full_Pathway@compareClusterResult$ID, "name" = full_Pathway@compareClusterResult$Description)
select_Pathway_terms <- df_no_duplicates$ID
selectPathway <- selectPathway$ID[!selectPathway$ID %in% select_Pathway_terms]
granules_Pathway <- dropGO(full_Pathway, term = selectPathway)
g <- as.data.frame(granules_Pathway)

write.csv (g,"Reactome_brain_intersection.csv")


#xxPathway <- pairwise_termsim(full_Pathway)  

xxPathway2 <- pairwise_termsim(granules_Pathway)  

x <- as.data.frame(xxPathway2)

dotplot(xxPathway2)

treeplot(xxPathway2, hclust_method = "average")

treeplot(xxPathway2, showCategory = 20, cluster.params = list(n = 7))



##########################


full_GOALL_brain <- merge_result (list("GOBP"=xxGOBP2,
                                       "GOCC"=xxGOCC2
))

xxGOALL <- pairwise_termsim(full_GOALL_brain)  

xxGOALL@compareClusterResult <- xxGOALL@compareClusterResult[ , -1]

teste <- as.data.frame(xxGOALL)



treeplot(xxGOALL, showCategory = 20, cluster.params = list(n = 7))

selected_pathways <- c( "vesicle-mediated transport in synapse",
                        "dendrite morphogenesis",
                        "granulocyte chemotaxis",
                        "positive regulation of cytokine production",
                        "type II interferon production",
                        "chemokine-mediated signaling pathway",
                        "cellular response to chemokine",
                        "adaptive immune response",
                        "cell killing",
                        "cellular response to interleukin-1",
                        "GABA-ergic synapse",
                        "postsynaptic membrane",
                        "neuron to neuron synapse",
                        "postsynaptic density",
                        "postsynaptic specialization",
                        "asymmetric synapse"
)


custom_colors <- c("#E1AF00","#E03F44","#4F1E4D")


custom_colors <- c("#997700","#E1AF00")


custom_colors <- c("#E1AF00","#E03F44","#4F1E4D","#A11B58")


custom_colors <- c("#E03F44","#A11B58","#4F1E4D","#997700","#E1AF00")

options(enrichplot.colours = c("#501E4D",, "#9F1A5B", "#B31657","#DE2E43"))

custom_colors <- c("#4F1E4D","#791F59","#B31657","#E03F44","#E1AF00")


treeplot(xxGOALL,showCategory=selected_pathways,cluster.params = list(n = 5), 
         group_color = custom_colors,
         offset_tiplab = 10, hclust_method = "average") 

t4 <- treeplot(xxGOALL,showCategory=selected_pathways,cluster.params = list(n = 5), 
               group_color = custom_colors,
               offset_tiplab = 10, hclust_method = "average") 



#width=1850&height=1016

ggsave("fig5c4.png", plot = t4, width = 19, height = 10.16, dpi = 600, limitsize = FALSE)
