#################################################################################
# Used to get unique protein coding genes (min pval) from autosomal chromosomes #
#################################################################################

library("tidyverse")
require (biomaRt)

setwd("~/Desktop/CS_Paper/Todos_datasets_originais/Updated_microarray/old_young_csb_wt")

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), 
                        filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)

#######################################


df <- read.csv("Csb_old_young.csv")
filtered_df <- df[df$Gene.symbol != "", ]
filtered_df$P.Value <- as.numeric(filtered_df$P.Value)
unique_df2 <- as.data.frame(unique(filtered_df$Gene.symbol))


unique_df <- filtered_df[order(filtered_df$Gene.symbol, filtered_df$P.Value), ]
unique_df <- unique_df[!duplicated(unique_df$Gene.symbol), ]

merged_data <- merge(genes, unique_df, by.x = "external_gene_name", by.y = "Gene.symbol")

write.csv(merged_data,"unique_Csb_old_young.csv")

#######################################
#######################################
#######################################

# Some required further processing, such as removing multiple names


library (dplyr)

setwd ("~/Desktop/CS_Paper/Todos_datasets_originais/Updated_microarray")

genes <- read.csv ("human_protein_coding_genes.csv")

GSE58068 <- read.csv ("gene_all_GSE58068_CSB_CS1AN.csv")

df <- GSE58068[complete.cases(GSE58068$Gene.symbol), ]

df <- subset(df, Gene.symbol != "")

setwd ("~/Desktop/CS_Paper/Todos_datasets_originais/Updated_microarray/GSE58068")

write.csv(df, "all_GSE58068_CSB_CS1AN.csv")

df$Gene.symbol <- sapply(strsplit(as.character(df$Gene.symbol), "///"), `[`, 1)

unique <- as.data.frame(unique(df$Gene.symbol))

df <- df %>%
  group_by(Gene.symbol) %>%
  filter(adj.P.Val == max(adj.P.Val)) %>%
  distinct(Gene.symbol, .keep_all = TRUE) %>%
  ungroup()

merged_data <- merge(genes, df, by.x = "external_gene_name", by.y = "Gene.symbol")

write.csv(merged_data ,"GSE58068_protein_coding_all.csv")
