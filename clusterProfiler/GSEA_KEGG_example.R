library (dplyr)
library (ggplot2)
library (grid)
library (clusterProfiler)
library (enrichplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library (tidydr)
library (ggtreeExtra)
library(GOSemSim)
library (multienrichjam)
library(DOSE)


options(enrichplot.colours = c("red","blue"))


############
# Selected #
############


#########
# Human #
#########

setwd ("~/Desktop/CS_Paper/GSEA/human")

organism = "org.Hs.eg.db"


df = read.csv("all_GSE58067_CS_Patients.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$logFC
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CS_Patients_KEGG_Wang_2014 <- gseKEGG(geneList=gene_list, 
                                      exponent = 1,
                                      eps = 0,
                                      minGSSize = 5, 
                                      maxGSSize = 1000, 
                                      pvalueCutoff = 0.05, 
                                      verbose = TRUE, 
                                      organism = "hsa", 
                                      pAdjustMethod = "BH")
CS_Patients_KEGG_Wang_2014 <- DOSE::setReadable(CS_Patients_KEGG_Wang_2014, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CS_P <- as.data.frame(CS_Patients_KEGG_Wang_2014)



df = read.csv("all_CSA_NeuroB.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSA_NeuroB_Wang_KEGG <- gseKEGG(geneList=gene_list, 
                                exponent = 1,
                                eps = 0,
                                minGSSize = 5, 
                                maxGSSize = 1000, 
                                pvalueCutoff = 0.05, 
                                verbose = TRUE, 
                                organism = "hsa", 
                                pAdjustMethod = "BH")
CSA_NeuroB_Wang_KEGG <- DOSE::setReadable(CSA_NeuroB_Wang_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CSA_Neuro <- as.data.frame(CSA_NeuroB_Wang_KEGG)


df = read.csv("all_GSE57923_CSB_CS1AN.csv", header=TRUE)
df <- df %>%
  group_by(Gene.symbol) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$logFC
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSB_CS1AN_Egly_KEGG <- gseKEGG(geneList=gene_list, 
                               exponent = 1,
                               eps = 0,
                               minGSSize = 5, 
                               maxGSSize = 1000, 
                               pvalueCutoff = 0.05, 
                               verbose = TRUE, 
                               organism = "hsa", 
                               pAdjustMethod = "BH")
CSB_CS1AN_Egly_KEGG <- DOSE::setReadable(CSB_CS1AN_Egly_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CSB_CS1AN_Egly <- as.data.frame(CSB_CS1AN_Egly_KEGG)



df = read.csv("all_GSE58068_CSB_CS1AN.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$logFC
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSB_CS1AN_Wang_KEGG <- gseKEGG(geneList=gene_list, 
                               exponent = 1,
                               eps = 0,
                               minGSSize = 5, 
                               maxGSSize = 1000, 
                               pvalueCutoff = 0.05, 
                               verbose = TRUE, 
                               organism = "hsa", 
                               pAdjustMethod = "BH")
CSB_CS1AN_Wang_KEGG <- DOSE::setReadable(CSB_CS1AN_Wang_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CSB_CS1AN_Wang <- as.data.frame(CSB_CS1AN_Wang_KEGG)


df = read.csv("all_CS1AN_Okur.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSB_CS1AN_Okur_KEGG <- gseKEGG(geneList=gene_list, 
                               exponent = 1,
                               eps = 0,
                               minGSSize = 5, 
                               maxGSSize = 1000, 
                               pvalueCutoff = 0.05, 
                               verbose = TRUE, 
                               organism = "hsa", 
                               pAdjustMethod = "BH")
CSB_CS1AN_Okur_KEGG <- DOSE::setReadable(CSB_CS1AN_Okur_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CSB_CS1AN_Okur <- as.data.frame(CSB_CS1AN_Okur_KEGG)



df = read.csv("all_iPSC.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSB_iPSC_Liu_KEGG <- gseKEGG(geneList=gene_list, 
                             exponent = 1,
                             eps = 0,
                             minGSSize = 5, 
                             maxGSSize = 1000, 
                             pvalueCutoff = 0.05, 
                             verbose = TRUE, 
                             organism = "hsa", 
                             pAdjustMethod = "BH")
CSB_iPSC_Liu_KEGG <- DOSE::setReadable(CSB_iPSC_Liu_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CSB_iPSC <- as.data.frame(CSB_iPSC_Liu_KEGG)



df = read.csv("all_MSC.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSB_MSC_Liu_KEGG <- gseKEGG(geneList=gene_list, 
                            exponent = 1,
                            eps = 0,
                            minGSSize = 5, 
                            maxGSSize = 1000, 
                            pvalueCutoff = 0.05, 
                            verbose = TRUE, 
                            organism = "hsa", 
                            pAdjustMethod = "BH")
CSB_MSC_Liu_KEGG <- DOSE::setReadable(CSB_MSC_Liu_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CSB_MSC <- as.data.frame(CSB_MSC_Liu_KEGG)



df = read.csv("all_CSB_NeuroB.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSB_NeuroB_Wang_KEGG <- gseKEGG(geneList=gene_list, 
                                exponent = 1,
                                eps = 0,
                                minGSSize = 5, 
                                maxGSSize = 1000, 
                                pvalueCutoff = 0.05, 
                                verbose = TRUE, 
                                organism = "hsa", 
                                pAdjustMethod = "BH")
CSB_NeuroB_Wang_KEGG <- DOSE::setReadable(CSB_NeuroB_Wang_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CSB_Neuro <- as.data.frame(CSB_NeuroB_Wang_KEGG)



#########
# Mouse #
#########


setwd ("~/Desktop/CS_Paper/GSEA/mouse")

organism = "org.Mm.eg.db"

df = read.csv("all_CSA_1m.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSA_Mouse_1m_KEGG <- gseKEGG(geneList=gene_list, 
                             exponent = 1,
                             eps = 0,
                             minGSSize = 5, 
                             maxGSSize = 1000, 
                             pvalueCutoff = 0.05, 
                             verbose = TRUE, 
                             organism = "mmu", 
                             pAdjustMethod = "BH")
CSA_Mouse_1m_KEGG <- DOSE::setReadable(CSA_Mouse_1m_KEGG, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID') 
CSA_mouse_1m <- as.data.frame(CSA_Mouse_1m_KEGG)


df = read.csv("all_CSA_12m.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSA_Mouse_12m_KEGG <- gseKEGG(geneList=gene_list, 
                              exponent = 1,
                              eps = 0,
                              minGSSize = 5, 
                              maxGSSize = 1000, 
                              pvalueCutoff = 0.05, 
                              verbose = TRUE, 
                              organism = "mmu", 
                              pAdjustMethod = "BH")
CSA_Mouse_12m_KEGG <- DOSE::setReadable(CSA_Mouse_12m_KEGG, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID') 
CSA_mouse_12m<- as.data.frame(CSA_Mouse_12m_KEGG)



df = read.csv("all_CSA_24m.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ENTREZID
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

CSA_Mouse_24m_KEGG <- gseKEGG(geneList=gene_list, 
                              exponent = 1,
                              eps = 0,
                              minGSSize = 5, 
                              maxGSSize = 1000, 
                              pvalueCutoff = 0.05, 
                              verbose = TRUE, 
                              organism = "mmu", 
                              pAdjustMethod = "BH")
CSA_Mouse_24m_KEGG <- DOSE::setReadable(CSA_Mouse_24m_KEGG, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID') 
CSA_mouse_24m <- as.data.frame(CSA_Mouse_24m_KEGG)



############
# Filtered #
############







########### Rat

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/Rat/")


organism = "org.Rn.eg.db"

# reading in data from deseq2
df = read.csv("all_Rat.csv", header=TRUE)
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

rat_gse <- gseKEGG(geneList=gene_list, 
                   exponent = 1,
                   eps = 0,
                   minGSSize = 5, 
                   maxGSSize = 1000, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   organism = "rat", 
                   pAdjustMethod = "BH")
rat_gse <- DOSE::setReadable(rat_gse, OrgDb = 'org.Rn.eg.db', keyType = 'ENTREZID') 
rat <- as.data.frame(rat_gse)

########### Nectin

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/Nectin")

organism = "org.Mm.eg.db"

df <- read.csv ("all_Nectin.csv")
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


mouse1_csb_gse <- gseKEGG(geneList=gene_list, 
                          exponent = 1,
                          eps = 0,
                          minGSSize = 5, 
                          maxGSSize = 1000, 
                          pvalueCutoff = 0.05, 
                          verbose = TRUE, 
                          organism = "mmu", 
                          pAdjustMethod = "BH")
mouse1_csb_gse <- DOSE::setReadable(mouse1_csb_gse, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID') 
mouse_csb_sh <- as.data.frame(mouse1_csb_gse)

########### CSB mouse - adult

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/GSE62194")

organism = "org.Mm.eg.db"

df <- read.csv ("all_GSE62194_adult_CSB.csv")
df <- df %>%
  group_by(Gene.symbol) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

adult_csb_mouse_gse <- gseKEGG(geneList=gene_list, 
                               exponent = 1,
                               eps = 0,
                               minGSSize = 5, 
                               maxGSSize = 1000, 
                               pvalueCutoff = 0.05, 
                               verbose = TRUE, 
                               organism = "mmu", 
                               pAdjustMethod = "BH")
adult_csb_mouse_gse <- DOSE::setReadable(adult_csb_mouse_gse, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID') 
adult_csb_mouse <- as.data.frame(adult_csb_mouse_gse)

########### CSB mouse - young

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/GSE62194")

organism = "org.Mm.eg.db"

df <- read.csv ("all_GSE62194_young_CSB.csv")
df <- df %>%
  group_by(Gene.symbol) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

young_csb_mouse_gse <- gseKEGG(geneList=gene_list, 
                               exponent = 1,
                               eps = 0,
                               minGSSize = 5, 
                               maxGSSize = 1000, 
                               pvalueCutoff = 0.05, 
                               verbose = TRUE, 
                               organism = "mmu", 
                               pAdjustMethod = "BH")
young_csb_mouse_gse <- DOSE::setReadable(young_csb_mouse_gse, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID') 
young_csb_mouse <- as.data.frame(young_csb_mouse_gse)


########### CSB mouse - old

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/GSE62194")

organism = "org.Mm.eg.db"

df <- read.csv ("all_GSE62194_old_CSB.csv")
df <- df %>%
  group_by(Gene.symbol) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

old_csb_mouse_gse <- gseKEGG(geneList=gene_list, 
                             exponent = 1,
                             eps = 0,
                             minGSSize = 5, 
                             maxGSSize = 1000, 
                             pvalueCutoff = 0.05, 
                             verbose = TRUE, 
                             organism = "mmu", 
                             pAdjustMethod = "BH")
old_csb_mouse_gse <- DOSE::setReadable(old_csb_mouse_gse, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID') 
old_csb_mouse <- as.data.frame(old_csb_mouse_gse)


########### GSE8939

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/GSE8939")

organism = "org.Hs.eg.db"

df <- read.csv ("GSE8939.csv")
df$Gene.symbol <- sub("///.*", "", df$Gene.symbol)
df <- df %>%
  group_by(Gene.symbol) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


CS1AN_1_gse <- gseKEGG(geneList=gene_list, 
                       exponent = 1,
                       eps = 0,
                       minGSSize = 5, 
                       maxGSSize = 1000, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       organism = "hsa", 
                       pAdjustMethod = "BH")
CS1AN_1_gse <- DOSE::setReadable(CS1AN_1_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CS1AN_1 <- as.data.frame(CS1AN_1_gse)


########### GSE3407

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/GSE3407")

organism = "org.Hs.eg.db"

df <- read.csv ("GSE3407.csv")
df$Gene.symbol <- sub("///.*", "", df$Gene.symbol)
df <- df %>%
  group_by(Gene.symbol) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


CS1AN_2_gse <- gseKEGG(geneList=gene_list, 
                       exponent = 1,
                       eps = 0,
                       minGSSize = 5, 
                       maxGSSize = 1000, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       organism = "hsa", 
                       pAdjustMethod = "BH")
CS1AN_2_gse <- DOSE::setReadable(CS1AN_2_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 
CS1AN_2 <- as.data.frame(CS1AN_2_gse)



########################


########################

########### C. elegans - csa

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/GSE144556")

organism = "org.Ce.eg.db"

df <- read.csv ("all_GSE144556_CSA_filtered.csv")
df$GENE_SYMBOL <- sub("///.*", "", df$GENE_SYMBOL)
df <- df %>%
  group_by(GENE_SYMBOL) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

a <- bitr(df$GENE_SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Ce.eg.db")

df <- merge(df, a, by.x = "GENE_SYMBOL", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

csa_c_elegans_gse <- gseKEGG(geneList=gene_list, 
                             exponent = 1,
                             eps = 0,
                             keyType =  "ncbi-geneid",
                             minGSSize = 5, 
                             maxGSSize = 1000, 
                             pvalueCutoff = 0.05, 
                             verbose = TRUE, 
                             organism = "cel", 
                             pAdjustMethod = "BH")
csa_c_elegans_gse <- DOSE::setReadable(csa_c_elegans_gse , OrgDb = 'org.Ce.eg.db', keyType = 'ENTREZID') 


csa_c_elegans <- as.data.frame(csa_c_elegans_gse)


########### C. elegans - csb


organism = "org.Ce.eg.db"

df <- read.csv ("all_GSE144556_CSB_filtered.csv")
df$GENE_SYMBOL <- sub("///.*", "", df$GENE_SYMBOL)
df <- df %>%
  group_by(GENE_SYMBOL) %>%
  arrange(
    # First arrange by whether positive or negative
    desc(sign(logFC)),  # positive first
    # Then by absolute value (higher first for positive, lower first for negative)
    if_else(logFC > 0, -logFC, logFC)  # trick to get correct ordering
  ) %>%
  slice_head(n = 1) %>%
  ungroup()
a <- bitr(df$GENE_SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Ce.eg.db")

df <- merge(df, a, by.x = "GENE_SYMBOL", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


csb_c_elegans_gse <- gseKEGG(geneList=gene_list, 
                             exponent = 1,
                             eps = 0,
                             keyType =  "ncbi-geneid",
                             minGSSize = 5, 
                             maxGSSize = 1000, 
                             pvalueCutoff = 0.05, 
                             verbose = TRUE, 
                             organism = "cel", 
                             pAdjustMethod = "BH")
csb_c_elegans_gse <- DOSE::setReadable(csb_c_elegans_gse , OrgDb = 'org.Ce.eg.db', keyType = 'ENTREZID') 


csb_c_elegans <- as.data.frame(csb_c_elegans_gse)

########### GSE124208 - NSC

setwd ("~/Desktop/CS_Paper/GSEA_filtered/Filtered/GSE124208")

organism = "org.Hs.eg.db"

df <- read.csv ("all_NSC.csv")
a <- bitr(df$Gene.symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(df, a, by.x = "Gene.symbol", by.y = "SYMBOL")

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$ENTREZID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


Liu_NSC <- gseKEGG(geneList=gene_list, 
                   exponent = 1,
                   eps = 0,
                   minGSSize = 5, 
                   maxGSSize = 1000, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   organism = "hsa", 
                   pAdjustMethod = "BH")
Liu_NSC <- DOSE::setReadable(Liu_NSC, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID') 

###########################################


#########
# Human #
#########


setwd ("~/Desktop/CS_Paper/GSEA_All/human")



full_KEGG_human <- merge_result (list("CS_Patients_KEGG_Wang_2014"=CS_Patients_KEGG_Wang_2014,
                                      "CSA_NeuroB_Wang_KEGG" =CSA_NeuroB_Wang_KEGG,
                                      "CSB_CS1AN_Wang_KEGG"=CSB_CS1AN_Wang_KEGG, 
                                      "CSB_CS1AN_Okur_KEGG"=CSB_CS1AN_Okur_KEGG,"CSB_iPSC_Liu_KEGG"=CSB_iPSC_Liu_KEGG,
                                      "CSB_MSC_Liu_KEGG"=CSB_MSC_Liu_KEGG, 
                                      "CSB_NSC_Liu_KEGG" = Liu_NSC,
                                      "CSB_NeuroB_Wang_KEGG"=CSB_NeuroB_Wang_KEGG,
                                      "CSB_CS1AN_Egly_KEGG" = CSB_CS1AN_Egly_KEGG,
                                      "CSB_CS1AN_GSE8939"=CS1AN_1_gse, "CSB_CS1AN_3407"=CS1AN_2_gse
))


saveRDS (full_KEGG_human,"GSEA_KEGG_ALL_human_2.rds")



#########
# Mouse #
#########

setwd ("~/Desktop/CS_Paper/GSEA_All/mouse")


full_KEGG_mouse <- merge_result (list("CSA_Mouse_1m_KEGG"=CSA_Mouse_1m_KEGG,
                                "CSA_Mouse_12m_KEGG"=CSA_Mouse_12m_KEGG, "CSA_Mouse_24m_KEGG"=CSA_Mouse_24m_KEGG,
                                "shCSB_mouse"=mouse1_csb_gse,
                                "CSB_adult_mouse"=adult_csb_mouse_gse, "CSB_young_mouse"=young_csb_mouse_gse, 
                                "CSB_old_mouse" =old_csb_mouse_gse
))

saveRDS (full_KEGG_mouse,"GSEA_KEGG_ALL_mouse_2.rds")



#######
# All #
#######


setwd ("~/Desktop/CS_Paper/GSEA_All/all")


full_KEGG <- merge_result (list("CS_Patients_KEGG_Wang_2014"=CS_Patients_KEGG_Wang_2014, "CSA_Mouse_1m_KEGG"=CSA_Mouse_1m_KEGG,
                                "CSA_Mouse_12m_KEGG"=CSA_Mouse_12m_KEGG, "CSA_Mouse_24m_KEGG"=CSA_Mouse_24m_KEGG, "CSA_NeuroB_Wang_KEGG" =CSA_NeuroB_Wang_KEGG,
                                "CSB_CS1AN_Wang_KEGG"=CSB_CS1AN_Wang_KEGG, "CSB_CS1AN_Okur_KEGG"=CSB_CS1AN_Okur_KEGG,"CSB_iPSC_Liu_KEGG"=CSB_iPSC_Liu_KEGG,
                                "CSB_MSC_Liu_KEGG"=CSB_MSC_Liu_KEGG, 
                                "CSB_NSC_Liu_KEGG" = Liu_NSC,
                                "CSB_NeuroB_Wang_KEGG"=CSB_NeuroB_Wang_KEGG,
                                "CSB_CS1AN_Egly_KEGG" = CSB_CS1AN_Egly_KEGG,"CSB_CS1AN_GSE8939"=CS1AN_1_gse, "CSB_CS1AN_3407"=CS1AN_2_gse,
                                "shCSB_mouse"=mouse1_csb_gse,
                                "CSB_adult_mouse"=adult_csb_mouse_gse, "CSB_young_mouse"=young_csb_mouse_gse, 
                                "CSB_old_mouse" =old_csb_mouse_gse, "CSB_Rat"=rat_gse, 
                                "CSB_celegans" = csb_c_elegans_gse, "CSA_celegans" = csa_c_elegans_gse
))


saveRDS (full_KEGG,"GSEA_KEGG_ALL_all_2.rds")




#########
# CS1AN #
#########


setwd ("~/Desktop/CS_Paper/GSEA_All/cs1an/")

full_KEGG <- merge_result (list("CSB_CS1AN_Wang_KEGG"=CSB_CS1AN_Wang_KEGG, "CSB_CS1AN_Okur_KEGG"=CSB_CS1AN_Okur_KEGG,
                                "CSB_CS1AN_Egly_KEGG" = CSB_CS1AN_Egly_KEGG,"CSB_CS1AN_GSE8939"=CS1AN_1_gse, "CSB_CS1AN_3407"=CS1AN_2_gse
))


saveRDS (full_KEGG,"GSEA_KEGG_ALL_CS1AN_2.rds")



################
# Brain/Neuron #
################

setwd ("~/Desktop/CS_Paper/GSEA_All/brain_neuron/")

full_KEGG <- merge_result (list("CS_Patients_KEGG_Wang_2014"=CS_Patients_KEGG_Wang_2014, "CSA_Mouse_1m_KEGG"=CSA_Mouse_1m_KEGG,
                                "CSA_Mouse_12m_KEGG"=CSA_Mouse_12m_KEGG, "CSA_Mouse_24m_KEGG"=CSA_Mouse_24m_KEGG, 
                                "CSA_NeuroB_Wang_KEGG" =CSA_NeuroB_Wang_KEGG,
                                "CSB_NeuroB_Wang_KEGG"=CSB_NeuroB_Wang_KEGG,
                                "CSB_Rat"=rat_gse, "shCSB_mouse"=mouse1_csb_gse,
                                "CSB_adult_mouse"=adult_csb_mouse_gse, "CSB_young_mouse"=young_csb_mouse_gse, 
                                "CSB_old_mouse" =old_csb_mouse_gse
))

saveRDS (full_KEGG,"GSEA_KEGG_ALL_brain_2.rds")



