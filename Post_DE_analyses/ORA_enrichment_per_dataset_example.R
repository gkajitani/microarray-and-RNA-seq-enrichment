library (dplyr)
library (ggplot2)
library (grid)
library (clusterProfiler)
library (enrichplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library (tidydr)
library (ggtreeExtra)
library(GOSemSim)
library (multienrichjam)

setwd ("~/Desktop/CS_Paper/") 

# "All_genes_protein.csv" contained columns of differentially expressed protein-coding genes (filtered using biomaRt) per dataset;
df <- read.csv("All_genes_protein.csv")

# "All_genes_protein_background.csv" contained columns of all detected proten-coding genes (filtered using biomaRt) per dataset;
df_backgroud <- read.csv("All_genes_protein_background.csv")

######################
# Biological Process #
######################

CS_Patients_GOBP_Wang_2014 <- enrichGO(gene         = df$CS_Patients_2014_Wang_All,
                                       universe = df_backgroud$CS_Patients_2014_Wang_All,
                                       OrgDb         = org.Hs.eg.db,
                                       keyType       = 'SYMBOL',
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.05,
                                       readable=TRUE)


CS_P <- as.data.frame(CS_Patients_GOBP_Wang_2014)



CSA_NeuroB_Wang_GOBP <- enrichGO(gene         = df$CSA_NeuroB_2019_Wang,
                                 universe = df_backgroud$CSA_NeuroB_2019_Wang,
                                 OrgDb         = org.Hs.eg.db,
                                 keyType       = 'SYMBOL',
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05)


CS_A4 <- as.data.frame(CSA_NeuroB_Wang_GOBP)
# Enriqueceu pra ERK1/2


CSB_CS1AN_Wang_GOBP <- enrichGO(gene         = df$CSB_CS1AN_2014_Wang,
                                universe = df_backgroud$CSB_CS1AN_2014_Wang,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = 'SYMBOL',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05)

CS_B1 <- as.data.frame(CSB_CS1AN_Wang_GOBP)


CSB_CS1AN_Okur_GOBP <- enrichGO(gene         = df$CSB_CS1AN_2020_OkurA,
                                universe = df_backgroud$CSB_CS1AN_2020_OkurA,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = 'SYMBOL',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05)

CS_B2 <- as.data.frame(CSB_CS1AN_Okur_GOBP)

CSB_iPSC_Liu_GOBP <- enrichGO(gene         = df$CSB_Liu_iPSC_2018,
                              universe = df_backgroud$CSB_Liu_iPSC_2018,
                              OrgDb         = org.Hs.eg.db,
                              keyType       = 'SYMBOL',
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.05)


CS_B3 <- as.data.frame(CSB_iPSC_Liu_GOBP)



CSB_MSC_Liu_GOBP <- enrichGO(gene         = df$CSB_Liu_MSC_2018,
                             universe = df_backgroud$CSB_Liu_MSC_2018,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = 'SYMBOL',
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05)

CS_B4 <- as.data.frame(CSB_MSC_Liu_GOBP)


CSB_NeuroB_Wang_GOBP <- enrichGO(gene         = df$CSB_NeuroB_2019_Wang,
                                 universe = df_backgroud$CSB_NeuroB_2019_Wang,
                                 OrgDb         = org.Hs.eg.db,
                                 keyType       = 'SYMBOL',
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05)


CS_B5 <- as.data.frame(CSB_NeuroB_Wang_GOBP)


CSB_CS1AN_Egly_GOBP <- enrichGO(gene         = df$CSB_CS1AN_2014_Egly_Fibroblast_CS1AN,
                                universe = df_backgroud$CSB_CS1AN_2014_Egly_Fibroblast_CS1AN,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = 'SYMBOL',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05)


CS_B6 <- as.data.frame(CSB_CS1AN_Egly_GOBP)


full_GOBP <- merge_result (list("1_CSB_Wang"=CSB_CS1AN_Wang_GOBP,
                                "2_CSB_Okur"=CSB_CS1AN_Okur_GOBP,"3_CSB_iPSC"=CSB_iPSC_Liu_GOBP,
                                "4_CSB_MSC"=CSB_MSC_Liu_GOBP, "5_CSB_NeuroB"=CSB_NeuroB_Wang_GOBP,
                                "6_CS_P"=CS_Patients_GOBP_Wang_2014, "7_CSA_NeuroB" =CSA_NeuroB_Wang_GOBP,
                                "CSB_CS1AN_Egly_GOBP" = CSB_CS1AN_Egly_GOBP
))


saveRDS (full_GOBP,"Humans_GOBP_background2.rds")
