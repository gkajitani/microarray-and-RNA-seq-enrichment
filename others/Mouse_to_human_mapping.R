library (tidyverse)
require (biomaRt)

setwd ("~/Desktop/CS_Paper/Todos_datasets_originais")

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}

df <- read.csv("All_genes_protein.csv")

a1 <- convert_mouse_to_human(df$CSA_mouse_01_m_AllGenes)
a2 <- convert_mouse_to_human(df$CSA_mouse_12_m_AllGenes)
a3 <- convert_mouse_to_human(df$CSA_mouse_24_m_AllGenes)

a1_2 <- as.data.frame(a1)
a2_2 <- as.data.frame(a2)
a3_2 <- as.data.frame(a3)

write.csv (a1_2, "CSA_mouse_01_m_human.csv")
write.csv (a2_2, "CSA_mouse_12_m_human.csv")
write.csv (a3_2, "CSA_mouse_24_m_human.csv")

##############
# Background #

background <- read.csv("background_all_mouse.csv")
background_un <- unique(as.data.frame(background$Gene))

back_human <- convert_mouse_to_human(background_un$`background$Gene`)
back_human_2 <- as.data.frame(back_human)
write.csv(back_human_2, "background_mouse_to_human.csv")

####################


human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

getLDS(attributes = c("external_gene_name"),
       filters = "external_gene_name", values = c("Tp53", "Stat1", "Pdcd1"), mart = mouse,
       attributesL = c("external_gene_name"), martL = human)

getLDS(attributes = c("external_gene_name"),
       filters = "external_gene_name", values = c("TP53", "MIR155HG", "STAT1", "PDCD1"), mart = human,
       attributesL = c("external_gene_name"), martL = mouse)

?useEnsembl

a <- getLDS(attributes = c("external_gene_name"),
            filters = "external_gene_name", values = c("Tp53", "Stat1", "Pdcd1"), mart = mouse,
            attributesL = c("external_gene_name"), martL = human)
