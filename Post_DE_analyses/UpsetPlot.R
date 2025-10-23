library (UpSetR)
library (tidyverse)
library (ggplot2)
library (grid)
library(data.table)
library(cowplot)

setwd ("~/Desktop/CS_Paper/Todos_datasets_originais") 

AllGenes <- read.csv("All_genes_protein_human_final.csv")

AllGenes[AllGenes == ""] <- NA



AllGenesList <- as.list(AllGenes)


remove_empty_values <- function(x) {
  if (is.list(x)) {
    x <- x[sapply(x, function(y) !is.null(y) && length(y) > 0)]
    x <- lapply(x, remove_empty_values)
  }
  return(x)
}

# Apply the function to remove empty values
AllGenesList <- remove_empty_values(AllGenesList)


# Define a function to replace empty strings with NA
replace_empty_with_na <- function(x) {
  ifelse(x == "", NA, x)
}



AllGenesList <- lapply(AllGenesList, function(sublist) sublist[!is.na(sublist)])


AllGenesList$CS_Patients <- AllGenesList$GSE58067_CS_Patients
AllGenesList$GSE58067_CS_Patients <- NULL

AllGenesList$CSA_SH.SY5Y <- AllGenesList$GSE129002_CSA_SH.SY5Y
AllGenesList$GSE129002_CSA_SH.SY5Y <- NULL

AllGenesList$CSB_SH.SY5Y <- AllGenesList$GSE129002_CSB_SH.SY5Y
AllGenesList$GSE129002_CSB_SH.SY5Y <- NULL

AllGenesList$CSB_MSC <- AllGenesList$GSE124208_CSB_MSC
AllGenesList$GSE124208_CSB_MSC <- NULL

AllGenesList$CSB_iPSC <- AllGenesList$GSE124208_CSB_iPSC
AllGenesList$GSE124208_CSB_iPSC <- NULL

AllGenesList$CSB_CS1AN_3 <- AllGenesList$GSE122736_CSB_CS1AN
AllGenesList$GSE122736_CSB_CS1AN <- NULL

AllGenesList$CSB_CS1AN_2 <- AllGenesList$GSE58068_CSB_CS1AN
AllGenesList$GSE58068_CSB_CS1AN <- NULL

AllGenesList$CSB_CS1AN_1 <- AllGenesList$GSE57923_CSB_CS1AN
AllGenesList$GSE57923_CSB_CS1AN <- NULL


upset(fromList (AllGenesList), 
      nintersects = 13,
      nsets = 11,
      empty.intersections = NULL,
      keep.order = T, 
      sets = c("Csa_mouse_24m", "Csa_mouse_12m", "Csa_mouse_1m", "CS_Patients","CSA_SH.SY5Y", 
               "CSB_SH.SY5Y", "CSB_MSC", "CSB_iPSC", "CSB_CS1AN_3",
               "CSB_CS1AN_2", "CSB_CS1AN_1"),
      order.by = "degree", 
      decreasing = TRUE, 
      number.angles = 0, 
      text.scale = 2, 
      point.size = 3.3, 
      line.size = 1.6,
      mb.ratio = c(0.5, 0.5),
)

#width=896&height=804



png(filename = "Fig1A.png", width = 9.5, height = 9, units = "in", res = 600)
upset(fromList (AllGenesList), 
      nintersects = 13,
      nsets = 11,
      empty.intersections = NULL,
      keep.order = T, 
      sets = c("Csa_mouse_24m", "Csa_mouse_12m", "Csa_mouse_1m", "CS_Patients","CSA_SH.SY5Y", 
               "CSB_SH.SY5Y", "CSB_MSC", "CSB_iPSC", "CSB_CS1AN_3",
               "CSB_CS1AN_2", "CSB_CS1AN_1"),
      order.by = "degree", 
      decreasing = TRUE, 
      number.angles = 0, 
      text.scale = 2, 
      point.size = 3.3, 
      line.size = 1.6,
      mb.ratio = c(0.5, 0.5),
)
dev.off()


jpeg(filename = "Fig1A.jpeg", width = 9.5, height = 9, units = "in", res = 600)
upset(fromList (AllGenesList), 
      nintersects = 13,
      nsets = 11,
      empty.intersections = NULL,
      keep.order = T, 
      sets = c("Csa_mouse_24m", "Csa_mouse_12m", "Csa_mouse_1m", "CS_Patients","CSA_SH.SY5Y", 
               "CSB_SH.SY5Y", "CSB_MSC", "CSB_iPSC", "CSB_CS1AN_3",
               "CSB_CS1AN_2", "CSB_CS1AN_1"),
      order.by = "degree", 
      decreasing = TRUE, 
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)
dev.off()





upset(fromList (AllGenesList), 
      nintersects = 13,
      nsets = 11,
      empty.intersections = NULL,
      keep.order = T, 
      sets = c("CSA_mouse_24m", "CSA_mouse_12m", "CSA_mouse_1m", "GSE58067_CS_Patients","GSE129002_CSA_SH.SY5Y", 
               "GSE129002_CSB_SH.SY5Y", "GSE124208_CSB_MSC", "GSE124208_CSB_iPSC", "GSE122736_CSB_CS1AN",
               "GSE58068_CSB_CS1AN", "GSE57923_CSB_CS1AN"),
      order.by = "degree", 
      decreasing = TRUE, 
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)

upset(fromList (AllGenesList), 
      nintersects = 13,
      nsets = 11,
      empty.intersections = NULL,
      keep.order = T, 
      sets = c("CSA_mouse_24m", "CSA_mouse_12m", "CSA_mouse_1m", "GSE58067_CS_Patients","GSE129002_CSA_SH.SY5Y", 
               "GSE129002_CSB_SH.SY5Y", "GSE124208_CSB_MSC", "GSE124208_CSB_iPSC", "GSE122736_CSB_CS1AN",
               "GSE58068_CSB_CS1AN", "GSE57923_CSB_CS1AN"),
      order.by = "degree", 
      decreasing = TRUE, 
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)

# 952 900

png(filename = "plot.png", width = 9.5, height = 9, units = "in", res = 600)

dev.off()


AllGenesList <- lapply(AllGenesList, function(sublist) sort(sublist))

upset(fromList(AllGenesList), 
      nintersects = 13,
      nsets = 11,
      empty.intersections = NULL,
      order.by = "degree", 
      decreasing = TRUE, 
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)




grid.text("All genes intersection",x = 0.65, y=0.98, gp=gpar(fontsize=10))



############## Só humanos ###############


setwd ("~/Desktop/CS_Paper/Todos_datasets_originais") 

AllGenes <- read.csv("All_genes_protein_only_human.csv")

AllGenes[AllGenes == ""] <- NA



AllGenesList <- as.list(AllGenes)


remove_empty_values <- function(x) {
  if (is.list(x)) {
    x <- x[sapply(x, function(y) !is.null(y) && length(y) > 0)]
    x <- lapply(x, remove_empty_values)
  }
  return(x)
}

# Apply the function to remove empty values
AllGenesList <- remove_empty_values(AllGenesList)


# Define a function to replace empty strings with NA
replace_empty_with_na <- function(x) {
  ifelse(x == "", NA, x)
}



AllGenesList <- lapply(AllGenesList, function(sublist) sublist[!is.na(sublist)])



upset(fromList (AllGenesList), 
      nintersects = 23,
      nsets = 11,
      empty.intersections = NULL,
      order.by = "degree", 
      decreasing = T, 
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)
grid.text("All genes intersection",x = 0.65, y=0.98, gp=gpar(fontsize=10))



