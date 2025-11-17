library (tidyverse)
library("wesanderson")

setwd ("~/Desktop/CS_Paper/Figure_3/csv")

#############
# Figure 3A #
#############

pal <- wes_palette("Zissou1", 10, type = "continuous")

pal <- rev(pal)

# "PAPER_fig3a_all.csv" contained a .csv with Database (GO:BP/CC/MF, KEGG, etc), Description (Ontology/Pathway), GeneRatio, Count, p.adjust from enrichment analysis of intersection
# of differentially expressed genes present in several datasets; The same goes for the other plots in Figure 3.
data <- read.csv ("PAPER_fig3a_all.csv")

# Helper function to add line breaks after 3 words
add_line_breaks <- function(text, n = 3) {
  words <- str_split(text, " ")[[1]]
  paste(str_wrap(paste(words, collapse = " "), width = n * 10), collapse = "\n")
}

# Apply the function to the Description column
data$Description <- sapply(data$Description, add_line_breaks)

# Reverse the order of the Description factor levels
data$Description <- factor(data$Description, levels = rev(data$Description))




#darjeeling
#custom_colors <- c("#FF0000", "#F98400", "#F2AD00")

#zissou
#custom_colors <- c("#F21A00", "#E1AF00", "#EBCC2A")

# Create the bar plot
ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "Gene ontologies - All datasets",
       x = "Number of associated genes",
       y = "")



fig3a <- ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "Gene ontologies - All datasets",
       x = "Number of associated genes",
       y = "")


ggsave("fig3A.png", plot = fig3a, width = 9.48, height = 5.51, dpi = 600, limitsize = FALSE)



#############
# Figure 3B #
#############

data <- read.csv ("PAPER_fig3b_all.csv")

# Helper function to add line breaks after 3 words
add_line_breaks <- function(text, n = 3) {
  words <- str_split(text, " ")[[1]]
  paste(str_wrap(paste(words, collapse = " "), width = n * 10), collapse = "\n")
}

# Apply the function to the Description column
data$Description <- sapply(data$Description, add_line_breaks)

# Reverse the order of the Description factor levels
data$Description <- factor(data$Description, levels = rev(data$Description))

# Create the bar plot
ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "KEGG and Reactome pathways - All datasets",
       x = "Number of associated genes",
       y = "")



fig3b <- ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "KEGG and Reactome pathways - All datasets",
       x = "Number of associated genes",
       y = "")


ggsave("fig3B.png", plot = fig3b, width = 9.48, height = 5.51, dpi = 600, limitsize = FALSE)



#############
# Figure 3C #
#############

data <- read.csv ("PAPER_fig3c_human.csv")

# Helper function to add line breaks after 3 words
add_line_breaks <- function(text, n = 3) {
  words <- str_split(text, " ")[[1]]
  paste(str_wrap(paste(words, collapse = " "), width = n * 10), collapse = "\n")
}

# Apply the function to the Description column
data$Description <- sapply(data$Description, add_line_breaks)

# Reverse the order of the Description factor levels
data$Description <- factor(data$Description, levels = rev(data$Description))


# Create the bar plot
ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "Ontologies and pathways - Human datasets",
       x = "Number of associated genes",
       y = "")

fig3c <- ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "Ontologies and pathways - Human datasets",
       x = "Number of associated genes",
       y = "")


ggsave("fig3C.png", plot = fig3c, width = 9.48, height = 5.51, dpi = 600, limitsize = FALSE)


#############
# Figure 3D #
#############


data <- read.csv ("PAPER_fig3d_mouse.csv")

# Helper function to add line breaks after 3 words
add_line_breaks <- function(text, n = 3) {
  words <- str_split(text, " ")[[1]]
  paste(str_wrap(paste(words, collapse = " "), width = n * 10), collapse = "\n")
}

# Apply the function to the Description column
data$Description <- sapply(data$Description, add_line_breaks)

# Reverse the order of the Description factor levels
data$Description <- factor(data$Description, levels = rev(data$Description))

# Create the bar plot
ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "Ontologies and pathways - Mouse datasets",
       x = "Number of associated genes",
       y = "")

fig3d <- ggplot(data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colours = pal)+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "lightgray", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # Minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.ticks.y = element_line(color = "black", size = 0.5)  # Axis ticks
  ) +
  labs(title = "Ontologies and pathways - Mouse datasets",
       x = "Number of associated genes",
       y = "")

ggsave("fig3D.png", plot = fig3d, width = 9.48, height = 5.51, dpi = 600, limitsize = FALSE)

