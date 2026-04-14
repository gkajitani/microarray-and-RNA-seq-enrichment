library(ggplot2)
library(reshape2)
library(dplyr)
library (gplots)
library(RColorBrewer)
library(edgeR)

# Set the working directory and read the data
setwd ("~/Desktop/CS_Paper")
degs <- read.csv("DEGs_todos_paper.csv")

# Melt the data for plotting
data1 <- melt(degs, id.vars = "Genes")

# Generate the heatmap with borders on the tiles
heatmap_plot <- data1 %>%
  ggplot(aes(x = variable, y = Genes, fill = value)) +
  geom_tile(color = "black", size = 0.2) +  # Adds the heatmap tiles with borders
  scale_fill_distiller(palette = "RdBu", direction = -1, na.value = "#bbbbbb") +  # Lighter gray for NA values
  labs(x = "Model", y = "Gene", fill = "Log2Fold Change") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),  # Rotate labels vertically and align left
        axis.title.x = element_blank(),  # Optionally remove the x-axis title
        axis.ticks.x = element_blank(),  # Optionally remove the x-axis ticks
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        legend.title = element_text(size = 8)) +  # Adjust the font size of the legend title
  scale_x_discrete(position = "top") +  # Move the x-axis labels to the top
  coord_fixed()  # Make tiles square

print(heatmap_plot)

heatmap_plot <- data1 %>%
  ggplot(aes(x = variable, y = Genes, fill = value)) +
  geom_tile(color = "black", size = 0.2) +  # Adds the heatmap tiles with borders
  scale_fill_gradient2(low = "#104cb5", mid = "white", high = "#cc0000", midpoint = 0, na.value = "#bbbbbb") +  # Custom gradient colors
  labs(x = "Model", y = "Gene", fill = "Log2Fold Change") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),  # Rotate labels vertically and align left
        axis.title.x = element_blank(),  # Optionally remove the x-axis title
        axis.ticks.x = element_blank(),  # Optionally remove the x-axis ticks
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        legend.title = element_text(size = 8)) +  # Adjust the font size of the legend title
  scale_x_discrete(position = "top") +  # Move the x-axis labels to the top
  coord_fixed()  # Make tiles square

print(heatmap_plot)


data1 <- data1 %>%
  mutate(value_category = cut(value,
                              breaks = c(-Inf, -2, -1, -0.001, 0.001, 1, 2, Inf),
                              labels = c("Very Low", "Low", "Slightly Low", "Zero", "Slightly High", "High", "Very High")))

# Generate the heatmap with borders on the tiles and custom color scale
heatmap_plot <- data1 %>%
  ggplot(aes(x = variable, y = Genes, fill = value_category)) +
  geom_tile(color = "black", size = 0.2) +  # Adds the heatmap tiles with borders
  scale_fill_manual(values = c("Very Low" = "#B2182B",  # Dark Red
                               "Low" = "#EF8A62",      # Light Red
                               "Slightly Low" = "#FDDBC7",  # Pinkish
                               "Zero" = "white",        # White for zero
                               "Slightly High" = "#D1E5F0",  # Light Blue
                               "High" = "#67A9CF",      # Medium Blue
                               "Very High" = "#2166AC"),  # Dark Blue
                    na.value = "#bbbbbb") +  # Lighter gray for NA values
  labs(x = "Model", y = "Gene", fill = "Log2Fold Change") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),  # Rotate labels vertically and align left
        axis.title.x = element_blank(),  # Optionally remove the x-axis title
        axis.ticks.x = element_blank(),  # Optionally remove the x-axis ticks
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        legend.title = element_text(size = 8)) +  # Adjust the font size of the legend title
  scale_x_discrete(position = "top") +  # Move the x-axis labels to the top
  coord_fixed()  # Make tiles square

# Print the plot
print(heatmap_plot)

data2 <- degs
row.names(data2) <- data2$Genes
data2$Genes <- NULL

data2[is.na(data2)] <- 0

m <- as.matrix(data2)

heatmap.2(m,col=brewer.pal(11,"RdBu"),scale="row", trace="none")
