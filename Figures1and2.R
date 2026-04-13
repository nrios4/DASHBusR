# Figures1and2.R

# Clear environment
rm(list = ls())

# Set working directory to the root directory of the archive, e.g.,
# setwd("C:/Users/riosn/OneDrive/Desktop/DASHBusRCodeR1")

# load libraries
library(igraph)
library(ggplot2)

bus_graph = readRDS("bus_graph.rds")
bus_data_X_Y = readRDS("bus_data_X_Y_format.rds")

X_edges = bus_data_X_Y[,-(32:34)]
day_edges = X_edges[275,]

# Set a default color and width for all edges so they are visible
E(bus_graph)$color = "grey"
E(bus_graph)$width = 1

highlighted_edge_inds = which(day_edges > 0)
for(i in 1:length(highlighted_edge_inds)){
  E(bus_graph)[highlighted_edge_inds[i]]$color = "red"
  E(bus_graph)[highlighted_edge_inds[i]]$width = 5
}


set.seed(1234)
plot(bus_graph, vertex.label.cex = 0.8,vertex.size = 3,vertex.label.dist = 1,vertex.label.degree = pi/4,
     layout = layout_with_kk)

X = as.matrix(X_edges)
barplot_df = data.frame(
  edge = colnames(X),
  frequency = colSums(X)
)

ggplot(barplot_df, aes(x = edge, y = frequency)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
