setwd("~/SUSPOLL/Code")

#install.packages('igraph')
library(igraph)

data = read.csv('DependenciesTable.csv',row.names=1)
#data[which(data=='NA')]=0
adjacencymatrix = as.matrix(data)

# build the graph object
network <- graph_from_adjacency_matrix(adjacencymatrix,mode='undirected')

# plot it
#plot(network)
plot(network, layout=layout.circle)
