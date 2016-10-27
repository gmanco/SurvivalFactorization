library(igraph)
library(Matrix)
library(gplots)

#### This file generates cascades from a network

#Some utility functions
source("/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/generate_cascades.R")
source("/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/output_format.R")

#set the path and choose the file to build cascades upon
setwd("/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/random")
edgeListFile = "random.txt"
dat = read.table(file = edgeListFile, header = FALSE, sep = "\t") # choose an adjacency matrix from a .csv file

# generate the graph
# NOTICE: the node names are those denoted in the edgeListFile
g = graph.data.frame(dat, directed = TRUE)
numNodes = vcount(g)

#Set cascade parameters
numCascades = 1024
numFactors = 8
Tmax = 100
Tmin = 5

#OPTIONAL: plot the graph to verify
#igraph.options(vertex.size=0.2, vertex.label=NA, edge.arrow.size=0.1)
#plot(g, layout=layout.kamada.kawai)


# Generate the influence matrices to be used for the cascade generation
out <- generateInfluenceMatrices(g,numNodes,numFactors)

S = out[[1]]
A = out[[2]]

# Cascades are generated as a matrix numCascades x numNodes. Each cell contains the activation time. cells with value 0 are not active
cascades = generate_cascades(g,numNodes,numCascades,numFactors,Tmax,Tmin,S,A)

# Check the nodes who did not become active in any cascade. These nodes should be removed from the network
nodes_not_in_cascades = apply(t(cascades), 1, function(row) all(row ==0 ))

# a vector for node renaming
names = rep(0,numNodes)

# Now we associate a name to each node, by skipping those nodes to remove
k = 1
for (i in 1:numNodes){
  if (!nodes_not_in_cascades[i]){
    names[i] = k
    k = k + 1
  }
}

# the new graph now has new names
gm <- set.vertex.attribute(g, "name", value=names)

# and the nodes to remove are deleted from the graph
ids = 1:numNodes
gm <-  delete.vertices(gm,V(gm)[ids[nodes_not_in_cascades]])

# and next from the cascade matrix
newcascades = cascades[,!nodes_not_in_cascades]

## We are done: let's generate the output files
# Relational format
cascades_relational = paste(edgeListFile, ".remapped.cascades", sep = "")
generate_relational_format(newcascades,cascades_relational)

# Netrate format
cascade_netrate = paste(edgeListFile, ".remapped.cascades.netrate", sep = "")
network_netrate = paste(edgeListFile, ".remapped.network.netrate", sep = "")
generate_netrate_format(newcascades,gm,cascade_netrate,network_netrate)

# the graph of two hops (for testing the accuracy)
network_hops = paste(edgeListFile, ".remapped.2_hops", sep = "")
generate_two_hops(gm,network_hops,1)




