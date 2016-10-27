library(igraph)

edges = read.table(file="/Users/ettore/git/SurvivalFactorization/resources/datasets/twitter/network",header = TRUE,sep = "\t")
users <- 1:max(edges)

g <- graph.empty()
g <- add.vertices(g,nv=length(users),attr=list(name=users))

# we need to turn edgeList into a vector (and using names instead of indexes)
edgeListVec <- as.vector(t(as.matrix(edges)))
g <- add.edges(g,edgeListVec)
d.in.g = degree(g,mode='in')
d.out.g = degree(g,mode='out')

A = as.matrix(read.table(file = "/Users/ettore/Desktop/champ matrix/small/twitter_16f.model.A",sep = "\t"))
S = as.matrix(read.table(file = "/Users/ettore/Desktop/champ matrix/small/twitter_16f.model.S",sep = "\t"))

maxA <- cbind(1:nrow(A), apply(A,1,max))
maxS <- cbind(1:nrow(S), apply(S,1,max))

plot(as.vector(d.in.g), maxA[1:(nrow(A)-1),2], log='xy')
plot(maxS[1:(nrow(A)-1), 2], log(as.vector(d.out.g)+1))

cor(maxA[1:(nrow(A)-1),2], as.vector(d.in.g))
cor(maxS[1:(nrow(A)-1), 2], as.vector(d.out.g))
