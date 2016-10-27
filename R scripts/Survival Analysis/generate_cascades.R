require('igraph')
require('Matrix')

# The algorithm for generating cascades works like follows: 
# 1.    For each node n and each factor k we generate S[n,k] and A[n,k]
# 2.    We randomly choose the initial nodes for each cascade
# 3.    for each cascade to generate, we proceed like follows: 
# 3.1   randomly pick a node n among those active nodes and remove it from the list of active
# 3.2   try to activate all his neighbors
# 3.3   for each neighbor m
# 3.3.1 generate a delay D using S[m,k]*A[n,k]. Add the delay to the activation time of n to obtain the activation time of m
# 3.3.2 if the activation time is greater than Tmax, then we skip it
# 3.3.3 otherwise we add it to the list of active nodes
# 4     the algorithm terminates when the list of active nodes is empty
generate_cascades <- function(g,numNodes,numCascades,numFactors,Tmax,Tmin,S, A){
  in.d.g = degree(g, mode = "in")
  in.d.g = in.d.g / sum(in.d.g)
  
  out.d.g = degree(g, mode = "out")
  out.d.g = out.d.g / sum(out.d.g)
  
  #FIXME: users with few outgoing edges are more susceptibles, users with more incoming edges are more authoritative
  
  cascades = Matrix(0,
                    nrow = numCascades,
                    ncol = numNodes,
                    sparse = TRUE)
  
#  initial_nodes = sample(numNodes,numCascades, replace = TRUE,prob = in.d.g)
  totnodes = 1:numNodes
  
  
  for (c in 1:numCascades) {
    triggers = rep(0, numNodes)
    K = sample(numFactors, 1)
    
#    currnode = initial_nodes[c]
    
    currnode = sample(numNodes,1,prob = A[,K]/sum(A[,K]))
    
    
    timestamp = runif(1, 0, Tmin)
    cascades[c, currnode] = timestamp
    canExpand = TRUE
    while (canExpand) {
      cond = cascades[c,] > 0 & triggers == 0
      currnodes = totnodes[cond]
      if (length(currnodes) > 0) {
        if (length(currnodes) > 1) {
          p = A[currnodes, K] / sum(A[currnodes, K])
          currnode = sample(currnodes, size = 1, prob = p)
        } else {
          currnode = currnodes[1]
        }
        triggers[currnode] = 1
        nb = as.numeric(neighbors(g, currnode, mode = "in"))
        for (nextnode in nb) {
          if (cascades[c, nextnode] == 0) {
            rate = S[nextnode, K] * A[currnode, K]
            
            timestamp = cascades[c, currnode] + rexp(1, rate)
            if (timestamp < Tmax) {
              cascades[c, nextnode] <- timestamp
            }
          }
        }
      } else {
        canExpand = FALSE
      }
    }
  }
  return(cascades)
}

generateInfluenceMatrices <- function(g,numNodes,numFactors){
  # Generate the matrix with default (low) values
  S = matrix(runif(numNodes * numFactors, min = 0.001, max = 0.01), ncol = numFactors)
  
  A = matrix(runif(numNodes * numFactors, min = 0.001, max = 0.01), ncol = numFactors)
  
  for (k in 1:numFactors){
    #For each factor we choose a number of nodes whose influence is boosted
    # The nodes are split into batches 
    numNodesToChoose = ceiling(numNodes/(2*numFactors))
    nodesToChoose = sample(numNodes,numNodesToChoose)
    
    # boost the influence degree for those nodes
    A[nodesToChoose,k] = runif(nodesToChoose, min = 0.1, max = 1)
    
    # Collect all incoming nodes
    neighborsOfNodesToChoose = c() 
    for (i in 1:length(nodesToChoose)) {
      neighborsOfNodesToChoose = c(neighborsOfNodesToChoose,as.numeric(neighbors(g,nodesToChoose[i],mode="in")))
    }
    neighborsOfNodesToChoose = unique(neighborsOfNodesToChoose)
  
    # and boost their susceptibility as well
    if (length(neighborsOfNodesToChoose) > 0){
      S[neighborsOfNodesToChoose,k] = runif(neighborsOfNodesToChoose, min = 0.1, max = 1)
    }
    
  }
  
  output<-list(S,A)
  return(output)
}



