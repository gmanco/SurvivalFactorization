setwd("/Users/ettore/git/SurvivalFactorization/resources/datasets/ba")

orFile = "synth.BA.communities.nn50.nc64.k2.clusters"
mdFile = "models/ba_2f.clusters"
orData <- read.csv(file = orFile, sep="\t", header = F)
orData = cbind(1:dim(orData)[1],orData)
names(orData)=c("cascade", "cluster")
#df <- aggregate(cascade ~ cluster, data=orData, length)
#barplot(df, main="orig", xlab = "clusters", ylab = "freq")
#box()
#grid (NULL,NULL, lty = 6, col = "cornsilk2")


barplot(as.matrix(orData)[,2]%%2 + 1, main="orig", xlab = "cascade", ylab = "cluster")


n = 70
i=0
i=1

for (i in 0:n) {
  mdData <- read.csv(file = paste(mdFile, "_it_", i, sep = ""), sep="\t", header = F)
  mdData <- cbind(1:dim(mdData)[1],mdData)
  names(mdData)=c("cascade", "cluster")#, "c1", "c2")
  #df <- aggregate(cascade ~ cluster, data=mdData, length)
  #barplot(t(df), main=paste("m",i, sep = ""), xlab = "clusters", ylab = "freq")
  #box()
  #grid (NULL,NULL, lty = 6, col = "cornsilk2")
  par(mfrow= c(1,2))
  barplot(as.matrix(orData)[,2]%%2 + 1, main="orig", xlab = "cascade", ylab = "cluster")
  barplot(as.matrix(mdData)[,2], main=paste("m",i, sep = ""), xlab = "cascade", ylab = "cluster")
}

Af <- as.matrix(read.csv(file = "models/ba_2f.model.A", sep="\t", header = F))
Sf <- as.matrix(read.csv(file = "models/ba_2f.model.S", sep="\t", header = F))


E = read.table('synth.BA.communities.nn50.nc64.k2.network',sep = '\t',header = TRUE)

M = sparseMatrix(i = E[,1], j = E[,2], x=nrow(E), dims=c(50,50))

M_hat = Sf %*% t(Af)
min_tol = 1e-1

pr = c(0,0)
#pr[2] = sum(sum(M_hat>min_tol & M>min_tol))/sum(sum(M>min_tol));
#pr[1] = sum(sum(M>min_tol))/sum(sum(M_hat>min_tol));
pr[2] = sum(M_hat>min_tol & M>min_tol)/sum(M>min_tol)
pr[1] = sum(M_hat>min_tol & M>min_tol)/sum(M_hat>min_tol)
pr

sum(M > min_tol)
sum(M_hat > min_tol)
sum(M_hat > min_tol & M > min_tol)

which(M_hat>min_tol, arr.in=TRUE)

plot(sort(as.vector((M_hat - min(M_hat)) / (max(M_hat - min(M_hat))) + 1)), log = "y")
     
     
F = read.table('synth.BA.communities.nn50.nc64.k2.cascades',sep = '\t',header = TRUE)
C = sparseMatrix(i = F[,1], j = F[,2], x=F[,3], dims=c(50,64))
C[26,]
K = which(C[26,] > 0)
which(C[33,] > 0)
