library(ROCR)

pred <- function(A, S, test) {
  rows <- nrow(test)
  K <- ncol(A)
  res <- matrix(0, nrow=rows, ncol=2)
  
  for (i in 1:rows) {
    u <- test[i,1]
    v <- test[i,2]
    link <- test[i,3]
    
    score <- sum(S[u,]*A[v,])
    
    res[i,1] <- score
    res[i,2] <- link
  }
  
  res <- as.data.frame(res)
  names(res) <- c('Prediction','Class')

  return(res)
}

mod <- '/Users/ettore/Desktop/champ matrix/synth-nr-comparison/s1/synth-cwn_s1_9f.pred'
net <- '/Users/ettore/Desktop/champ matrix/synth-nr-comparison/s1/s1.communities.nn1000.nc2048.k9.network.2_hops_inv.pred'
out <- '/Users/ettore/Desktop/champ matrix/synth-nr-comparison/s1/roc_s1.pdf'
A_file <- "/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s1/s1.communities.nn1000.nc2048.k9.A"
S_file <- "/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s1/s1.communities.nn1000.nc2048.k9.S"
Test_file <- '/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s1/s1.communities.nn1000.nc2048.k9.network.2_hops'
i = 1

A = as.matrix(read.table(A_file, sep = "\t", header = FALSE))
S = as.matrix(read.table(S_file, sep = "\t", header = FALSE))
test = read.table(Test_file, sep = "\t", header = T)

result_mod <- read.table(mod, header = TRUE)
result_net <- read.table(net, header = TRUE)

result_mod <- read.csv(file = mod, header = TRUE, sep="\t")
result_net <- read.csv(file = net, sep=" ", header=F)
names(result_net) <- c('c1','c2','Prediction','Class')
result_orig <- pred(A, S, test)

pred_mod <- prediction(result_mod$Prediction,result_mod$Class)
perf_mod <- performance(pred_mod,"tpr","fpr")
auc_mod <- performance(pred_mod,"auc")

pred_net <- prediction(result_net$Prediction,result_net$Class)
perf_net <- performance(pred_net,"tpr","fpr")
auc_net <- performance(pred_net,"auc")
x = unlist(perf_net@x.values)
y = unlist(perf_net@y.values)
index = which(!is.na(y))

pred_orig <- prediction(result_orig$Prediction,result_orig$Class)
perf_orig <- performance(pred_orig,"tpr","fpr")
auc_orig <- performance(pred_orig,"auc")
x2 = unlist(perf_orig@x.values)
y2 = unlist(perf_orig@y.values)
index2 = which(!is.na(y2))

s = c(paste("S_",i," (AUC = ",sprintf("%.3f",auc_mod@y.values),")", sep = ""),
      paste("NetRate (AUC = ",sprintf("%.3f",auc_net@y.values),")", sep = ""),
      paste("Original (AUC = ",sprintf("%.3f",auc_orig@y.values),")", sep = "")
      )

pdf(file=out)
par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
plot(perf_mod,col=1,type='l',lty=1,lwd=3,xlab='FPR',ylab='TPR')
lines(x[index],y[index],col=2,lty=1,lwd=3)
lines(x2[index2],y2[index2],col=3,lty=1,lwd=3)
abline(0,1,col="grey",lty=3,lwd=2)
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="bottomright",s,col=1:3,lty=1,lwd=3,cex=0.8,bg="white")
dev.off()
