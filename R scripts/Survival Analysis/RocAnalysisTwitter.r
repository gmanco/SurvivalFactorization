library(ROCR)

setwd("/Users/ettore/Desktop/twitter 1 stats final/full")
png(file="Twitter_network_reconstruction_AUC.png")

par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex=1.3,cex.axis=1.3,cex.lab=1.3)

exps = 6
#s <- rep("",exps + 1)
s <- rep("",exps)

for(i in 1:exps) {
  print(i)
  topics = 2^i
  data <- read.csv(file = paste("twitter2_", topics, "f.pred", sep = ""),sep="\t",header=T)
  
  pred <- prediction(data$Prediction,data$Class)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  s[i] = paste(topics," Topics (AUC = ",sprintf("%.3f",auc@y.values),")", sep = "")

  if (i == 1)
    plot(perf,col=i,lty=1,lwd=3,main='Roc Curve (Twitter-Large)')
  else
    plot(perf,col=i,lty=1,lwd=3,add=TRUE)
}

#data <- read.csv(file = "netrate_twitter_network_2_hops.pred",sep=" ",header=F)
#names(data) <- c('c1','c2','Prediction','Class')
#pred <- prediction(1/data$Prediction,data$Class)
#perf <- performance(pred,"tpr","fpr")
#auc <- performance(pred,"auc")
#s[exps + 1] = paste("Netrate (AUC = ",sprintf("%.3f",auc@y.values),")", sep = "")
#plot(perf,col=exps + 1,lty=1,lwd=3,add=TRUE)

abline(0,1,col="grey",lty=3,lwd=2)
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)

#legend(x="bottomright",s,col=c(1:(exps+1)),lty=1,lwd=3,cex=0.8, bg="white")
legend(x="bottomright",s,col=c(1:exps),lty=1,lwd=3,cex=0.8, bg="white")

dev.off()
