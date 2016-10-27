library(ROCR)

setwd("/Users/ettore/Desktop/Java projects/Survival Factorization/resources/datasets/weibo/preds")
pdf(file="AUC.pdf")
  
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
  
  s <- rep("",7)
  
  data <- read.csv(file = "weibo_dpu_2f.pred",sep="\t")

  pred <- prediction(data$Prediction,data$Class)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  f1 <- performance(pred,"f")
  s[1] = paste("Weibo 2 Topics (AUC = ",sprintf("%.3f",auc@y.values),")")
  plot(perf,col=1, lty=1,lwd=3)
  


data <- read.csv(file = "weibo_dpu_4f.pred",sep="\t")

pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
f2 <- performance(pred,"f")
s[2] = paste("Weibo 4 Topics (AUC = ",sprintf("%.3f",auc@y.values),")")
plot(perf,col=2, lty=1,lwd=3,add=TRUE)


data <- read.csv("weibo_dpu_8f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
f3 <- performance(pred,"f")
s[3] = paste("Weibo 8 Topics (AUC = ",sprintf("%.3f",auc@y.values),")")
plot(perf,col=3, lty=1,lwd=3,add=TRUE)


data <- read.csv("weibo_dpu_16f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
f4 <- performance(pred,"f")
s[4] = paste("Weibo 16 Topics (AUC = ",sprintf("%.3f",auc@y.values),")")
plot(perf,col=4, lty=1,lwd=3,add=TRUE)

data <- read.csv("weibo_dpu_32f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
f5 <- performance(pred,"f")
s[5] = paste("Weibo 32 Topics (AUC = ",sprintf("%.3f",auc@y.values),")")
plot(perf,col=5, lty=1,lwd=3,add=TRUE)

data <- read.csv("weibo_dpu_64f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
f6 <- performance(pred,"f")
s[6] = paste("Weibo 64 Topics (AUC = ",sprintf("%.3f",auc@y.values),")")
plot(perf,col=6, lty=1,lwd=3,add=TRUE)

data <- read.csv("weibo_dpu_128f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
f7 <- performance(pred,"f")
s[7] = paste("Weibo 128 Topics (AUC = ",sprintf("%.3f",auc@y.values),")")
plot(perf,col=7, lty=1,lwd=3,add=TRUE)


  
  abline(0,1,col="grey",lty=3,lwd=2)
  abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
  abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
  
  
  legend(x="bottomright",c(s[1],s[2],s[3],s[4],s[5],s[6],s[7]),col=c(1,2,3,4,5,6,7),lty=1,lwd=3,cex=0.5)


dev.off()

png(file="F1.png")

par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex=1.3,cex.axis=1.3,cex.lab=1.3)

plot(f6,col=6, lty=1,lwd=4,ylim = c(0,1), log="x")
plot(f2,col=2, lty=1,lwd=4,add=TRUE)
plot(f3,col=3, lty=1,lwd=4,add=TRUE)
plot(f3,col=4, lty=1,lwd=4,add=TRUE)
plot(f3,col=5, lty=1,lwd=4,add=TRUE)
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="topright",c(s[6]),col=c(6),lty=1,lwd=4,cex=1)

dev.off()

png(file="tmp.png")

data <- read.csv("weibo_dpu_64f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)

perf <- performance(pred,"f")
plot(perf, log="x")

dev.off()
