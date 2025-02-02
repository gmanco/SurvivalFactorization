library(ROCR)

setwd("/Users/ettore/Desktop/Java projects/Survival Factorization/resources/datasets/weibo/preds")
png(file="WeiboDunfFM.png")

par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex=1.3,cex.axis=1.3,cex.lab=1.3)

s <- rep("",10)

data <- read.csv(file = "weibo_dunf_2f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[1] = paste("Weibo 2 Topics")
plot(perf,col=1, lty=1,lwd=3,log="x",ylim=c(0,0.16))

data <- read.csv(file = "weibo_dunf_4f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[2] = paste("Weibo 4 Topics")
plot(perf,col=2, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_8f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[3] = paste("Weibo 8 Topics")
plot(perf,col=3, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_16f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[4] = paste("Weibo 16 Topics")
plot(perf,col=4, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_32f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[5] = paste("Weibo 32 Topics")
plot(perf,col=5, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_64f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[6] = paste("Weibo 64 Topics")
plot(perf,col=6, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_128f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[7] = paste("Weibo 128 Topics")
plot(perf,col=7, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_256f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[8] = paste("Weibo 256 Topics")
plot(perf,col=8, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_512f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[9] = paste("Weibo 512 Topics")
plot(perf,col=9, lty=1,lwd=3,log="x",add=TRUE)

data <- read.csv(file = "weibo_dunf_1024f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[10] = paste("Weibo 1024 Topics")
plot(perf,col=10, lty=1,lwd=3,log="x",add=TRUE)

abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="topleft",c(s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10]),col=c(1,2,3,4,5,6,7,8,9,10),lty=1,lwd=4,cex=1)

dev.off()
