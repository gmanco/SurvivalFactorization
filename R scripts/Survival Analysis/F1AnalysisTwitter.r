library(ROCR)

setwd("/Users/ettore/git/SurvivalFactorization/resources/datasets/twitter/preds")

s <- rep("",8)

data <- read.csv(file = "twitter_2f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf1 <- performance(pred,"f")
# precrec1 <- performance(pred,"prec", "rec")
s[1] = paste("Twitter 2 Topics")


data <- read.csv(file = "twitter_4f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf2 <- performance(pred,"f")
s[2] = paste("Twitter 4 Topics")


data <- read.csv(file = "twitter_8f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf3 <- performance(pred,"f")
s[3] = paste("Twitter 8 Topics")


data <- read.csv(file = "twitter_16f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf4 <- performance(pred,"f")
s[4] = paste("Twitter 16 Topics")


data <- read.csv(file = "twitter_32f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf5 <- performance(pred,"f")
s[5] = paste("Twitter 32 Topics")


data <- read.csv(file = "twitter_64f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf6 <- performance(pred,"f")
s[6] = paste("Twitter 64 Topics")


data <- read.csv(file = "twitter_128f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf7 <- performance(pred,"f")
s[7] = paste("Twitter 128 Topics")

data <- read.csv(file = "netrate_twitter_clean.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf8 <- performance(pred,"f")
s[8] = paste("Netrate Twitter")

rm(data)
rm(pred)


png(file="TwitterFM.png")

par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
plot(perf1,col=1,lty=1,lwd=3,log="x",ylim=c(0,0.51))#,xlim=c(1e-45,1e-7))
plot(perf2,col=2,lty=1,lwd=3,log="x",add=TRUE)
plot(perf3,col=3,lty=1,lwd=3,log="x",add=TRUE)
plot(perf4,col=4,lty=1,lwd=3,log="x",add=TRUE)
plot(perf5,col=5,lty=1,lwd=3,log="x",add=TRUE)
plot(perf6,col=6,lty=1,lwd=3,log="x",add=TRUE)
plot(perf7,col=7,lty=1,lwd=3,log="x",add=TRUE)
plot(perf8,col=8,lty=1,lwd=3,log="x",add=TRUE)
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="bottomleft",c(s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8]),col=c(1,2,3,4,5,6,7,8),lty=1,lwd=4,cex=1)

dev.off()
