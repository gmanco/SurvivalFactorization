library(ROCR)

setwd("/Users/ettore/git/SurvivalFactorization/resources/datasets/random/preds")

s <- rep("",7)

data <- read.csv(file = "random_2f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf1 <- performance(pred,"f")
precrec1 <- performance(pred,"prec", "rec")
s[1] = paste("Twitter 2 Topics")


data <- read.csv(file = "random_4f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf2 <- performance(pred,"f")
precrec2 <- performance(pred,"prec", "rec")
s[2] = paste("Twitter 4 Topics")


data <- read.csv(file = "random_8f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf3 <- performance(pred,"f")
precrec3 <- performance(pred,"prec", "rec")
s[3] = paste("Twitter 8 Topics")


data <- read.csv(file = "random_16f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf4 <- performance(pred,"f")
precrec4 <- performance(pred,"prec", "rec")
s[4] = paste("Twitter 16 Topics")


data <- read.csv(file = "random_32f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf5 <- performance(pred,"f")
precrec5 <- performance(pred,"prec", "rec")
s[5] = paste("Twitter 32 Topics")


data <- read.csv(file = "random_64f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf6 <- performance(pred,"f")
precrec6 <- performance(pred,"prec", "rec")
s[6] = paste("Twitter 64 Topics")


data <- read.csv(file = "random_128f.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf7 <- performance(pred,"f")
precrec7 <- performance(pred,"prec", "rec")
s[7] = paste("Twitter 128 Topics")

rm(data)
rm(pred)

png(file="randomPR.png")

par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)

plot(precrec1,col=1,lty=1,lwd=3,log="x")#,ylim=c(0.2,0.31),xlim=c(1e-150,1e-5))
plot(precrec2,col=2,lty=1,lwd=3,log="x",add=TRUE)
plot(precrec3,col=3,lty=1,lwd=3,log="x",add=TRUE)
plot(precrec4,col=4,lty=1,lwd=3,log="x",add=TRUE)
plot(precrec5,col=5,lty=1,lwd=3,log="x",add=TRUE)
plot(precrec6,col=6,lty=1,lwd=3,log="x",add=TRUE)
plot(precrec7,col=7,lty=1,lwd=3,log="x",add=TRUE)

plot(perf1,col=1,lty=1,lwd=3,log="x")#,ylim=c(0.2,0.31),xlim=c(1e-150,1e-5))
plot(perf2,col=2,lty=1,lwd=3,log="x",add=TRUE)
plot(perf3,col=3,lty=1,lwd=3,log="x",add=TRUE)
plot(perf4,col=4,lty=1,lwd=3,log="x",add=TRUE)
plot(perf5,col=5,lty=1,lwd=3,log="x",add=TRUE)
plot(perf6,col=6,lty=1,lwd=3,log="x",add=TRUE)
plot(perf7,col=7,lty=1,lwd=3,log="x",add=TRUE)

abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="topright",c(s[1],s[2],s[3],s[4],s[5],s[6],s[7]),col=c(1,2,3,4,5,6,7),lty=1,lwd=4,cex=1)

dev.off()
