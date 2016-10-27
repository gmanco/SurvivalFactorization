library(ROCR)

setwd("/Users/ettore/git/SurvivalFactorization/resources/datasets/twitter/preds")

s <- rep("",1)

data <- read.csv(file = "netrate_twitter_clean.pred",sep="\t")
pred <- prediction(data$Prediction,data$Class)
perf <- performance(pred,"f")
s[1] = paste("Netrate Twitter")

rm(data)
rm(pred)

pdf(file="NetrateTwitterFM.pdf")

par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
plot(perf,col=1,lty=1,lwd=3,log="x",ylim=c(0,0.51))#,xlim=c(1e-45,1e-7))
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x=1e-30,y=0.451,c(s[1]),col=c(1),lty=1,lwd=4,cex=1)

dev.off()
