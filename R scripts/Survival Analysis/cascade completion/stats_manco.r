library(ROCR)

n <- 6

splits = c(50,60,70,80)

for (split in splits){
predFile <- paste('/Users/manco/Google Drive/SurvivalFactorization/twitter 1/time_prediction/TwitterSmall_pred/twitter_result_',split,'_',sep="")
outputFile <- paste('/Users/manco/Google Drive/SurvivalFactorization/twitter 1/time_prediction/TwitterSmall_pred/twitter_small_res_',split,'_',sep="")

fmFile <- paste(outputFile, "FM.pdf", sep = "")
rocFile <- paste(outputFile, "ROC.pdf", sep = "")
prFile <- paste(outputFile, "PR.pdf", sep = "")

main_title=paste("Split=",split,"%",sep="")


s = rep("", n)
s_auc = rep("", n)

pdf(file=fmFile)
zF <- '#Topics\tF1-Measure\tCutoff'
par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)

for (i in 1:n) {
  k = 2^i
  result <- read.table(paste(predFile, k, 'f.pred', sep = ""), header = TRUE)
  result$predictedts[which(is.infinite(result$predictedts))] = 1e40
    raw_pred <- as.data.frame(cbind(1 / result$predictedts, as.numeric(result$truets > 0)))
  names(raw_pred) = c("Prediction", "Class")
#  raw_pred$Prediction[which(is.infinite(raw_pred$Prediction))] = 1e40
  pred <- prediction(raw_pred$Prediction,raw_pred$Class)
  
  fm <- performance(pred,"f",x.measure = "cutoff")
  

  xlimit = 10^c(-12,-4)
  
  x = unlist(fm@x.values)
  y = unlist(fm@y.values)
  index = which(!is.na(y))
  
  
  if (i==1) {
    plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,log='x',xlab='Threshold',ylab='F1-Measure',xlim=xlimit,ylim=c(0,0.71),xaxt='n',main=main_title)
#    plot(fm,col=i,type='l',lty=1,lwd=3,ylim=c(0,1),xaxt = 'n')
    axis(1,at=axTicks(1),labels=1/axTicks(1))
           #formatC(axTicks(2),digits=0,format="f"))
  } else
    lines(x[index],y[index],col=i,lty=1,lwd=3)
#    plot(fm,col=i,type='l',lty=1,lwd=3,add=TRUE)

  bestFmInd <- which.max(fm@"y.values"[[1]])
  best = fm@"y.values"[[1]][bestFmInd]
  s[i] <- paste(k," Topics (F = ",sprintf("%.3f",best),")",sep = "")
  zF = paste(zF, paste(
    k, sprintf("%.9f",fm@"y.values"[[1]][bestFmInd]),fm@"x.values"[[1]][bestFmInd], sep = "\t"
  ), sep = "\r\n")
}

abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=10^seq(-12,-4,1),col="grey",lwd=0.5)
legend(x="bottomleft",s,col=1:n,lty=1,lwd=3,cex=0.8,bg="white")
dev.off()


pdf(file=rocFile)

for (i in 1:n) {
  k = 2^i
  result <- read.table(paste(predFile, k, 'f.pred', sep = ""), header = TRUE)
  raw_pred <- as.data.frame(cbind(1 / result$predictedts, as.numeric(result$truets > 0)))
  names(raw_pred) = c("Prediction", "Class")
  pred <- prediction(raw_pred$Prediction,raw_pred$Class)

  roc <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")

  s_auc[i] <- paste(k," Topics (AUC = ",sprintf("%.3f",auc@y.values),")",sep = "")

  x = unlist(roc@x.values)
  y = unlist(roc@y.values)
  index = which(!is.na(y))
  if (i==1) {
    par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
    plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,xlab='FPR',ylab='TPR',main=main_title)
  } else
    lines(x[index],y[index],col=i,lty=1,lwd=3)
}

abline(0,1,col="grey",lty=3,lwd=2)
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="bottomright",s_auc,col=1:n,lty=1,lwd=3,cex=0.8,bg="white")
dev.off()


pdf(file=prFile)
par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)

s = rep("", n)


for (i in 1:n) {
  k = 2^i
  result <- read.table(paste(predFile, k, 'f.pred', sep = ""), header = TRUE)
  raw_pred <- as.data.frame(cbind(1 / result$predictedts, as.numeric(result$truets > 0)))
  names(raw_pred) = c("Prediction", "Class")
  pred <- prediction(raw_pred$Prediction,raw_pred$Class)

  pr <- performance(pred,"prec","rec")

  x = unlist(pr@x.values)
  y = unlist(pr@y.values)
  index = which(!is.na(y))
  s[i] <- paste(k," Topics",sep = "")
  
  if (i==1) {
#    plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,xlab='Recall',ylab='Precision',xlim=c(0,1),ylim=c(0,1),main=main_title)
    plot(pr,col=i,type='l',lty=1,lwd=3,xlab='Recall',ylab='Precision',xlim=c(0,1),ylim=c(0,1),main=main_title)
  } else
    plot(pr,col=i,type='l',lty=1,lwd=3,xlab='Recall',ylab='Precision',xlim=c(0,1),ylim=c(0,1),main=main_title,add=TRUE)
#    lines(x[index],y[index],col=i,lty=1,lwd=3)
}

dev.set(dev.next())
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="bottomright",s,col=1:7,lty=1,lwd=3,cex=0.8, bg="white")
dev.off()


cat(zF)
}