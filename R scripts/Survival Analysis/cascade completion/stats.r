n <- 6
predFile <- '/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s1/synth-cwn_s1_result_50_'
outputFile <- '/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s1/synth-cwn_s1_result_50_'

fmFile <- paste(outputFile, "FM.pdf", sep = "")
rocFile <- paste(outputFile, "ROC.pdf", sep = "")
prFile <- paste(outputFile, "PR.pdf", sep = "")

s = rep("", n)
s_auc = rep("", n)

pdf(file=fmFile)
zF <- '#Topics\tF1-Measure\tCutoff'

for (i in 1:n) {
  k = 2^i
  result <- read.table(paste(predFile, k, 'f.pred', sep = ""), header = TRUE)
  raw_pred <- as.data.frame(cbind(1 / result$predictedts, as.numeric(result$truets > 0)))
  names(raw_pred) = c("Prediction", "Class")
  raw_pred$Prediction[which(is.infinite(raw_pred$Prediction))] = 1e40
  pred <- prediction(raw_pred$Prediction,raw_pred$Class)
  
  fm <- performance(pred,"f",x.measure = "cutoff")
  
  s[i] <- paste("#Topics", k)

  x = unlist(fm@x.values)
  y = unlist(fm@y.values)
  index = which(!is.na(y))
  if (i==1) {
    par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
    plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,log='x',xlab='Threshold',ylab='F1-Measure')#,ylim=c(0,0.71),xlim=c(1e-15,1e-6))
  } else
    lines(x[index],y[index],col=i,lty=1,lwd=3)

  bestFmInd <- which.max(fm@"y.values"[[1]])
  zF = paste(zF, paste(
    k, sprintf("%.9f",fm@"y.values"[[1]][bestFmInd]),fm@"x.values"[[1]][bestFmInd], sep = "\t"
  ), sep = "\r\n")
}

abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=10^seq(-27,-7,5),col="grey",lwd=0.5)
legend(x="bottomright",s,col=1:n,lty=1,lwd=3,cex=0.8,bg="white")
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
    plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,xlab='FPR',ylab='TPR')
  } else
    lines(x[index],y[index],col=i,lty=1,lwd=3)
}

abline(0,1,col="grey",lty=3,lwd=2)
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="bottomright",s_auc,col=1:n,lty=1,lwd=3,cex=0.8,bg="white")
dev.off()


pdf(file=prFile)

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
  if (i==1) {
    par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
    plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,xlab='Recall',ylab='Precision',xlim=c(0,1),ylim=c(0,1))
  } else
    lines(x[index],y[index],col=i,lty=1,lwd=3)
}

dev.set(dev.next())
abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
legend(x="bottomright",s,col=1:7,lty=1,lwd=3,cex=0.8, bg="white")
dev.off()


cat(zF)
