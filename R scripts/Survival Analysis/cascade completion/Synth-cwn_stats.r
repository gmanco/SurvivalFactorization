library(ROCR)

splits = c(50,60,70,80)
topics <- c(9,7,11,6)


for (split in splits){
  predFile <- '/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/s'
  outputFile <- paste('/Users/ettore/git/SurvivalFactorization/resources/datasets/synth-cwn/synth-cwn_s1_result_',split,"_", sep="")
  
  fmFile <- paste(outputFile, "FM.pdf", sep = "")
  rocFile <- paste(outputFile, "ROC.pdf", sep = "")
  prFile <- paste(outputFile, "PR.pdf", sep = "")
  
  main_title=paste("Split=",split,"%",sep="")
  
  i = 0
  
  pdf(file=fmFile)
  s = rep("", length(topics))
  
  for (k in topics) {
    i = i + 1
    
    predFile2 <- paste(predFile,i,"/synth-cwn_s",i,'_result_',split,"_", sep="")
    result <- read.table(paste(predFile2, k, 'f.pred', sep = ""), header = TRUE)
    result$predictedts[which(is.infinite(result$predictedts))] = 1e40
    raw_pred <- as.data.frame(cbind(1 / result$predictedts, as.numeric(result$truets > 0)))
    names(raw_pred) = c("Prediction", "Class")
    
    pred <- prediction(raw_pred$Prediction,raw_pred$Class)
    fm <- performance(pred,"f",x.measure = "cutoff")
    
    xlimit = 10^c(-6,-1)
    
    x = unlist(fm@x.values)
    y = unlist(fm@y.values)
    index = which(!is.na(y))
    
    if (i==1) {
      par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
      plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,log='x',xlab='Threshold',ylab='F1-Measure',xlim=xlimit,ylim=c(0,1),xaxt='n',main=main_title)
      axis(1,at=axTicks(1),labels=1/axTicks(1))
    } else
      lines(x[index],y[index],col=i,lty=1,lwd=3)
    
    bestFmInd <- which.max(fm@"y.values"[[1]])
    best = fm@"y.values"[[1]][bestFmInd]
    s[i] <- paste("S_", i, " (F = ",sprintf("%.3f",best),")",sep = "")
    zF = paste(zF, paste(
      k, sprintf("%.9f",fm@"y.values"[[1]][bestFmInd]),fm@"x.values"[[1]][bestFmInd], sep = "\t"
    ), sep = "\r\n")
  }
  
  abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
  abline(v=10^seq(-6,-1,1),col="grey",lwd=0.5)
  legend(x="bottomleft",s,col=1:length(topics),lty=1,lwd=3,cex=0.8,bg="white")
  dev.off()
  
  
  
  i = 0
  
  pdf(file=rocFile)
  
  for (k in topics) {
    i = i + 1
    
    predFile2 <- paste(predFile,i,"/synth-cwn_s",i,'_result_',split,"_", sep="")
    f <- paste(predFile2, k, 'f.pred', sep = "")
    result <- read.table(f, header = TRUE)
    result$predictedts[which(is.infinite(result$predictedts))] = 1e40
    raw_pred <- as.data.frame(cbind(1 / result$predictedts, as.numeric(result$truets > 0)))
    names(raw_pred) = c("Prediction", "Class")
    
    pred <- prediction(raw_pred$Prediction,raw_pred$Class)
    roc <- performance(pred,"tpr","fpr")
    auc <- performance(pred,"auc")
    
    s_auc[i] <- paste("S_",i," (AUC = ",sprintf("%.3f",auc@y.values),")",sep = "")
    
    x = unlist(roc@x.values)
    y = unlist(roc@y.values)
    index = which(!is.na(y))
    
    if (i==1) {
      par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
      plot(x[index],y[index],col=i,type='l',lty=1,lwd=3,xlab='FPR',ylab='TPR',main=main_title)
    } else {
      lines(x[index],y[index],col=i,lty=1,lwd=3)
    }
  }
  
  abline(0,1,col="grey",lty=3,lwd=2)
  abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
  abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
  legend(x="bottomright",s_auc,col=1:length(topics),lty=1,lwd=3,cex=0.8,bg="white")
  dev.off()
  
  
  i = 0
  s = rep("", length(topics))
  
  pdf(file=prFile)
  
  for (k in topics) {
    i = i + 1
    
    predFile2 <- paste(predFile,i,"/synth-cwn_s",i,'_result_',split,"_", sep="")
    result <- read.table(paste(predFile2, k, 'f.pred', sep = ""), header = TRUE)
    result$predictedts[which(is.infinite(result$predictedts))] = 1e40
    raw_pred <- as.data.frame(cbind(1 / result$predictedts, as.numeric(result$truets > 0)))
    names(raw_pred) = c("Prediction", "Class")
    
    pred <- prediction(raw_pred$Prediction,raw_pred$Class)
    pr <- performance(pred,"prec","rec")
    
    x = unlist(pr@x.values)
    y = unlist(pr@y.values)
    index = which(!is.na(y))
    s[i] <- paste("S_",i,sep = "")
    
    if (i==1) {
      par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
      plot(pr,col=i,type='l',lty=1,lwd=3,xlab='Recall',ylab='Precision',xlim=c(0,1),ylim=c(0,1),main=main_title)
    } else {
      lines(x[index],y[index],col=i,lty=1,lwd=3)
    }
  }
  
  abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
  abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
  legend(x="bottomleft",s,col=1:length(topics),lty=1,lwd=3,cex=0.8, bg="white")
  dev.off()
}
