library(ROCR)

optRoc.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = x^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

setwd("/Users/ettore/Desktop/synth/core-periphery/data")

zF = "nn\tnc\tk\tModel\tF1-Measure\tcutoff"

dataset = "core-periphery.communities"
nodes = c(813)
cascades = c(1024)
topics = c(2,4,8,16,32,64,128)

NRbool = 0

# dataFile <- "random.communities.nn842.nc2048.k48.cascades_48_128f.model.pred"
# dataFileNR <- "random.communities.nn842.nc2048.k48.network.2_hops.netrate.inv.out"

for (nn in nodes) {
  for (nc in cascades) {
    for (k in topics) {
      dataFile <- paste(dataset, ".nn", nn, ".nc", nc, ".k15.cascades_15_", k, "f.model.pred", sep = "")
      dataFileNR <- paste(dataset, ".nn", nn, ".nc", nc, ".k15.network.2_hops.netrate.inv.out", sep = "")

      fmFile <- paste(dataset, ".nn", nn, ".nc", nc, ".k", k, ".fm.pdf", sep = "")
      rocFile <- paste(dataset, ".nn", nn, ".nc", nc, ".k", k, ".roc.pdf", sep = "")
      prFile <- paste(dataset, ".nn", nn, ".nc", nc, ".k", k, ".pr.pdf", sep = "")

      data <- read.csv(file = dataFile,sep="\t")
      pred <- prediction(data$Prediction,data$Class)

      fm <- performance(pred,"f",x.measure="cutoff")
      pr <- performance(pred,"prec", "rec")
      roc <- performance(pred,"tpr","fpr")
      auc <- performance(pred,"auc")
      
      if (NRbool == 1) {
        data <- read.csv(file = dataFileNR,sep=" ",header=F)
        names(data) = c("src","dat", "Prediction", "Class")
        pred <- prediction(data$Prediction,data$Class)

        fmNR <- performance(pred,"f",x.measure="cutoff")
        prNR <- performance(pred,"prec", "rec")
        rocNR <- performance(pred,"tpr","fpr")
        aucNR <- performance(pred,"auc")
      }

      bestFmInd <- which.max(fm@"y.values"[[1]])
      zF = paste(zF, paste(nn, nc, k, "Model", fm@"y.values"[[1]][bestFmInd], fm@"x.values"[[1]][bestFmInd], sep = "\t"), sep = "\r\n")
  
      if (NRbool == 1) {
        bestFmInd <- which.max(fmNR@"y.values"[[1]])
        zF = paste(zF, paste(nn, nc, k, "NetRate", fmNR@"y.values"[[1]][bestFmInd], fmNR@"x.values"[[1]][bestFmInd], sep = "\t"), sep = "\r\n")
      }

      pdf(file=fmFile)

      par(mar=c(5,5,2,2),xaxs="i",yaxs="i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
      plot(fm,col=1,lty=1,lwd=3,log="x")
      if (NRbool == 1) { plot(fmNR,col=2,lty=1,lwd=3,log="x",add=TRUE) }
      abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
      abline(v=seq(0,1,0.1),col="grey",lwd=0.5)
      
      if (NRbool == 1) {
      legend(x="bottomleft",c("Model", "NetRate"),col=c(1,2),lty=1,lwd=4,cex=1)
      } else {
        legend(x="bottomleft",c("Model"),col=c(1),lty=1,lwd=4,cex=1)
      }

      dev.off()

      pdf(file=rocFile)

      m1 <- paste("Model (AUC = ",sprintf("%.3f",auc@y.values),")", sep = "")

      if (NRbool == 1) {
        m2 <- paste("NetRate (AUC = ",sprintf("%.3f",aucNR@y.values),")", sep = "")
      }

      par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
      plot(roc,col=1, lty=1,lwd=3)
      if (NRbool == 1) { plot(rocNR,col=2, lty=1,lwd=3,add=TRUE) }
      abline(0,1,col="grey",lty=3,lwd=2)
      abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
      abline(v=seq(0,1,0.1),col="grey",lwd=0.5)

      if (NRbool == 1) {
        legend(x="bottomright",c(m1,m2),col=c(1,2),lty=1,lwd=3,cex=0.8, bg="white")
      } else {
        legend(x="bottomright",c(m1),col=c(1),lty=1,lwd=3,cex=0.8, bg="white")
      }

      dev.off()

      pdf(file=prFile)

      par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex=1.3,cex.axis=1.3,cex.lab=1.3)
      plot(pr,col=1, lty=1,lwd=3)
      if (NRbool == 1) { plot(prNR,col=2, lty=1,lwd=3,add=TRUE) }
      abline(h=seq(0,1,0.1),col="grey",lwd=0.5)
      abline(v=seq(0,1,0.1),col="grey",lwd=0.5)

      if (NRbool == 1) {
        legend(x="bottomright",c("Model","NetRate"),col=c(1,2),lty=1,lwd=3,cex=0.8, bg="white")
      } else {
        legend(x="bottomright",c("Model"),col=c(1),lty=1,lwd=3,cex=0.8, bg="white")
      }
      
      dev.off()
    }
  }
}

cat(zF)
