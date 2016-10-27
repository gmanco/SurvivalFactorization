library(ROCR)

optRoc.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = x^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

setwd("/Users/ettore/Desktop/Java projects/Survival Factorization/resources/datasets/weibo/preds")

zRoc = "#topics\tsensitivity\tspecificity\tcutoff"
zF = "#topics\tF1-Measure\tcutoff"

i = 2

for(i in 1:10) {
  topics = 2^i
  data <- read.csv(file = paste("weibo_dunf_", topics, "f.pred", sep = ""),sep="\t")
  pred <- prediction(data$Prediction,data$Class)
  
  perf <- performance(pred,"tpr","fpr")
  cut <- optRoc.cut(perf, pred)
  zRoc = paste(zRoc, paste(topics, cut[1,1], cut[2,1], cut[3,1], sep = "\t"), sep = "\r\n")
  
  perf <- performance(pred,"f",x.measure="cutoff")
  bestAccInd <- which.max(perf@"y.values"[[1]])
  zF = paste(zF, paste(topics, perf@"y.values"[[1]][bestAccInd], perf@"x.values"[[1]][bestAccInd], sep = "\t"), sep = "\r\n")
}

cat(zRoc, "\n")
cat(zF, "\n")