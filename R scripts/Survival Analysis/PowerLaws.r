setwd("/Users/ettore/Desktop/Java projects/Survival Factorization/resources/datasets/weibo/")
data <- read.table(file = "weibo_dpu_activations.txt", sep = "\t", header = T)
data <- data[,2:3]

pdf(file="CascadeFrequency_weibo.pdf")

count <- function(x) { length(na.omit(x)) } 
counts.users <- data.frame(aggregate(data$ItemId,by=list(data$UserId),FUN="count"))
names(counts.users) = c('UserId','count')
sorted = counts.users[order(counts.users$count,decreasing = TRUE),]
recounts.users <- data.frame(aggregate(sorted$UserId,by=list(sorted$count),FUN="count"))
names(recounts.users) = c('Cascades','Frequeny')

plot(recounts.users,log="xy")

dev.off()

pdf(file="UserFrequency_weibo.pdf")

counts.cascades <- data.frame(aggregate(data$UserId,by=list(data$ItemId),FUN="count"))
names(counts.cascades) = c('ItemId','count')
sorted = counts.cascades[order(counts.cascades$count,decreasing = TRUE),]
recounts.cascades <- data.frame(aggregate(sorted$ItemId,by=list(sorted$count),FUN="count"))
names(recounts.cascades) = c('Users','Frequeny')
plot(recounts.cascades,log="xy")

dev.off()

setwd("/Users/ettore/Desktop/Java projects/Survival Factorization/resources/datasets/twitter/")
data <- read.table(file = "activations.txt", sep = "\t", header = T)
data <- data[,1:2]

pdf(file="CascadeFrequency_twitter.pdf")

count <- function(x) { length(na.omit(x)) } 
counts.users <- data.frame(aggregate(data$ItemId,by=list(data$UserId),FUN="count"))
names(counts.users) = c('UserId','count')
sorted = counts.users[order(counts.users$count,decreasing = TRUE),]
recounts.users <- data.frame(aggregate(sorted$UserId,by=list(sorted$count),FUN="count"))
names(recounts.users) = c('Cascades','Frequeny')

plot(recounts.users,log="xy")

dev.off()

pdf(file="UserFrequency_twitter.pdf")

counts.cascades <- data.frame(aggregate(data$UserId,by=list(data$ItemId),FUN="count"))
names(counts.cascades) = c('ItemId','count')
sorted = counts.cascades[order(counts.cascades$count,decreasing = TRUE),]
recounts.cascades <- data.frame(aggregate(sorted$ItemId,by=list(sorted$count),FUN="count"))
names(recounts.cascades) = c('Users','Frequeny')
plot(recounts.cascades,log="xy")

dev.off()

