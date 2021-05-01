library(TDA)
library(caret)
library(tictoc)
library(doParallel)
library(neuralnet)
library(e1071)
library(class)

setwd("~/Projects/cytoplasmic_streaming/data/original/actin_13092019/")

load('../../altered/alldata_13092019.Rdata')

alldata <- lapply(alldata, '[[', 'diagram')

set.seed(1)
folds <- createFolds(summary_stats$label, k=10)

summary_stats <- data.frame(do.call(rbind, lapply(alldata, function(pd){
  pd0 <-pd[pd[,1]==0,2:3]
  pd1 <-pd[pd[,1]==1,2:3]
  
  m <- apply(pd[,2:3], MARGIN = 2, mean)
  sd <- apply(pd[,2:3], MARGIN = 2, sd)
  
  m0 <- apply(pd0, MARGIN = 2, mean)
  sd0 <- apply(pd0, MARGIN = 2, sd)
  
  m1 <- apply(pd1, MARGIN = 2, mean)
  sd1 <- apply(pd1, MARGIN = 2, sd)
  
  c(mean=m,sd=sd, mean0=m0,sd0=sd, mean1=m1,sd1=sd1)
})))


summary_stats_sum <- data.frame(do.call(rbind, lapply(alldata, function(pd){
  c(sum=nrow(pd[(pd[,1]==0 & pd[,3]>mean(pd[,3])),]))
})))

summary_stats$label <- rep(1:3, each=50)

summary_stats <- summary_stats[,!names(summary_stats)%in%c('mean0.Birth', 'sd0.Birth')] # THese are constants that should probably just be removed...

head(summary_stats)



mean(sapply(folds, function(fold){
  svmMod <- svm(data=summary_stats[-fold,], label ~ .)
  sum(round(predict(svmMod, summary_stats[fold,]))==summary_stats$label[fold])/length(fold)
}))


mean(sapply(folds, function(fold){
  nnMod <- neuralnet::neuralnet(data=summary_stats[-fold,], label ~ .)
  sum(round(predict(nnMod, summary_stats[fold,]))==summary_stats$label[fold])/length(fold)
}))

mean(sapply(folds, function(fold){
  knnMod <- knn(train=summary_stats[-fold,-11],test =summary_stats[fold,-11], cl = summary_stats$label[-fold], k=3)
  sum(knnMod==summary_stats$label[fold])/length(fold)
}))


mean(sapply(folds, function(fold){
  
  glmod1 <- glm(data=summary_stats[-fold,], as.numeric(label==1) ~ mean.Birth + mean.Death + sd.Birth +
                sd.Death, family = binomial())
  glmod2 <- glm(data=summary_stats[-fold,], as.numeric(label==2) ~ mean.Birth + mean.Death + sd.Birth +
                  sd.Death  , family = binomial())
  glmod3 <- glm(data=summary_stats[-fold,], as.numeric(label==3) ~ mean.Birth + mean.Death + sd.Birth +
                  sd.Death, family = binomial())
  
  preds <- cbind(
    predict(glmod1, summary_stats[fold,], type='response'),
    predict(glmod2, summary_stats[fold,], type='response'),
    predict(glmod3, summary_stats[fold,], type='response'))
  
  sum(apply(preds, MARGIN=1, function(r) which(r==max(r)))==summary_stats$label[fold])/length(fold)
  }))


glmod <- glm(data=summary_stats, as.numeric(label==2) ~ mean.Birth + mean.Death + sd.Birth +
                sd.Death, family = binomial())
summary(glmod)




pca <- princomp(summary_stats[,1:10])

comps <- data.frame(pca$scores[,1:5])
comps$label <- summary_stats$label

mean(sapply(folds, function(fold){
  
  glmod1 <- glm(data=comps[-fold,], label==1 ~ ., family = binomial())
  glmod2 <- glm(data=comps[-fold,], label==2 ~ ., family = binomial())
  glmod3 <- glm(data=comps[-fold,], label==3 ~ ., family = binomial())
  
  preds <- cbind(
    predict(glmod1, comps[fold,], type='response'),
    predict(glmod2, comps[fold,], type='response'),
    predict(glmod3, comps[fold,], type='response'))
  
  sum(apply(preds, MARGIN=1, function(r) which(r==max(r)))==summary_stats$label[fold])/length(fold)
}))



lapply(seq(1:5), function(i){
  set.seed(i)
  folds <- createFolds(summary_stats$label, k=10)
  mean(sapply(folds, function(fold){
    nnMod <- neuralnet::neuralnet(data=comps[-fold,], label ~ .)
    sum(round(predict(nnMod, comps[fold,]))==comps$label[fold])/length(fold)
  }))
})



