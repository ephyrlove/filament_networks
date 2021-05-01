library(e1071)
library(caret)
library(ggplot2)
library(REdaS)
options(stringsAsFactors = F)

wt_dat <- read.csv('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/WT-morphology-210203.csv')
mutant_dat <- read.csv('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/mutant-morphology-210203.csv')
mutant2_dat <- read.csv('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/xi2k-morphology-210210.csv')

wt_dat$class <- 'wt'
mutant_dat$class <- 'xik'
mutant2_dat$class <- 'xi2k'

# Drop columns not in common
mutant2_dat <- mutant2_dat[names(mutant2_dat)%in%names(mutant_dat)]

colLables <- strsplit(readLines('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/mutant-morphology-210203.csv',n = 1),',')[[1]]
colLables[13] <- "Ratio Left/Right\nMedian Dist"
colLables[14] <- "Ratio Outer/Inner\nMedian Dist"
colLables[4] <- "Mean Angle (rad)"

data <- data.frame(rbind(wt_dat,mutant_dat,mutant2_dat))
data$class <- as.factor(data$class)

folds <- createFolds(data$class, k=5)
data$fold <- NA

for(i in seq(folds)){
  data$fold[row.names(data)%in%folds[[i]]] <- i
}


summary(sapply(seq(length(folds)),function(f){
  mod <- svm(x=data[data$fold!=f,3:14],y=data[data$fold!=f,15])
  test <- data[data$fold==f,]
  test$preds <- predict(mod, test[,3:14])
  sum(test$preds==test$class)/nrow(test)
}))


data$Mean.Angle<- deg2rad(data$Mean.Angle)
data_split <- split(data,data$class)

plotDat <- data.frame(do.call(rbind, lapply(data_split, function(s) {
  df <- data.frame(do.call(rbind,lapply(s[3:14], function(col){
    data.frame(mean=mean(col),sd=sd(col))
  })))
  df$metric <- row.names(df)
  df$class <- s$class[1]
  return(df)
})))

png(filename = '~/Projects/cytoplasmic_streaming/fig/traditionalMetricsBarPlot.png',units = 'in',width = 10,height=5,res = 300)
ggplot(data=plotDat) + 
  geom_bar(aes(x=metric,y=mean,fill=class),stat='identity',position='dodge') + 
  geom_errorbar(aes(x=metric,ymin=mean-sd, ymax=mean+sd,fill=class), position='dodge') + 
  scale_x_discrete(breaks=names(data[3:14]), labels=colLables[3:14])+ theme(axis.text.x = element_text(angle = 25,vjust = .75)) + ylab('mean +/- SD')
  # ggtitle('Traditional metrics of actin organization by class')
dev.off()



######### PCA?
pcaDat <- prcomp(data[3:14],scale=T,center = T)
plot(pcaDat)
biplot(pcaDat)


data$PC1 <- pcaDat$x[,1]
data$PC2 <- pcaDat$x[,2]

ggplot(data=data) + geom_point(aes(x=PC1,y=PC2,color=class),size=2)


ggplot(data=data) + geom_density(aes(PC1,fill=class),alpha=.5)
ggplot(data=data) + geom_density(aes(PC2,fill=class),alpha=.5)



ggplot(data=data) + geom_density2d(aes(PC1,PC2,color=class))




