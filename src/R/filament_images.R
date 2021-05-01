library(TDA)
library(caret)
library(tictoc)
library(doParallel)
library(sp)
library(raster)
library(e1071)
library(data.table)

setwd("~/Projects/cytoplasmic_streaming/data/original/actin_07102019/")

loadData <- function(pattern){
  files <- dir(pattern, full.names = T)
  return(lapply(files, read.table))
}

rawData <- lapply(c("825cp/", "1650cp/", "3300cp/"), loadData)
rawData <- do.call('c', rawData)
sapply(rawData, nrow)

labels <- as.factor(rep(c('825','1650','3300'), each=50))



rawData_multiLines <- lapply(seq(rawData), function(i){
  dat <- rawData[[i]]
  SpatialLines(list(Lines(lapply(split(dat[,1:2], dat[,3]), Line), ID=i)))
})


# r <- raster()
# extent(r) <- c(0,20,0,20)
# ncol(r) <- 200
# nrow(r) <- 200
# 
# rasterizeML <- function(ml, r){
#   ras <- raster::rasterize(ml, r, 1)
#   raster::values(ras)[is.na(raster::values(ras))] <- 0
#   return(ras)
# }
# 
# cl <- makeCluster(4)  
# registerDoParallel(cl) 
# 
# rawData_rasters <- foreach(ml=rawData_multiLines) %dopar% rasterizeML(ml, r)
# 
# save(rawData_multiLines, rawData_rasters, file='../../altered/filament_images/linesAndRasters.Rdata')
load('../../altered/filament_images/linesAndRasters.Rdata')
  
rawData_mats <- lapply(rawData_rasters, function(ras){
  as.matrix(coordinates(ras)[values(ras)==1,])
})
summary(sapply(rawData_mats, nrow))

set.seed(1)
rawData_mats_3s <- lapply(rawData_mats, function(mat){
  samp <- sample(seq(nrow(mat)), size = 3000, replace=F)
  list(s1=mat[samp[1:1000],],s2=mat[samp[1001:2000],],s3=mat[samp[2001:3000],])
})

computeRips <- function(data) TDA::ripsDiag(data, 1, 10)

# tic()
# rips_s1 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s1)
# rips_s2 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s2)
# rips_s3 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s3)
# toc()

#save(rips_s1, rips_s2, rips_s3, file='../../altered/filament_images/rips_3s_filtrations.Rdata')
load('../../altered/filament_images/rips_3s_filtrations.Rdata')

allRips <- c(rips_s1, rips_s2, rips_s3)

metadat <- data.frame(labels=rep(labels,3), sampleNum=rep(1:3, each=length(rips_s1)), fold=NA, networkNum = rep(1:150,3))

set.seed(7)
folds <- createFolds(metadat$labels, k=10)
for(i in seq(folds)){
  metadat$fold[as.integer(rownames(metadat))%in%folds[[i]]] <- i
}

alldata <- lapply(allRips, '[[', 'diagram')

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

summary_stats_forplot <- summary_stats #Save all columns for plotting later on

summary_stats <- summary_stats[,!names(summary_stats)%in%c('mean0.Birth', 'sd0.Birth')] # THese are constants that should probably just be removed...

summary_stats$label <- metadat$labels
summary_stats$sampleNum <- metadat$sampleNum
summary_stats$networkNum <- metadat$networkNum
summary_stats$fold <- metadat$fold


xValidReuslts <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats[summary_stats$fold!=fold,], as.factor(label) ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
  test <- summary_stats[summary_stats$fold==fold,]
  test$preds <- predict(svmMod, test)
  results <- c()
  for(net in unique(test$networkNum)){
    t <- table(test[test$networkNum==net,'preds'])
    results <- c(results,test$label[test$networkNum==net][1]==names(which(t==max(t))))
  }
  xValidReuslts <- c(xValidReuslts,sum(results)/length(results))
}
round(mean(xValidReuslts)*100,2)


#### Check skewed class example
set.seed(1)
results <- lapply(seq(.1,.9,.1), function(reducFact){
  results <- c()
  for(lev in levels(summary_stats$label)){
    networksToCull <- unique(summary_stats$networkNum[summary_stats$label==lev])
    outOfSample <- sample(networksToCull,length(networksToCull)*reducFact)
    # summary_stats_skewed <- summary_stats[!summary_stats$networkNum%in%outOfSample | summary_stats$sampleNum!=1,]
    summary_stats_skewed <- summary_stats[!summary_stats$networkNum%in%outOfSample,]
    summary_stats_skewed$fold <- NA
    
    folds <- createFolds(summary_stats_skewed$label, k=5)
    for(i in seq(folds)){
      summary_stats_skewed$fold[seq(1:nrow(summary_stats_skewed))%in%folds[[i]]] <- i
    }
    
    
    a <- round(mean(sapply(na.exclude(unique(summary_stats_skewed$fold)), function(fold){
      svmMod <- svm(data=summary_stats_skewed[!summary_stats_skewed$fold==fold,], label ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
      test <- summary_stats_skewed[summary_stats_skewed$fold==fold,]
      test <- test[which(test$sampleNum==1),]
      sum(predict(svmMod, test)==test$label)/nrow(test)*100
    })),2)
    results <- c(results, a)
  }
  return(results)
})


rdf <- data.frame(do.call(rbind,results))
names(rdf) <- levels(summary_stats$label)
rdf <- data.frame(accuracy=c(rdf$`1650`,rdf$`3300`,rdf$`825`), less=rep(seq(.1,.9,.1),3), crosslinkers=rep(c('1650','3300','825'),each=nrow(rdf)))

ggplot()+geom_line(data=rdf,aes(x=less,y=accuracy,group=crosslinkers,color=crosslinkers))+scale_x_continuous(breaks=seq(.1,.9,.1)) + 
  theme(panel.grid.minor.x = element_blank()) + scale_y_continuous(breaks=seq(0,100,1))


#####








### Was Distance

#Cut top 20
cutTop20 <- function(ripserPD, cut=.8){
  pd0 <- ripserPD$diagram[ripserPD$diagram[,'dimension']==0,]
  pd1 <- ripserPD$diagram[ripserPD$diagram[,'dimension']==1,]
  
  persistence0 <- pd0[,'Death']-pd0[,'Birth']
  persistence1 <- pd1[,'Death']-pd1[,'Birth']
  
  pd0_cut <- quantile(persistence0, cut)
  pd1_cut <- quantile(persistence1, cut)
  
  pd0 <- pd0[persistence0>=pd0_cut,]
  pd1 <- pd1[persistence1>=pd1_cut,]
  
  ripserPD$diagram <- rbind(pd0,pd1)
  
  return(ripserPD)
}
top20Rips <- lapply(allRips, cutTop20)


loopOver <- function(dat, i, dim){
  c(rep(0, i),
    sapply((i+1):450, function(j){
      TDA::wasserstein(dat[[i]]$diagram,dat[[j]]$diagram,p=2,dimension =dim)^(0.5)
    }))
}

cl <- makeCluster(60)  
registerDoParallel(cl) 


tic()
a <- foreach(i=seq(449)) %dopar% loopOver(top20Rips, i, 0)
b <- foreach(i=seq(449)) %dopar% loopOver(top20Rips, i, 1 )
toc()

save(a,b, file = '~/Projects/cytoplasmic_streaming/data/altered/wasserstein_distances_top20_3bs.Rdata')
load('~/Projects/cytoplasmic_streaming/data/altered/wasserstein_distances_top20_3bs.Rdata')


was0 <- do.call(rbind, a)
was1 <- do.call(rbind, b)
was0 <- rbind(was0, rep(0,length(a[[1]])))
was1 <- rbind(was1, rep(0,length(b[[1]])))

was0 <- t(was0)+was0
was1 <- t(was1)+was1


### Test wasserstein 3xsampled distance
netNums <- unique(metadat$networkNum)
was0_means <- matrix(0,nrow = 150, ncol=150)
for(i in netNums){
  rn <- as.numeric(row.names(metadat[metadat$networkNum==i,]))
  for(j in netNums){
    rn2 <- as.numeric(row.names(metadat[metadat$networkNum==j,]))
    was0_means[i,j] <- mean(was0[rn2,rn])
  }
}
diag(was0_means) <- 0

was1_means <- matrix(0,nrow = 150, ncol=150)
for(i in netNums){
  rn <- as.numeric(row.names(metadat[metadat$networkNum==i,]))
  for(j in netNums){
    rn2 <- as.numeric(row.names(metadat[metadat$networkNum==j,]))
    was1_means[i,j] <- mean(was1[rn2,rn])
  }
}
diag(was1_means) <- 0

ss<-50
actuallabels <- as.factor(rep(c(1,2,3),each=ss))

# Distance from one sample to one class
# Input: item - label of the sample, class - No. of the class, dismatrix - distance matrix
Dis_to_class <- function(item,class,dismatrix,classTrain){
  if(class == 1){
    average_dist=mean(dismatrix[item,classTrain[[1]]])
  }
  if(class == 2){
    average_dist=mean(dismatrix[item,classTrain[[2]]])
  }
  if(class == 3){
    average_dist=mean(dismatrix[item,classTrain[[3]]])
  }
  return(average_dist)
}

# classify the sample
# Input: item, w=weight of 0 feature
classify <- function(item,w,dist_0,dist_1,classTrain){
  Dis <- matrix(0,2,3)
  for(i in 1:3){
    Dis[1,i] <- w*Dis_to_class(item,i,dist_0,classTrain)
    Dis[2,i] <- (1-w)*Dis_to_class(item,i,dist_1,classTrain)
  }
  bsum <- colSums(Dis)
  label <- which(bsum == min(bsum))
  return(label)
}

# cross validation: 9 folds vs  1 fold
# returns vector of accuracy rates
kFoldClassify <- function(classify, dist_0, dist_1, actuallabels, ss, Nf){
  # create 10 folds
  set.seed(1)
  folds <- createFolds(actuallabels,k = Nf)
  
  accrate <- rep(NA, 21)
  for (j in 0:20) {
    dataIndex <- seq(actuallabels)
    testLabels <- rep(NA, length(actuallabels))
    Ts <- ss-(ss/Nf)
    for (i in 1:Nf) {
      Testing  <- dataIndex[folds[[i]]]
      Training <- dataIndex[-folds[[i]]]
      c1Training <- Training[Training<=ss]
      c2Training <- Training[Training>ss & Training<=(ss*2)]
      c3Training <- Training[Training>(ss*2)]
      testLabels[folds[[i]]] <- sapply(Testing,w=0.05*j,dist_0=dist_0, dist_1=dist_1,
                                       classTrain=list(c1Training,c2Training,c3Training),classify)
    }
    testLabels <- factor(sapply(testLabels, function(x) c('825','1650','3300')[x]), levels=levels(actuallabels))
    accrate[j+1]=sum(testLabels==actuallabels)/length(actuallabels)
  }
  return(accrate)
}

accrate <- kFoldClassify(classify, dist_0=was0_means, dist_1=was1_means, metadat$labels[1:150], ss=50, Nf=10)
max(accrate)
plot(seq(0,1, by=.05), accrate, type='o', ylab='Accuracy Rate', xlab='Dim-0 Weight')



###  FIgure for paper showing algo ###
p1 <- ggplot() + 
  geom_point(data=rawData[[1]],aes(V1,V2,color=as.factor(V3))) + 
  scale_color_discrete() +theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

# Make ggplot-friendly version of multiline
i=0
df = data.frame(V1=NULL,V2=NULL,V3=NULL)
for(coords in coordinates(rawData_multiLines[[1]])[[1]]){
  df <- rbind(df, cbind(coords,group=i,order=seq(nrow(coords))))
  i=i+1
}


p2 <- ggplot() + 
  geom_path(data=fortify(SpatialLinesDataFrame(rawData_multiLines[[1]],dat=data.frame(row.names = 0:99))),aes(x=long,y=lat, group=group, color=as.factor(group))) +
  scale_color_discrete() +theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())


# Make ggplot-friendly version of raster
rdr_spdf <- as(rawData_rasters[[1]], "SpatialPixelsDataFrame")
rdr_df <- as.data.frame(rdr_spdf)

p3 <- ggplot() +  
  geom_tile(data=rdr_df, aes(x=x, y=y, fill=as.factor(layer),alpha=as.factor(layer))) + scale_fill_manual(values = c('#ffffff','#000000')) +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

rawData_mats_3s_exe <- data.frame(do.call(rbind, rawData_mats_3s[[1]]))
rawData_mats_3s_exe$Sample <- as.factor(rep(1:3,each=1000))


p4 <- ggplot() +
  geom_point(data=rawData_mats_3s_exe, aes(x=x,y=y, color=Sample)) + scale_color_manual(values = c('red','blue','green')) +
  theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())


#Make version of PDs for ggplot
rips_exe_df <- data.frame(
  dimension=c(rips_s1[[1]]$diagram[,'dimension'],rips_s2[[1]]$diagram[,'dimension'],rips_s3[[1]]$diagram[,'dimension']),
  Birth=c(rips_s1[[1]]$diagram[,'Birth'],rips_s2[[1]]$diagram[,'Birth'],rips_s3[[1]]$diagram[,'Birth']),
  Death=c(rips_s1[[1]]$diagram[,'Death'],rips_s2[[1]]$diagram[,'Death'],rips_s3[[1]]$diagram[,'Death']),
  Sample=c(rep(1, nrow(rips_s1[[1]]$diagram)),rep(2, nrow(rips_s2[[1]]$diagram)),rep(3, nrow(rips_s3[[1]]$diagram))))


p5 <- ggplot() + 
  geom_point(data=rips_exe_df, aes(x=Birth,y=Death, ,shape=as.factor(dimension), color=as.factor(Sample))) + 
  scale_color_manual(values = c('red','blue','green')) +
  theme(legend.position="none")

all <- data.frame(summary_stats_forplot[c(1,151,301),c('mean.Birth',"mean.Death", "sd.Birth", "sd.Death")])
names(all) <- c('mb','md','sdb','sdd')
all$Sample <-1:3
zeroth <- data.frame(summary_stats_forplot[c(1,151,301),c('mean0.Birth',"mean0.Death", "sd0.Birth", "sd0.Death")])
names(zeroth) <- c('mb','md','sdb','sdd')
zeroth$Sample <-1:3
oneth <- data.frame(summary_stats_forplot[c(1,151,301),c('mean1.Birth',"mean1.Death", "sd1.Birth", "sd1.Death")])
names(oneth) <- c('mb','md','sdb','sdd')
oneth$Sample <-1:3

ssexe <- rbind(all,zeroth,oneth)
ssexe$features <- rep(c('all', '0th', '1th'), each=3)

p6 <- ggplot() + geom_point(data=ssexe, aes(x=mb,y=md,color=as.factor(Sample), shape=features, size=sdb)) + 
  scale_color_manual(values = c('red','blue','green')) + xlab('Birth') + ylab('Death') +
  theme(legend.position="none")


grid.arrange(p1 +  geom_text(aes(x=2,y=19,label='A'),cex=10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
             p2 +  geom_text(aes(x=2,y=19,label='B'),cex=10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
             p3 +  geom_text(aes(x=2,y=19,label='C'),cex=10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
             p4 +  geom_text(aes(x=2,y=19,label='D'),cex=10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
             p5 +  geom_text(aes(x=.2,y=9.6,label='E'),cex=10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
             p6 +  geom_text(aes(x=.05,y=.76,label='F'),cex=10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
             nrow=2)

#########################################################
