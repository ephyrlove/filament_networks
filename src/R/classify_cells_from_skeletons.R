library(tiff)
library(raster)
library(spatstat)
library(doParallel)
library(tictoc)
library(caret)
library(e1071)
library(vegan)
library(EBImage)


scale_list_mutant = list()
scale_list_mutant[[1]] = c(0.044,	0.044)
scale_list_mutant[[2]] = c(0.044,	0.044)
scale_list_mutant[[3]] = c(0.043,	0.043)
scale_list_mutant[[4]] = c(0.043,	0.043)
scale_list_mutant[[5]] = c(0.038,	0.038)
scale_list_mutant[[6]] = c(0.038,	0.038)
scale_list_mutant[[7]] = c(0.038,	0.038)
scale_list_mutant[[8]] = c(0.044,	0.044)
scale_list_mutant[[9]] = c(0.044,	0.044)
scale_list_mutant[[10]] = c(0.044,	0.044)
scale_list_mutant[[11]] = c(0.044,	0.044)
scale_list_mutant[[12]] = c(0.044,	0.044)
scale_list_mutant[[13]] = c(0.044,	0.044)
scale_list_mutant[[14]] = c(0.044,	0.044)
scale_list_mutant[[15]] = c(0.052,	0.052)
scale_list_mutant[[16]] = c(0.052,	0.052)
scale_list_mutant[[17]] = c(0.044,	0.044)
scale_list_mutant[[18]] = c(0.044,	0.044)
scale_list_mutant[[19]] = c(0.044,	0.044)
scale_list_mutant[[20]] = c(0.044,	0.044)
scale_list_mutant[[21]] = c(0.044,	0.044)
scale_list_mutant[[22]] = c(0.044,	0.044)
scale_list_mutant[[23]] = c(0.049,	0.049)
scale_list_mutant[[24]] = c(0.043,	0.043)
scale_list_mutant[[25]] = c(0.044,	0.044)
scale_list_mutant[[26]] = c(0.044,	0.044)
scale_list_mutant[[27]] = c(0.044,	0.044)
scale_list_mutant[[28]] = c(0.044,	0.044)
scale_list_mutant[[29]] = c(0.052,	0.052)
scale_list_mutant[[30]] = c(0.043,	0.043)
scale_list_mutant[[31]] = c(0.043,	0.043)
scale_list_mutant[[32]] = c(0.049,	0.049)
scale_list_mutant[[33]] = c(0.043,	0.043)
scale_list_mutant[[34]] = c(0.042,	0.042)
scale_list_mutant[[35]] = c(0.044,	0.044)
scale_list_mutant[[36]] = c(0.051,	0.051)
scale_list_mutant[[37]] = c(0.043,	0.043)
scale_list_mutant[[38]] = c(0.044,	0.044)
scale_list_mutant[[39]] = c(0.043,	0.043)
scale_list_mutant[[40]] = c(0.047,	0.047)
scale_list_mutant[[41]] = c(0.042,	0.042)
scale_list_mutant[[42]] = c(0.042,	0.042)
scale_list_mutant[[43]] = c(0.044,	0.044)
scale_list_mutant[[44]] = c(0.042,	0.042)
scale_list_mutant[[45]] = c(0.043,	0.043)
scale_list_mutant[[46]] = c(0.042,	0.042)
scale_list_mutant[[47]] = c(0.049,	0.049)
scale_list_mutant[[48]] = c(0.047,	0.047)
scale_list_mutant[[49]] = c(0.044,	0.044)
scale_list_mutant[[50]] = c(0.051,	0.051)


scale_list_wt = list()
scale_list_wt[[1]] = c(0.043,	0.043)
scale_list_wt[[2]] = c(0.043,	0.043)
scale_list_wt[[3]] = c(0.043,	0.043)
scale_list_wt[[4]] = c(0.044,	0.044)
scale_list_wt[[5]] = c(0.043,	0.043)
scale_list_wt[[6]] = c(0.037,	0.037)
scale_list_wt[[7]] = c(0.05,	0.05)
scale_list_wt[[8]] = c(0.043,	0.043)
scale_list_wt[[9]] = c(0.043,	0.043)
scale_list_wt[[10]] = c(0.043,	0.043)
scale_list_wt[[11]] = c(0.044,	0.044)
scale_list_wt[[12]] = c(0.044,	0.044)
scale_list_wt[[13]] = c(0.044,	0.044)
scale_list_wt[[14]] = c(0.044,	0.044)
scale_list_wt[[15]] = c(0.044,	0.044)
scale_list_wt[[16]] = c(0.044,	0.044)
scale_list_wt[[17]] = c(0.044,	0.044)
scale_list_wt[[18]] = c(0.044,	0.044)
scale_list_wt[[19]] = c(0.044,	0.044)
scale_list_wt[[20]] = c(0.044,	0.044)
scale_list_wt[[21]] = c(0.044,	0.044)
scale_list_wt[[22]] = c(0.044,	0.044)
scale_list_wt[[23]] = c(0.044,	0.044)
scale_list_wt[[24]] = c(0.047,	0.047)
scale_list_wt[[25]] = c(0.048,	0.048)
scale_list_wt[[26]] = c(0.049,	0.049)
scale_list_wt[[27]] = c(0.047,	0.047)
scale_list_wt[[28]] = c(0.047,	0.047)
scale_list_wt[[29]] = c(0.05,	0.05)
scale_list_wt[[30]] = c(0.048,	0.048)
scale_list_wt[[31]] = c(0.048,	0.048)
scale_list_wt[[32]] = c(0.043,	0.043)
scale_list_wt[[33]] = c(0.049,	0.049)
scale_list_wt[[34]] = c(0.05,	0.05)
scale_list_wt[[35]] = c(0.047,	0.047)
scale_list_wt[[36]] = c(0.05,	0.05)
scale_list_wt[[37]] = c(0.051,	0.051)
scale_list_wt[[38]] = c(0.044,	0.044)
scale_list_wt[[39]] = c(0.05,	0.05)
scale_list_wt[[40]] = c(0.043,	0.043)
scale_list_wt[[41]] = c(0.042,	0.042)
scale_list_wt[[42]] = c(0.047,	0.047)
scale_list_wt[[43]] = c(0.04,	0.04)
scale_list_wt[[44]] = c(0.04,	0.04)
scale_list_wt[[45]] = c(0.051,	0.051)
scale_list_wt[[46]] = c(0.051,	0.051)
scale_list_wt[[47]] = c(0.044,	0.044)
scale_list_wt[[48]] = c(0.043,	0.043)
scale_list_wt[[49]] = c(0.043,	0.043)
scale_list_wt[[50]] = c(0.043,	0.043)



image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy_filtered/xik-50/', full.names = T)

mutants <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'xik')[[1]][3],'-')[[1]][1])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  values(im_ras) <- (values(im_ras)-min(values(im_ras)))/(max(values(im_ras))-min(values(im_ras)))
  return(im_ras)
})


all_image_names <- image_names

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy_filtered/WT-50/', full.names = T)

wildTypes <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'WT')[[1]][3],'-')[[1]][1])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  values(im_ras) <- (values(im_ras)-min(values(im_ras)))/(max(values(im_ras))-min(values(im_ras)))
  return(im_ras)
})



xshrink <- .1
yshrink <- .05
threshold <- 1-.33

smallestX <- min(c(sapply(wildTypes, function(im) extent(im)[2]), sapply(mutants, function(im) extent(im)[2]))) * (1-xshrink)
smallestY <- min(c(sapply(wildTypes, function(im) extent(im)[4]), sapply(mutants, function(im) extent(im)[4]))) * (1-yshrink)


wildTypes_threshold <- lapply(wildTypes, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  newExtent <- extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
  im <- raster::crop(im,newExtent)
  extent(im) <- newExtent
  return(im)
})

mutants_threshold <- lapply(mutants, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  newExtent <- extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
  im <- raster::crop(im,newExtent)
  extent(im) <- newExtent
  return(im)
})



# TRY WITHOUT CROPPING





cells <- append(wildTypes_threshold,mutants_threshold)

rawData_mats <- lapply(cells, function(ras){
  mat <- data.frame(as.matrix(coordinates(ras)))
  mat$value <- values(ras)
  return(mat)
})


pSS<-floor(min(sapply(rawData_mats, function(mat) sum(mat[,3]>0)))/3)

set.seed(1)

resamples<-3
totalSampled <- resamples*pSS

rawData_mats_3s <- mclapply(rawData_mats, mc.cores=20, function(mat){
  samp <- sample(seq(nrow(mat)), size = totalSampled, replace=F, prob = mat$value)
  sampSeq <- seq(0,totalSampled,by=pSS)
  l <- list()
  for(i in seq(length(sampSeq)-1)){
    l[[i]]<-mat[samp[sampSeq[i]:sampSeq[i+1]],c('x','y')]
  }
  return(l)
})

computeRips <- function(data) {
  lapply(data, function(d){
    TDA::ripsDiag(d, 1, 10)
  })
}

cl <- makeCluster(20)
registerDoParallel(cl)

tic()
rips <- foreach(x=rawData_mats_3s) %dopar% computeRips(x)
toc()

rips <- unlist(rips,recursive = F)

metadat <- data.frame(labels=c(rep('WT', length(wildTypes_threshold)*resamples),rep('M', length(mutants_threshold)*resamples)),
                      sampleNum=rep(1:resamples, length(cells)),
                      fold=NA,
                      networkNum=rep(1:length(cells),each=resamples),
                      fileName=rep(sapply(strsplit(all_image_names, '//'),'[',2),each=resamples))

set.seed(2)
folds <- createFolds(metadat$labels, k=5)
for(i in seq(folds)){
  metadat$fold[row.names(metadat)%in%folds[[i]]] <- i
}

alldata <- lapply(rips, function(x) x[['diagram']])

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
summary_stats$fold <- metadat$fold
summary_stats$networkNum <- metadat$networkNum
summary_stats$fileName <- metadat$fileName
summary_stats$sampleNum <- metadat$sampleNum





summary_stats$pc1 <- princomp(summary_stats[c("mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,1]
summary_stats$pc2 <- princomp(summary_stats[c("mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,2]




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
summary(xValidReuslts)








summary_stats_1perCell <- summary_stats[summary_stats$sampleNum==1,]
set.seed(1)
folds <- createFolds(summary_stats_1perCell$label, k=5)
for(i in seq(folds)){
  summary_stats_1perCell$fold[row.names(summary_stats_1perCell)%in%folds[[i]]] <- i
}

xValidReuslts <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats_1perCell[summary_stats_1perCell$fold!=fold,], as.factor(label) ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
  test <- summary_stats_1perCell[summary_stats_1perCell$fold==fold,]
  test$preds <- predict(svmMod, test)
  xValidReuslts <- c(xValidReuslts,sum(diag(table(pred=test$preds,label=test$label)))/length(test$preds))
}

mean(xValidReuslts)



wt_params <- read.csv('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/WT parameters.csv')
mutant_params <- read.csv('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/xik parameters.csv')
wt_params$label <- 'WT'
mutant_params$label <- 'M'
params <- rbind(wt_params,mutant_params)

explanatoryCols <- c("Total.Filament.Length", "Filament.Density", "Occupancy", 
                     "X..Branches", "Branch.Freq.", "Branch.Density", "Mean.Angle", 
                     "Parallelness", "Mean.Intensity..Filaments.", "StdDev", "Coefficient.of.Variation",
                     "Skewness", "Threshold", "pixels..50", "pixels.51.100", "pixels.101.150", "pixels.151.200", 
                     "pixels.201.255", "X75th..ile", "X90th..ile", "Black.Pixels", "Saturated.Pixels",
                     "Black.Pixels.Fraction", "Saturated.Pixels.Fraction", "Maximal.Distance", "Median.Distance", 
                     "Ratio.lower.higher.Median", "Ratio.Upper.Lower.Median", "Ratio.Outer.Inner.Median",
                     "Ratio.lower.higher.max", "Ratio.Upper.Lower.max", "Ratio.Outer.Inner.max")


sumExplanCols <- c("mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")

summary_stats_andParams <- data.frame(params[explanatoryCols],summary_stats_1perCell[c(sumExplanCols,'label','batch')])


fullPCAmod <- prcomp(summary_stats_andParams[c(sumExplanCols,explanatoryCols)],center = T,scale. = T)

plot(fullPCAmod)




####### Make cluster plots ###########

kmod <- kmeans(x=fullPCAmod$x[,1:3],centers = 2)

t <- table(kmod$cluster,summary_stats_andParams$label)


c1 <- as.numeric(which(kmod$cluster==1))
c2 <- as.numeric(which(kmod$cluster==2))

png('~/Projects/cytoplasmic_streaming/fig/allVariable_k2_cluster1.png', width = 8,height=length(c1)*2,units = 'in', res=150)
par(mfrow=c(length(c1),1),mar=c(0.2,0.1,0.2,0.1),bg='black')
for(i in seq(length(c1))){
  image(cells[[c1[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(summary_stats_andParams$label[[c1[i]]]=='M'){
    if(summary_stats_andParams$batch[[c1[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'red',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'red',cex=30,which='plot')
    }
  }else{
    if(summary_stats_andParams$batch[[c1[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'green',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'green',cex=30,which='plot')
    }
  }
}
dev.off()


png('~/Projects/cytoplasmic_streaming/fig/allVariable_k2_cluster2.png', width = 8,height=length(c2)*2,units = 'in', res=150)
par(mfrow=c(length(c2),1),mar=c(0.2,0.1,0.2,0.1),bg='black')
for(i in seq(length(c2))){
  image(cells[[c2[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(summary_stats_andParams$label[[c2[i]]]=='M'){
    if(summary_stats_andParams$batch[[c2[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'red',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'red',cex=30,which='plot')
    }
  }else{
    if(summary_stats_andParams$batch[[c2[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'green',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'green',cex=30,which='plot')
    }
  }
}
dev.off()




kmod <- kmeans(x=fullPCAmod$x[,1:3],centers = 3)

t <- table(kmod$cluster,summary_stats_andParams$label)


c1 <- as.numeric(which(kmod$cluster==1))
c2 <- as.numeric(which(kmod$cluster==2))
c3 <- as.numeric(which(kmod$cluster==3))

png('~/Projects/cytoplasmic_streaming/fig/allVariable_k3_cluster1.png', width = 8,height=length(c1)*2,units = 'in', res=150)
par(mfrow=c(length(c1),1),mar=c(0.2,0.1,0.2,0.1),bg='black')
for(i in seq(length(c1))){
  image(cells[[c1[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(summary_stats_andParams$label[[c1[i]]]=='M'){
    if(summary_stats_andParams$batch[[c1[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'red',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'red',cex=30,which='plot')
    }
  }else{
    if(summary_stats_andParams$batch[[c1[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'green',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'green',cex=30,which='plot')
    }
  }
}
dev.off()


png('~/Projects/cytoplasmic_streaming/fig/allVariable_k3_cluster2.png', width = 8,height=length(c2)*2,units = 'in', res=150)
par(mfrow=c(length(c2),1),mar=c(0.2,0.1,0.2,0.1),bg='black')
for(i in seq(length(c2))){
  image(cells[[c2[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(summary_stats_andParams$label[[c2[i]]]=='M'){
    if(summary_stats_andParams$batch[[c2[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'red',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'red',cex=30,which='plot')
    }
  }else{
    if(summary_stats_andParams$batch[[c2[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'green',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'green',cex=30,which='plot')
    }
  }
}
dev.off()


png('~/Projects/cytoplasmic_streaming/fig/allVariable_k3_cluster3.png', width = 8,height=length(c3)*2,units = 'in', res=150)
par(mfrow=c(length(c3),1),mar=c(0.2,0.1,0.2,0.1),bg='black')
for(i in seq(length(c3))){
  image(cells[[c3[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(summary_stats_andParams$label[[c3[i]]]=='M'){
    if(summary_stats_andParams$batch[[c3[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'red',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'red',cex=30,which='plot')
    }
  }else{
    if(summary_stats_andParams$batch[[c3[i]]]=='Batch 1'){
      box(lty = 'solid', col = 'green',cex=30,which='plot')
    }else{
      box(lty = 'dashed', col = 'green',cex=30,which='plot')
    }
  }
}
dev.off()




