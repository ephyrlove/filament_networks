library(tiff)
library(raster)
library(spatstat)
library(doParallel)
library(tictoc)
library(caret)
library(e1071)
library(vegan)

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









image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/mutant/', full.names = T)
image_names <- c(image_names,dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/xik mutant 28 additional pics/', full.names = T))
image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/xik mutant 28 additional pics/', full.names = T)

image_names <- image_names[grepl('tif',image_names)]


# Read in images and scale based on measurements provided

mutants <- lapply(image_names, function(fn){
  if(!grepl('#',fn)){
    imnum <- as.numeric(strsplit(strsplit(fn,'.tif')[[1]],'-')[[1]][2])
  }else{
    imnum <- as.numeric(strsplit(strsplit(fn,'#')[[1]][2],'_')[[1]][1])
  }
  im <- readTIFF(fn)
  im_ras <- raster(im)
  
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  return(im_ras)
})

all_image_names <- image_names

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/wildType/', full.names = T)
# image_names <- c(image_names,dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/WT 30 additional pics/', full.names = T))
# image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/WT 30 additional pics/', full.names = T)
image_names <- image_names[grepl('tif',image_names)]

all_image_names <- c(image_names,all_image_names)

wildTypes <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'#')[[1]][2],'_')[[1]][1])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  return(im_ras)
})


# Find the smallest extents in x&y and crop in a bit\

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

cells <- append(wildTypes_threshold,mutants_threshold)

rawData_mats <- lapply(cells, function(ras){
  mat <- data.frame(as.matrix(coordinates(ras)))
  mat$value <- values(ras)
  mat$value <- (mat$value - min(mat$value))/(max(mat$value)-min(mat$value)) # added min-max normalization
  return(mat)
})

set.seed(1)

resamples<-3
totalSampled <- resamples*1000

rawData_mats_3s <- mclapply(rawData_mats, mc.cores=20,function(mat){
  samp <- sample(seq(nrow(mat)), size = totalSampled, replace=F, prob = mat$value)
  sampSeq <- seq(0,totalSampled,by=1000)
  l <- list(length=resamples)
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
folds <- createFolds(paste(metadat$labels,metadat$sampleNum), k=5)
for(i in seq(folds)){
  metadat$fold[row.names(metadat)%in%folds[[i]]] <- i
}

alldata <- lapply(rips, '[[', 'diagram')

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


# write.csv(summary_stats,'~/Projects/cytoplasmic_streaming/data/altered/vectorizedPDsForAN.csv',row.names = F)





summary_stats$pc1 <- princomp(summary_stats[c("mean.Birth","mean.Death","sd.Birth","sd.Death","mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,1]
summary_stats$pc2 <- princomp(summary_stats[c("mean.Birth","mean.Death","sd.Birth","sd.Death","mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,2]




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



wssplot <- function(data, nc=15, seed=1234)
{
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc)
  {
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)
  }
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

wssplot(summary_stats_1perCell[,c('pc1','pc2')])

kmod <- kmeans(summary_stats_1perCell[,c('pc1','pc2')],centers=2)

t <- table(kmod$cluster,summary_stats_1perCell$label)

clusterAssigments <- names(apply(t, MARGIN = 2,function(r)which(r==max(r))))


wrongC1 <- kmod$cluster==1 & summary_stats_1perCell$label!=clusterAssigments[1]
wrongC2 <- kmod$cluster==2 & summary_stats_1perCell$label!=clusterAssigments[2]


allCells <- append(wildTypes,mutants)

c1 <- as.numeric(which(kmod$cluster==1))
c2 <- as.numeric(which(kmod$cluster==2))

# png('~/Projects/cytoplasmic_streaming/fig/cluster_xamples.png', width = 7,height=11,units = 'in', res=150)
par(mfrow=c(50,2),mar=c(0.1,0.1,0.1,0.1))
for(i in seq(50)){
  image(allCells[[c1[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(wrongC1[c1[i]]){
    mp <- as.numeric(bbox(allCells[[c1[i]]])[,2]/2)
    points(y=mp[2],x=mp[1],pch=4,cex=10,col='red')
  }
  image(allCells[[c2[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(wrongC2[c2[i]]){
    mp <- as.numeric(bbox(allCells[[c2[i]]])[,2]/2)
    points(y=mp[2],x=mp[1],pch=4,cex=10,col='red')
  }
}




