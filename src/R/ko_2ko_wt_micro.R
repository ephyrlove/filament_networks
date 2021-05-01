library(tiff)
library(raster)
library(spatstat)
library(doParallel)
library(tictoc)
library(caret)
library(e1071)

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

scale_list_dubMutant = list()
scale_list_dubMutant[[1]] = c(0.04426,	0.04426)
scale_list_dubMutant[[2]] = c(0.04426,	0.04426)
scale_list_dubMutant[[3]] = c(0.04426,	0.04426)
scale_list_dubMutant[[4]] = c(0.04466,	0.04466)
scale_list_dubMutant[[5]] = c(0.04466,	0.04466)
scale_list_dubMutant[[6]] = c(0.04426,	0.04426)
scale_list_dubMutant[[7]] = c(0.04362,	0.04362)
scale_list_dubMutant[[8]] = c(0.04362,	0.04362)
scale_list_dubMutant[[9]] = c(0.04309,	0.04309)
scale_list_dubMutant[[10]] = c(0.05913,	0.05913)
scale_list_dubMutant[[11]] = c(0.04414,	0.04414)
scale_list_dubMutant[[12]] = c(0.04375,	0.04375)
scale_list_dubMutant[[13]] = c(0.04247,	0.04247)
scale_list_dubMutant[[14]] = c(0.04375,	0.04375)
scale_list_dubMutant[[15]] = c(0.04375,	0.04375)
scale_list_dubMutant[[16]] = c(0.04375,	0.04375)
scale_list_dubMutant[[17]] = c(0.04499,	0.04499)
scale_list_dubMutant[[18]] = c(0.04499,	0.04499)
scale_list_dubMutant[[19]] = c(0.04363,	0.04363)
scale_list_dubMutant[[20]] = c(0.04363,	0.04363)
scale_list_dubMutant[[21]] = c(0.04353,	0.04353)
scale_list_dubMutant[[22]] = c(0.04353,	0.04353)
scale_list_dubMutant[[23]] = c(0.04436,	0.04436)
scale_list_dubMutant[[24]] = c(0.04468,	0.04468)
scale_list_dubMutant[[25]] = c(0.03838,	0.03838)
scale_list_dubMutant[[26]] = c(0.04453,	0.04453)
scale_list_dubMutant[[27]] = c(0.043,	0.043)
scale_list_dubMutant[[28]] = c(0.043,	0.043)
scale_list_dubMutant[[29]] = c(0.043,	0.043)


image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/mutant/', full.names = T)
image_names <- image_names[grepl('tif',image_names)]


# Read in images and scale based on measurements provided

mutants <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'.tif')[[1]],'-')[[1]][2])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  return(im_ras)
})

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/wildType/', full.names = T)
image_names <- image_names[grepl('tif',image_names)]

wildTypes <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'#')[[1]][2],'_')[[1]][1])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  return(im_ras)
})

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/xi2k/', full.names = T)
image_names <- image_names[grepl('tif',image_names)]

doubleMutants <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'-')[[1]][2],'.tif'))
  im <- readTIFF(fn)
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_dubMutant[[imnum]][1],0,(nrow(im)-1)*scale_list_dubMutant[[imnum]][2])
  return(im_ras)
})



# Find the smallest extents in x&y and crop in a bit

xshrink <- .1
yshrink <- .05
threshold <- 1-.33

smallestX <- min(c(sapply(wildTypes, function(im) extent(im)[2]), sapply(mutants, function(im) extent(im)[2]))) * (1-xshrink)
smallestY <- min(c(sapply(wildTypes, function(im) extent(im)[4]), sapply(mutants, function(im) extent(im)[4]))) * (1-yshrink)


wildTypes_bboxes <- lapply(wildTypes, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
})

wildTypes_bboxes[[3]] <- wildTypes_bboxes[[3]] + c(20,20,0,0)
wildTypes_bboxes[[6]] <- wildTypes_bboxes[[6]] + c(0,0,3,3)
wildTypes_bboxes[[8]] <- wildTypes_bboxes[[8]] + c(-35,-35,0,0)
wildTypes_bboxes[[9]] <- wildTypes_bboxes[[9]] + c(15,15,-1,-1)
wildTypes_bboxes[[13]] <- wildTypes_bboxes[[13]] + c(-5,-5,0,0)
wildTypes_bboxes[[14]] <- wildTypes_bboxes[[14]] + c(0,0,-3,-3)
wildTypes_bboxes[[15]] <- wildTypes_bboxes[[15]] + c(0,0,-3,-3)
wildTypes_bboxes[[16]] <- wildTypes_bboxes[[16]] + c(15,15,0,0)
wildTypes_bboxes[[19]] <- wildTypes_bboxes[[19]] + c(10,10,0,0)
wildTypes_bboxes[[20]] <- wildTypes_bboxes[[20]] + c(-10,-10,0,0)


mutant_bboxes <- lapply(mutants, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
})

mutant_bboxes[[1]] <-  mutant_bboxes[[1]] + c(25,25,0,0)
mutant_bboxes[[2]] <-  mutant_bboxes[[2]] + c(25,25,0,0)
mutant_bboxes[[3]] <-  mutant_bboxes[[3]] + c(15,15,0,0)
mutant_bboxes[[6]] <-  mutant_bboxes[[6]] + c(25,25,0,0)
mutant_bboxes[[7]] <-  mutant_bboxes[[7]] + c(15,15,0,0)
mutant_bboxes[[9]] <-  mutant_bboxes[[9]] + c(30,30,0,0)
mutant_bboxes[[12]] <-  mutant_bboxes[[12]] + c(20,20,0,0)
mutant_bboxes[[20]] <-  mutant_bboxes[[20]] + c(0,0,-2,-2)
mutant_bboxes[[22]] <-  mutant_bboxes[[22]] + c(30,30,0,0)


doubleMutants_bboxes <- lapply(doubleMutants, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
})



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

doubleMutants_threshold <- lapply(doubleMutants, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  newExtent <- extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
  im <- raster::crop(im,newExtent)
  extent(im) <- newExtent
  return(im)
})





# Remove bad images?
# WT 10,11
# wildTypes_threshold <- wildTypes_threshold[-c(10,11)]

cells <- append(wildTypes_threshold,mutants_threshold)
cells <- append(cells,doubleMutants_threshold)


rawData_mats <- lapply(cells, function(ras){
  mat <- data.frame(as.matrix(coordinates(ras)))
  mat$value <- values(ras)
  return(mat)
})

set.seed(1)

resamples<-3
totalSampled <- resamples*1000

rawData_mats_3s <- lapply(rawData_mats, function(mat){
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


metadat <- data.frame(labels=c(rep('WT', length(wildTypes_threshold)*resamples),
                               rep('M', length(mutants_threshold)*resamples),
                               rep('DM', length(doubleMutants_threshold)*resamples)),
                      sampleNum=rep(1:resamples, length(cells)),
                      fold=NA,
                      networkNum=rep(1:length(cells),each=resamples))

set.seed(1)

folds <- createFolds(metadat$labels, k=5)
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

# summary_stats <- summary_stats[summary_stats$label%in%c('M','DM'),]

#Limit to 20 of each
# netNums <- c(1:20,21:40,43:62)
# summary_stats <- summary_stats[summary_stats$networkNum%in%netNums,]
# table(summary_stats$label)

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






######### PCA?
pcaTopDat <- prcomp(summary_stats[c('mean0.Death','sd0.Death','mean1.Birth','mean1.Death','sd1.Birth','sd1.Death')],scale=T,center = T)
summary_stats$PC1 <- pcaTopDat$x[,1]
summary_stats$PC2 <- pcaTopDat$x[,2]

ggplot(data=summary_stats) + geom_density(aes(PC1,fill=label),alpha=.5)
ggplot(data=summary_stats) + geom_density(aes(PC2,fill=label),alpha=.5)















