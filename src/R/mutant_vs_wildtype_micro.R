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
  im <- as.im(as.matrix(im))
  im <- spatstat::blur(im, sigma = 4)
  im <- raster(as.matrix(im))
  extent(im) <- newExtent
  values(im) <- as.numeric(values(im)>=quantile(values(im),threshold)[1])
  return(im)
})

mutants_threshold <- lapply(mutants, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  newExtent <- extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
  im <- raster::crop(im,newExtent)
  im <- as.im(as.matrix(im))
  im <- spatstat::blur(im, sigma = 4)
  im <- raster(as.matrix(im))
  extent(im) <- newExtent
  values(im) <- as.numeric(values(im)>=quantile(values(im),threshold)[1])
  return(im)
})




# May be that we need to identify the "center of entropy" or some such to focus on, rather than the spatial center...
# Also maybe method to enhance the filamentness

cells <- append(wildTypes_threshold,mutants_threshold)

rawData_mats <- lapply(cells, function(ras){
  as.matrix(coordinates(ras)[values(ras)==1,])
})
summary(sapply(rawData_mats, nrow))

set.seed(1)
rawData_mats_3s <- lapply(rawData_mats, function(mat){
  samp <- sample(seq(nrow(mat)), size = 3000, replace=F)
  list(s1=mat[samp[1:1000],],s2=mat[samp[1001:2000],],s3=mat[samp[2001:3000],])
})

computeRips <- function(data) TDA::ripsDiag(data, 1, 10)

cl <- makeCluster(20)
registerDoParallel(cl)

tic()
rips_s1 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s1)
rips_s2 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s2)
rips_s3 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s3)
toc()



allRips <- c(rips_s1, rips_s2, rips_s3)

metadat <- data.frame(labels=rep(c(rep('WT', length(wildTypes)), rep('M', length(mutants))),3), sampleNum=rep(1:3, each=length(rips_s1)), fold=NA, networkNum = rep(1:length(cells),3))

set.seed(1)
folds <- createFolds(metadat$labels, k=5)
for(i in seq(folds)){
  metadat$fold[metadat$row%in%folds[[i]]] <- i
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

ind <- which(metadat$sampleNum==1)
round(mean(sapply(folds, function(fold){
  svmMod <- svm(data=summary_stats[-fold,], as.factor(label) ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
  test <- summary_stats[fold[fold%in%ind],]
  sum(predict(svmMod, test)==test$label)/sum(fold%in%ind)*100
})),2)






##################### Custom bounding boxes ##########################
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

wildTypes_threshold <- lapply(seq(wildTypes), function(i){
  im <- wildTypes[[i]]
  newExtent <- wildTypes_bboxes[[i]]
  im <- raster::crop(im,newExtent)
  im <- as.im(as.matrix(im))
  im <- spatstat::blur(im, sigma = 4)
  im <- raster(as.matrix(im))
  extent(im) <- newExtent
  values(im) <- as.numeric(values(im)>=quantile(values(im),threshold)[1])
  return(im)
})

mutants_threshold <- lapply(seq(mutants), function(i){
  im <- mutants[[i]]
  newExtent <- mutant_bboxes[[i]]
  im <- raster::crop(im,newExtent)
  im <- as.im(as.matrix(im))
  im <- spatstat::blur(im, sigma = 4)
  im <- raster(as.matrix(im))
  extent(im) <- newExtent
  values(im) <- as.numeric(values(im)>=quantile(values(im),threshold)[1])
  return(im)
})

# Remove bad images?
# WT 10,11
# wildTypes_threshold <- wildTypes_threshold[-c(10,11)]

cells <- append(wildTypes_threshold,mutants_threshold)

rawData_mats <- lapply(cells, function(ras){
  as.matrix(coordinates(ras)[values(ras)==1,])
})
summary(sapply(rawData_mats, nrow))

set.seed(1)
rawData_mats_3s <- lapply(rawData_mats, function(mat){
  samp <- sample(seq(nrow(mat)), size = 3000, replace=F)
  list(s1=mat[samp[1:1000],],s2=mat[samp[1001:2000],],s3=mat[samp[2001:3000],])
})

computeRips <- function(data) TDA::ripsDiag(data, 1, 10)

cl <- makeCluster(20)
registerDoParallel(cl)

tic()
rips_s1 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s1)
rips_s2 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s2)
rips_s3 <- foreach(x=rawData_mats_3s) %dopar% computeRips(x$s3)
toc()



allRips <- c(rips_s1, rips_s2, rips_s3)

metadat <- data.frame(labels=rep(c(rep('WT', length(wildTypes_threshold)), rep('M', length(mutants_threshold))),3), sampleNum=rep(1:3, each=42), fold=NA, networkNum = rep(1:length(cells),3))

set.seed(1)
folds <- createFolds(metadat$labels, k=5)
for(i in seq(folds)){
  metadat$fold[metadat$row%in%folds[[i]]] <- i
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

ind <- which(metadat$sampleNum==3)
round(mean(sapply(folds, function(fold){
  svmMod <- svm(data=summary_stats[-fold,], as.factor(label) ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
  test <- summary_stats[fold[fold%in%ind],]
  sum(predict(svmMod, test)==test$label)/sum(fold%in%ind)*100
})),2)



###### Try DPc #######
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
    sapply((i+1):120, function(j){
      TDA::wasserstein(dat[[i]]$diagram,dat[[j]]$diagram,p=2,dimension =dim)^(0.5)
    }))
}

cl <- makeCluster(30)  
registerDoParallel(cl) 


tic()
a <- foreach(i=seq(119)) %dopar% loopOver(top20Rips, i, 0)
b <- foreach(i=seq(119)) %dopar% loopOver(top20Rips, i, 1 )
toc()

# save(a,b, file = '~/Projects/cytoplasmic_streaming/data/altered/micros_wasD_top20_3s.Rdata')
load('~/Projects/cytoplasmic_streaming/data/altered/micros_wasD_top20_3s.Rdata')


was0 <- do.call(rbind, a)
was1 <- do.call(rbind, b)
was0 <- rbind(was0, rep(0,length(a[[1]])))
was1 <- rbind(was1, rep(0,length(b[[1]])))

was0 <- t(was0)+was0
was1 <- t(was1)+was1


### Test wasserstein 3xsampled distance
netNums <- unique(metadat$networkNum)
was0_means <- matrix(0,nrow = 120, ncol=120)
for(i in netNums){
  rn <- as.numeric(row.names(metadat[metadat$networkNum==i,]))
  for(j in netNums){
    rn2 <- as.numeric(row.names(metadat[metadat$networkNum==j,]))
    was0_means[i,j] <- mean(was0[rn2,rn])
  }
}
diag(was0_means) <- 0

was1_means <- matrix(0,nrow = 120, ncol=120)
for(i in netNums){
  rn <- as.numeric(row.names(metadat[metadat$networkNum==i,]))
  for(j in netNums){
    rn2 <- as.numeric(row.names(metadat[metadat$networkNum==j,]))
    was1_means[i,j] <- mean(was1[rn2,rn])
  }
}
diag(was1_means) <- 0

# Distance from one sample to one class
# Input: item - label of the sample, class - No. of the class, dismatrix - distance matrix
Dis_to_class <- function(item,class,dismatrix,classTrain){
  if(class == 1){
    average_dist=mean(dismatrix[item,classTrain[[1]]])
  }
  if(class == 2){
    average_dist=mean(dismatrix[item,classTrain[[2]]])
  }
  return(average_dist)
}

# classify the sample
# Input: item, w=weight of 0 feature
classify <- function(item,w,dist_0,dist_1,classTrain){
  Dis <- matrix(0,2,2)
  for(i in 1:2){
    Dis[1,i] <- w*Dis_to_class(item,i,dist_0,classTrain)
    Dis[2,i] <- (1-w)*Dis_to_class(item,i,dist_1,classTrain)
  }
  bsum <- colSums(Dis)
  label <- which(bsum == min(bsum))[1]
  return(label)
}

# cross validation: 9 folds vs  1 fold
# returns vector of accuracy rates
kFoldClassify <- function(classify, dist_0, dist_1, actuallabels, Nf){
  # create 10 folds
  set.seed(1)
  folds <- createFolds(actuallabels,k = Nf)
  
  accrate <- rep(NA, 21)
  for (j in 0:20) {
    dataIndex <- seq(actuallabels)
    testLabels <- rep(NA, length(actuallabels))
    for (i in 1:Nf) {
      Testing  <- dataIndex[folds[[i]]]
      Training <- dataIndex[-folds[[i]]]
      c1Training <- Training[which(actuallabels[-folds[[i]]]=='WT')]
      c2Training <- Training[which(actuallabels[-folds[[i]]]=='M')]
      testLabels[folds[[i]]] <- sapply(Testing,w=0.05*j,dist_0=dist_0, dist_1=dist_1,
                                       classTrain=list(c1Training,c2Training),classify)
    }
    testLabels <- factor(testLabels, labels=c('WT','M'))
    accrate[j+1]=sum(testLabels==actuallabels)/length(actuallabels)
  }
  return(accrate)
}

accrate <- kFoldClassify(classify, dist_0=was0_means, dist_1=was1_means, metadat$labels, Nf=5)
max(accrate)
plot(seq(0,1, by=.05), accrate, type='o', ylab='Accuracy Rate', xlab='Dim-0 Weight')






############################## Comparing with Andreas' thresholded images

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy_filtered/WT/4-filaments/', full.names = T)
image_names <- image_names[grepl('tif',image_names)]


wildTypes_T <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'.tif')[[1]],'-')[[1]][2])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  return(im_ras)
})

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy_filtered/xik/3-masks/', full.names = T)
image_names <- image_names[grepl('tif',image_names)]

mutants_T <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'.tif')[[1]][1],'509-')[[1]][2])
  im <- as.matrix(readTIFF(fn))
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  return(im_ras)
})














