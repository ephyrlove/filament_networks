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

all_image_names <- image_names

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/wildType/', full.names = T)
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


# Remove bad images?
# WT 10,11
# wildTypes_threshold <- wildTypes_threshold[-c(10,11)]

cells <- append(wildTypes_threshold,mutants_threshold)

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


metadat <- data.frame(labels=c(rep('WT', length(wildTypes_threshold)*resamples),rep('M', length(mutants_threshold)*resamples)),
                      sampleNum=rep(1:resamples, 42),
                      fold=NA,
                      networkNum=rep(1:length(cells),each=resamples),
                      fileName=rep(sapply(strsplit(all_image_names, '//'),'[',2),each=resamples))
                      
set.seed(2)
# folds <- createFolds(paste(metadat$labels,metadat$networkNum), k=5)
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
summary_stats$fileName <- metadat$fileName
summary_stats$sampleNum <- metadat$sampleNum


# write.csv(summary_stats,'~/Projects/cytoplasmic_streaming/data/altered/vectorizedPDsForAN.csv',row.names = F)






xValidReuslts <- c()
fnames <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats[summary_stats$fold!=fold,], as.factor(label) ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
  test <- summary_stats[summary_stats$fold==fold,]
  test$preds <- predict(svmMod, test)
  results <- c()
  for(net in unique(test$networkNum)){
    t <- table(test[test$networkNum==net,'preds'])
    results <- c(results,test$label[test$networkNum==net][1]==names(which(t==max(t))))
  }
  
  fnames <- c(fnames,unique(test[test$networkNum%in%which(!results),'fileName']))
  
  xValidReuslts <- c(xValidReuslts,sum(results)/length(results))
}
summary(xValidReuslts)
kable(data.frame(table(fnames)), format='rst')



### Outputting files for Dr. N ###
test <- summary_stats
svmMod <- svm(data=summary_stats, as.factor(label) ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
test$preds <- predict(svmMod, summary_stats)
results <- c()
for(net in unique(summary_stats$networkNum)){
  t <- table(test[test$networkNum==net,'preds'])
  results <- c(results,test$label[test$networkNum==net][1]==names(which(t==max(t))))
}

unique(test[test$networkNum%in%which(!results),'fileName'])


mutantVecs <- read.csv('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/mutant-morphology-210203.csv')
wtVecs <- read.csv('~/Projects/cytoplasmic_streaming/data/altered/wt_mutant_summary_data/WT-morphology-210203.csv')
cellData <- rbind(mutantVecs,wtVecs)


summary_stats_cellDat <- merge(x=summary_stats,y=cellData, by.x='fileName',by.y='File',all.x=T)

xValidReuslts <- c()
fnames <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats_cellDat[summary_stats_cellDat$fold!=fold,], as.factor(label) ~ mean0.Death+sd0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death+
                  Branch.Frequency+Mean.Angle+Parallelness+Occupancy+Filament.Density+Coefficient.of.Variation+Skewness)
  test <- summary_stats_cellDat[summary_stats_cellDat$fold==fold,]
  test$preds <- predict(svmMod, test)
  results <- c()
  for(net in unique(test$networkNum)){
    t <- table(test[test$networkNum==net,'preds'])
    results <- c(results,test$label[test$networkNum==net][1]==names(which(t==max(t))))
  }
  
  fnames <- c(fnames,unique(test[test$networkNum%in%which(!results),'fileName']))
  
  xValidReuslts <- c(xValidReuslts,sum(results)/length(results))
}
summary(xValidReuslts)
kable(data.frame(table(fnames)), format='rst')

xValidReuslts <- c()
fnames <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats_cellDat[summary_stats_cellDat$fold!=fold,], as.factor(label) ~ Branch.Frequency+Mean.Angle+Parallelness+Occupancy+Filament.Density+Skewness)
  test <- summary_stats_cellDat[summary_stats_cellDat$fold==fold,]
  test$preds <- predict(svmMod, test)
  results <- c()
  for(net in unique(test$networkNum)){
    t <- table(test[test$networkNum==net,'preds'])
    results <- c(results,test$label[test$networkNum==net][1]==names(which(t==max(t))))
  }
  
  fnames <- c(fnames,unique(test[test$networkNum%in%which(!results),'fileName']))
  
  xValidReuslts <- c(xValidReuslts,sum(results)/length(results))
}
summary(xValidReuslts)
kable(data.frame(table(fnames)), format='rst')




summary_stats_cellDat_1 <- summary_stats_cellDat[summary_stats_cellDat$sampleNum==1,]

summary_stats_cellDat_1pd <- summary_stats_cellDat_1[c("mean.Birth","mean.Death","sd.Birth","sd.Death","mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")]
summary_stats_cellDat_1vec <- summary_stats_cellDat_1[c("Branch.Frequency","Mean.Angle","Parallelness","Occupancy","Filament.Density","Skewness")]

pcaPD <- prcomp(summary_stats_cellDat_1pd, center = TRUE, scale. = FALSE)
pcaVec <- prcomp(summary_stats_cellDat_1vec, center = TRUE, scale. = FALSE)

summary_stats_cellDat_1$pdPC1 <- pcaPD$x[,1]
summary_stats_cellDat_1$pdPC2 <- pcaPD$x[,2]
summary_stats_cellDat_1$vecPC1 <- pcaVec$x[,1]
summary_stats_cellDat_1$vecPC2 <- pcaVec$x[,2]

plot(pcaPD)
plot(pcaVec)


plot(summary_stats_cellDat_1$pdPC1,summary_stats_cellDat_1$vecPC1)
plot(summary_stats_cellDat_1$pdPC1,summary_stats_cellDat_1$vecPC2)
plot(summary_stats_cellDat_1$pdPC2,summary_stats_cellDat_1$vecPC1)
plot(summary_stats_cellDat_1$pdPC2,summary_stats_cellDat_1$vecPC2)


summary(glm(summary_stats_cellDat_1$label=='WT'~summary_stats_cellDat_1$pdPC1+summary_stats_cellDat_1$pdPC2+summary_stats_cellDat_1$vecPC1+summary_stats_cellDat_1$vecPC2,family='binomial'))


summary(lm(summary_stats_cellDat_1$pdPC1~summary_stats_cellDat_1$vecPC1+summary_stats_cellDat_1$vecPC2))
summary(lm(summary_stats_cellDat_1$pdPC2~summary_stats_cellDat_1$vecPC2))
summary(lm(summary_stats_cellDat_1$vecPC1~summary_stats_cellDat_1$pdPC1))




summary_stats$pc1 <- princomp(summary_stats[c("mean.Birth","mean.Death","sd.Birth","sd.Death","mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,1]
summary_stats$pc2 <- princomp(summary_stats[c("mean.Birth","mean.Death","sd.Birth","sd.Death","mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,2]

xValidReuslts <- c()
fnames <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats[summary_stats$fold!=fold,], as.factor(label) ~ pc1+pc2)
  test <- summary_stats[summary_stats$fold==fold,]
  test$preds <- predict(svmMod, test)
  results <- c()
  for(net in unique(test$networkNum)){
    t <- table(test[test$networkNum==net,'preds'])
    results <- c(results,test$label[test$networkNum==net][1]==names(which(t==max(t))))
  }
  
  fnames <- c(fnames,unique(test[test$networkNum%in%which(!results),'fileName']))
  
  xValidReuslts <- c(xValidReuslts,sum(results)/length(results))
}
summary(xValidReuslts)
kable(data.frame(table(fnames)), format='rst')


summary_stats_cellDat$pc1 <- princomp(summary_stats_cellDat[c("mean.Birth","mean.Death","sd.Birth","sd.Death","mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,1]
summary_stats_cellDat$pc2 <- princomp(summary_stats_cellDat[c("mean.Birth","mean.Death","sd.Birth","sd.Death","mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,2]
summary_stats_cellDat$pc1v <- princomp(summary_stats_cellDat[c("Branch.Frequency","Mean.Angle","Parallelness","Occupancy","Filament.Density","Skewness")])$scores[,1]
summary_stats_cellDat$pc2v <- princomp(summary_stats_cellDat[c("Branch.Frequency","Mean.Angle","Parallelness","Occupancy","Filament.Density","Skewness")])$scores[,2]


xValidReuslts <- c()
fnames <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats_cellDat[summary_stats_cellDat$fold!=fold,], as.factor(label) ~ pc1+pc2+pc1v)
  test <- summary_stats_cellDat[summary_stats_cellDat$fold==fold,]
  test$preds <- predict(svmMod, test)
  results <- c()
  for(net in unique(test$networkNum)){
    t <- table(test[test$networkNum==net,'preds'])
    results <- c(results,test$label[test$networkNum==net][1]==names(which(t==max(t))))
  }
  
  fnames <- c(fnames,unique(test[test$networkNum%in%which(!results),'fileName']))
  
  xValidReuslts <- c(xValidReuslts,sum(results)/length(results))
}
summary(xValidReuslts)
kable(data.frame(table(fnames)), format='rst')





dbrdaMod <- vegan::dbrda(summary_stats_cellDat_1pd~.,data=summary_stats_cellDat_1vec)
plot(dbrdaMod)

biplot(pcaVec)
biplot(pcaPD)

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


png('~/Projects/cytoplasmic_streaming/fig/pcPlot.png',height = 5,width=6, res=300,units = 'in')
ggplot() + geom_point(aes(x=summary_stats_cellDat_1$pdPC1,y=summary_stats_cellDat_1$pdPC2,color=summary_stats_cellDat_1$label)) +
  xlab('PD Principal Component 1') + ylab('PD Principal Component 2') + labs(color = element_text('Treatment'))
dev.off()


png('~/Projects/cytoplasmic_streaming/fig/wss.png',height = 6,width=6, res=300,units = 'in')
wssplot(summary_stats_cellDat_1[,c('pdPC2','pdPC1')])
dev.off()


kmod <- kmeans(summary_stats_cellDat_1[,c('pdPC2','pdPC1')],centers=2)


png('~/Projects/cytoplasmic_streaming/fig/pcPlot_clustered.png',height = 5,width=6, res=300,units = 'in')
ggplot() + geom_point(aes(x=summary_stats_cellDat_1$pdPC1,y=summary_stats_cellDat_1$pdPC2,color=summary_stats_cellDat_1$label,shape=as.factor(kmod$cluster))) + geom_point(aes(y=kmod$centers[,1],x=kmod$centers[,2]),size=2,shape=8)+
  xlab('PD Principal Component 1') + ylab('PD Principal Component 2') + labs(color = element_text('Treatment'),shape=element_text('Cluster'))
dev.off()

allCells <- append(wildTypes,mutants)
c1 <- as.numeric(which(kmod$cluster==1))
c2 <- as.numeric(which(kmod$cluster==2))

t <- table(kmod$cluster,summary_stats_cellDat_1$label)
clusterAssigments <- names(apply(t, MARGIN = 2,function(r)which(r==max(r))))

wrongC1 <- kmod$cluster==1 & summary_stats_cellDat_1$label!=clusterAssigments[1]
wrongC2 <- kmod$cluster==2 & summary_stats_cellDat_1$label!=clusterAssigments[2]


m<-matrix(seq(1:44),ncol=2,byrow = T)
m[21,1] <-0
m[21,2] <-41
m[22,1] <-0
m[22,2] <-42
png('~/Projects/cytoplasmic_streaming/fig/cluster_xamples.png', width = 6,height=10,units = 'in', res=150)
layout(mat = m)
par(mar=c(0.1,0.1,0.1,0.1))
for(i in seq(22)){
  if(i<=20){
    image(allCells[[c1[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
    if(wrongC1[c1[i]]){
      mp <- as.numeric(bbox(allCells[[c1[i]]])[,2]/2)
      points(y=mp[2],x=mp[1],pch=4,cex=10,col='red')
    } 
  }
  image(allCells[[c2[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  if(wrongC2[c2[i]]){
    mp <- as.numeric(bbox(allCells[[c2[i]]])[,2]/2)
    points(y=mp[2],x=mp[1],pch=4,cex=10,col='red')
  }
}
dev.off()



sum(diag(table(kmod$cluster,summary_stats_cellDat_1$label)))/length(kmod$cluster)






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














