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






image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/mutant/', full.names = T)
image_names <- c(image_names,dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/xik mutant 28 additional pics/', full.names = T))

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
  values(im_ras) <- (values(im_ras)-min(values(im_ras)))/(max(values(im_ras))-min(values(im_ras)))
  return(im_ras)
})

all_image_names <- image_names

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/wildType/', full.names = T)
image_names <- c(image_names,dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/WT 30 additional pics/', full.names = T))
image_names <- image_names[grepl('tif',image_names)]

all_image_names <- c(image_names,all_image_names)

wildTypes <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'#')[[1]][2],'_')[[1]][1])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  values(im_ras) <- (values(im_ras)-min(values(im_ras)))/(max(values(im_ras))-min(values(im_ras)))
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
  
  im <- EBImage::as.Image(as.matrix(im))
  
  im <-medianFilter(im,5)
  
  f <- matrix(c(1,0,-1,2,0,-2,1,0,1),nrow=3,byrow = T)
  filt <- as.Image(f)
  filt_t <- as.Image(t(f))
  filtered <- filter2(im,filt)
  filtered_t <- filter2(im,filt_t)
  fullFilter <- sqrt(filtered^2 + filtered_t^2)
  
  
  fullFilterMat <- as.matrix(fullFilter@.Data)
  fullFilterMat <- (fullFilterMat-min(fullFilterMat))/(max(fullFilterMat)-min(fullFilterMat))
  
  newr <- raster(fullFilterMat)
  newr@extent <- extent(newExtent)
  
  return(newr)
})


mutants_threshold <- lapply(mutants, function(im){
  imMidPoint <- (as.vector(extent(im))/2)[c(2,4)]
  newExtent <- extent(imMidPoint[1] - smallestX/2 , imMidPoint[1] + smallestX/2, imMidPoint[2] - smallestY/2, imMidPoint[2] + smallestY/2)
  im <- raster::crop(im,newExtent)
  extent(im) <- newExtent
  
  im <- EBImage::as.Image(as.matrix(im))
  
  im <-medianFilter(im,5)
  
  f <- matrix(c(1,0,-1,2,0,-2,1,0,1),nrow=3,byrow = T)
  filt <- as.Image(f)
  filt_t <- as.Image(t(f))
  filtered <- filter2(im,filt)
  filtered_t <- filter2(im,filt_t)
  fullFilter <- sqrt(filtered^2 + filtered_t^2)
  
  
  fullFilterMat <- as.matrix(fullFilter@.Data)
  fullFilterMat <- (fullFilterMat-min(fullFilterMat))/(max(fullFilterMat)-min(fullFilterMat))
  
  newr <- raster(fullFilterMat)
  newr@extent <- extent(newExtent)
  
  return(newr)
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



summary_stats$pc1 <- princomp(summary_stats[c("mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,1]
summary_stats$pc2 <- princomp(summary_stats[c("mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")])$scores[,2]




xValidReuslts <- c()
for(fold in unique(metadat$fold)){
  svmMod <- svm(data=summary_stats[summary_stats$fold!=fold,], as.factor(label) ~ pc1+pc2)
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


summary_stats_1perCell$batch<-factor(NA,levels=c('Batch 1','Batch 2'))

summary_stats_1perCell[summary_stats_1perCell$label=='M','batch'] <- ifelse(summary_stats_1perCell$networkNum[summary_stats_1perCell$label=='M']<74,'Batch 1','Batch 2')
summary_stats_1perCell[summary_stats_1perCell$label=='WT','batch'] <- ifelse(summary_stats_1perCell$networkNum[summary_stats_1perCell$label=='WT']<21,'Batch 1','Batch 2')




ggplot(data=summary_stats_1perCell[summary_stats_1perCell$batch=='Batch 1',]) + geom_point(aes(x=pc1,y=pc2,color=label))
ggplot(data=summary_stats_1perCell[summary_stats_1perCell$batch=='Batch 2',]) + geom_point(aes(x=pc1,y=pc2,color=label))

ggplot(data=summary_stats_1perCell) + geom_point(aes(x=pc1,y=pc2,color=batch))


ggplot(data=summary_stats_1perCell) + geom_point(aes(x=pc1,y=pc2,color=label,size=batch),alpha=.8)




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

length(explanatoryCols)

params[]


ggplot(data=params) + geom_point(aes(x=Mean.Angle,y=Skewness,color=label))


pvs <- sapply(explanatoryCols,function(c) t.test(params[params$label=='M',c],params[params$label=='WT',c])$p.value)

pvs_adj <- p.adjust(pvs,method = 'bonferroni', n = length(pvs))

sum(pvs_adj<.05)

1-((1-.05)^32) # Family-wise error rate



pca <- prcomp(params[explanatoryCols],center = T,scale. = T)

params$pc1 <- pca$x[,1]
params$pc2 <- pca$x[,2]
params$pc3 <- pca$x[,3]

ggplot(data=params) + geom_point(aes(x=pc1,y=pc2,color=label,size=summary_stats_1perCell$batch))




t.test(params[params$label=='M','pc1'],params[params$label=='WT','pc1'])
t.test(params[params$label=='M','pc2'],params[params$label=='WT','pc2'])

t.test(summary_stats_1perCell[summary_stats_1perCell$label=='M','pc1'],summary_stats_1perCell[summary_stats_1perCell$label=='WT','pc1'])
t.test(summary_stats_1perCell[summary_stats_1perCell$label=='M','pc2'],summary_stats_1perCell[summary_stats_1perCell$label=='WT','pc2'])


sumExplanCols <- c("mean0.Death" ,"sd0.Death","mean1.Birth","mean1.Death","sd1.Birth","sd1.Death")

summary_stats_andParams <- data.frame(params[explanatoryCols],summary_stats_1perCell[c(sumExplanCols,'label')])



# Loading the library
library(glmnet)

x_vars <- model.matrix(label~. , summary_stats_andParams)[,-1]
y_var <- as.numeric(summary_stats_andParams$label=='M')
lambda_seq <- 10^seq(2, -2, by = -.1)

# Splitting the data into test and train
set.seed(1)
train = sample(1:nrow(x_vars), nrow(x_vars)*.75)
x_test = (-train)
y_test = y_var[x_test]

cv_output <- cv.glmnet(x_vars[train,], y_var[train],
                       alpha = 1, lambda = lambda_seq, 
                       nfolds = 5, family='binomial')

# identifying best lamda
best_lam <- cv_output$lambda.min
best_lam


lasso_best <- glmnet(x_vars[train,], y_var[train], alpha = 1, lambda = best_lam, family = 'binomial')
pred <- predict(lasso_best, s = best_lam, newx = x_vars[x_test,],type='response')

mean(y_test==(pred>.5))

coef(lasso_best)



fullPCAmod <- prcomp(summary_stats_andParams[c(sumExplanCols,explanatoryCols)],center = T,scale. = T)

plot(fullPCAmod)

kmod <- kmeans(x=fullPCAmod$x[,1:7],centers = 3)

t <- table(kmod$cluster,summary_stats_andParams$label, summary_stats_1perCell$batch)
t
sum(t[2]+t[3])/sum(t)






summary_stats_andParams <- data.frame(params[explanatoryCols],summary_stats_1perCell[c(sumExplanCols,'batch')])


x_vars <- model.matrix(label~. , summary_stats_andParams)[,-1]
y_var <- as.numeric(summary_stats_andParams$batch=='Batch 1')
lambda_seq <- 10^seq(2, -2, by = -.1)

# Splitting the data into test and train
set.seed(1)
train = sample(1:nrow(x_vars), nrow(x_vars)*.75)
x_test = (-train)
y_test = y_var[x_test]

cv_output <- cv.glmnet(x_vars[train,], y_var[train],
                       alpha = 1, lambda = lambda_seq, 
                       nfolds = 5, family='binomial')

# identifying best lamda
best_lam <- cv_output$lambda.min
best_lam


lasso_best <- glmnet(x_vars[train,], y_var[train], alpha = 1, lambda = .09, family = 'binomial')
pred <- predict(lasso_best, s = best_lam, newx = x_vars[x_test,],type='response')
mean(y_test==(pred>.5))


