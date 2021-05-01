library(TDA)
library(caret)
library(tictoc)
library(doParallel)
library(sp)
library(raster)
library(e1071)
library(gridExtra)

setwd("~/Projects/cytoplasmic_streaming/data/original/actin_07102019/")

loadData <- function(pattern){
  files <- dir(pattern, full.names = T)
  return(lapply(files, read.table))
}

rawData <- lapply(c("825cp/", "1650cp/", "3300cp/"), loadData)
rawData <- do.call('c', rawData)
sapply(rawData, nrow)

labels <- as.factor(rep(c('825','1650','3300'), each=50))

load('../../altered/filament_images/linesAndRasters.Rdata')

###################
# Write out points to csvs for ripser
# sample_out <- min(sapply(rawData_rasters, function(r) sum(values(r)==1))) # take same cardinality based on lowest #
# 
# set.seed(1)
# rawData_points <- lapply(rawData_rasters, function(r){
#   coord <- coordinates(r)[which(values(r)==1),]
#   coord[sample(1:nrow(coord),size = sample_out, replace = F),]
#   })
# 
# # tic()
# # pd <- TDA::ripsDiag(rawData_points[[1]], maxdimension = 1, maxscale = 10)
# # toc()
# # 7 min and got bad_alloc error
# 
# for(i in seq(rawData_points)){
#   write.csv(rawData_points[[i]], file = sprintf('../../altered/subsampled_actin_pcs/pc_%s.csv', i), row.names = F)
# }

files <- dir('../../altered/sampled_pds/', full.names = T)

all_rips_3000 <- lapply(files, readLines)

alldata <- lapply(all_rips_3000, function(pd){
  head_positions <- which(grepl("persistence intervals in dim", pd))+1
  pd0 <- pd[seq(from=head_positions[1], to=(head_positions[2]-2))]
  pd1 <- pd[seq(from=head_positions[2], to=length(pd))]
  pd0 <- data.frame(do.call(rbind, strsplit(gsub('\\)', '', sapply(strsplit(pd0, '\\['), '[', 2)), ',')))
  pd1 <- data.frame(do.call(rbind, strsplit(gsub('\\)', '', sapply(strsplit(pd1, '\\['), '[', 2)), ',')))
  names(pd0) <- c('Birth', 'Death')
  names(pd1) <- c('Birth', 'Death')
  pd0$betti <- 0
  pd1$betti <- 1
  pd_out <- rbind(pd0, pd1)
  pd_out$Birth <- as.numeric(pd_out$Birth)
  pd_out$Death <- as.numeric(pd_out$Death)
  pd_out$Death[is.na(pd_out$Death)] <- 10 # [0, ) should impute threshold value...
  pd_out <- pd_out[,c('betti', 'Birth', 'Death')]
  return(pd_out)
})

metadat <- data.frame(labels=labels, row=seq(labels), fold=NA)

set.seed(3)
folds <- createFolds(metadat$labels, k=10)
for(i in seq(folds)){
  metadat$fold[metadat$row%in%folds[[i]]] <- i
}

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

summary_stats <- summary_stats[,!names(summary_stats)%in%c('mean0.Birth', 'sd0.Birth')] # THese are constants that should probably just be removed...

summary_stats$label <- metadat$labels

mean(sapply(folds, function(fold){
  svmMod <- svm(data=summary_stats[-fold,], label ~ mean0.Death+mean1.Birth+mean1.Death+sd1.Birth+sd1.Death)
  test <- predict(svmMod, summary_stats[fold,])
  sum(test==summary_stats[fold,'label'])/length(fold)
}))





#### Distances
alldata <- lapply(all_rips_3000, function(pd){
  head_positions <- which(grepl("persistence intervals in dim", pd))+1
  pd0 <- pd[seq(from=head_positions[1], to=(head_positions[2]-2))]
  pd1 <- pd[seq(from=head_positions[2], to=length(pd))]
  pd0 <- data.frame(do.call(rbind, strsplit(gsub('\\)', '', sapply(strsplit(pd0, '\\['), '[', 2)), ',')))
  pd1 <- data.frame(do.call(rbind, strsplit(gsub('\\)', '', sapply(strsplit(pd1, '\\['), '[', 2)), ',')))
  names(pd0) <- c('Birth', 'Death')
  names(pd1) <- c('Birth', 'Death')
  pd0$dimension <- 0
  pd1$dimension <- 1
  pd_out <- rbind(pd0, pd1)
  pd_out$Birth <- as.numeric(pd_out$Birth)
  pd_out$Death <- as.numeric(pd_out$Death)
  pd_out$Death[is.na(pd_out$Death)] <- 10 # [0, ) should impute threshold value...
  pd_out <- pd_out[,c('dimension', 'Birth', 'Death')]
  return(as.matrix(pd_out))
})


loopOver <- function(dat, i, dim){
  c(rep(0, i),
    sapply((i+1):150, function(j){
      TDA::wasserstein(dat[[i]],dat[[j]],p=2,dimension =dim)^(0.5)
    }))
}


tic()
a <- foreach(i=seq(149)) %dopar% loopOver(alldata, i, 0)
b <- foreach(i=seq(149)) %dopar% loopOver(alldata, i, 1 )
toc()

