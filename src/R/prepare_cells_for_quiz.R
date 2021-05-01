options(stringsAsFactors = F)
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

mutants <- mutants[1:20]

image_names <- dir('~/Projects/cytoplasmic_streaming/data/original/microscopy/wildType/', full.names = T)
image_names <- image_names[grepl('tif',image_names)]

wildTypes <- lapply(image_names, function(fn){
  imnum <- as.numeric(strsplit(strsplit(fn,'#')[[1]][2],'_')[[1]][1])
  im <- readTIFF(fn)
  im_ras <- raster(im)
  extent(im_ras)<- c(0,(ncol(im)-1)*scale_list_mutant[[imnum]][1],0,(nrow(im)-1)*scale_list_mutant[[imnum]][2])
  return(im_ras)
})


all_cells <- c(wildTypes, mutants)

set.seed(1)
random_ids = sample(1:40,40,replace = F)

codesheet <- data.frame(vidId=random_ids, uid=1:40)
codesheet$lable <- rep(c('wt','mutant'),each=20)

setwd('~/Projects/cytoplasmic_streaming/data/altered/quiz cells/')
for(uid in codesheet$uid){
  visId <- codesheet$vidId[codesheet$uid==uid]
  tiff(paste0(visId,'.tif'),height=3,width=5,res = 300,units = 'in')
  raster::plot(all_cells[[uid]],axes=F,legend=F,box=F)
  dev.off()
}


tiff('all_images.tif',height=20,width=40,res = 80,units = 'in')
par(mfrow=c(5,8), mar=c(0,0,0,0))
for(vid in 1:40){
  uid <- codesheet$uid[codesheet$vidId==vid]
  raster::plot(all_cells[[uid]],axes=F,legend=F,box=F)
  text(x=3,y=3,labels=vid,size=4)
}
dev.off()

responses <- read.csv('~/Projects/cytoplasmic_streaming/data/original/responses.csv')
names(responses)[3:42] <- paste0('q',1:40)

accs <- c()
cellBio <- c()
for(i in 1:nrow(responses)){
  answers <- responses[i,paste0('q',1:40)]
  answers <- as.numeric(factor(answers,levels=c('Option 1', 'Option 2')))
  t <- table(codesheet$lable[order(codesheet$vidId)],answers)
  d1 <- sum(c(t[1,1],t[2,2]))/40
  d2 <- sum(c(t[1,2],t[2,1]))/40
  accs <- c(accs,max(d1,d2))
  cellBio <- c(cellBio, responses$Do.you..have.experience.with.microscopy.of.cells.[i])
}
mean(accs)
sd(accs)

