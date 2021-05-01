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


plot(wildTypes[[1]], legend=F,col=colorRamp(c('white','black')))

raster_dfs <- lapply(wildTypes, function(r) data.frame(coordinates(r),i=values(r)))


png('~/Projects/cytoplasmic_streaming/fig/Examples.png', width = 7,height=11,units = 'in', res=150)
par(mfrow=c(10,2),mar=c(0.1,0.1,0.1,0.1))
for(i in seq(10)){
  image(mutants[[i]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  image(wildTypes[[i]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
}
dev.off()


allCells <- append(wildTypes,mutants)

c1 <- as.numeric(which(kmod$cluster==1))
c2 <- as.numeric(which(kmod$cluster==2))

png('~/Projects/cytoplasmic_streaming/fig/cluster_xamples.png', width = 7,height=11,units = 'in', res=150)
par(mfrow=c(10,2),mar=c(0.1,0.1,0.1,0.1))
for(i in seq(10)){
  image(allCells[[c1[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  image(allCells[[c2[i]]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
}
dev.off()



set.seed(1)


png('~/Projects/cytoplasmic_streaming/fig/randomized_examples.png', width = 7,height=11,units = 'in', res=150)
par(mfrow=c(10,2),mar=c(0.1,0.1,0.1,0.1))
for(i in seq(10)){
  if(runif(1)>.5){
    image(mutants[[i]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
    image(wildTypes[[i]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  }else{
    image(wildTypes[[i]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
    image(mutants[[i]], col = gray.colors(10, start = 0, end = 1, gamma = 2, alpha = NULL),yaxt="n",xaxt='n')
  }
}
dev.off()








