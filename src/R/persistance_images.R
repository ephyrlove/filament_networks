library(TDA)
library(caret)
library(tictoc)
library(doParallel)
library(sp)
library(raster)
library(e1071)
library(data.table)
library(imager)

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


r <- raster()
extent(r) <- c(0,20,0,20)
ncol(r) <- 80
nrow(r) <- 80

rasterizeML <- function(ml, r){
  ras <- raster::rasterize(ml, r, fun='sum')
  raster::values(ras)[is.na(raster::values(ras))] <- 0
  raster::values(ras) <- raster::values(ras)/max(raster::values(ras))
  return(ras)
}

cl <- makeCluster(40)
registerDoParallel(cl)

rawData_rasters <- foreach(ml=rawData_multiLines) %dopar% rasterizeML(ml, r)

for(i in seq(rawData_rasters)){
  img <- as.cimg(rawData_rasters[[i]])
  save.image(img, file = sprintf('../../altered/filament_images/filament_images/net%s.png',i))
}


