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




mSkew <- sapply(mutants,function(x) skewness(values(x)))
wtSkew <- sapply(wildTypes,function(x) skewness(values(x)))

ggplot() + geom_density(aes(mSkew,fill=rep(c('first','additional'),c(22,28))), alpha=.5) + ggtitle('Mean of mutant images by sample') 
ggplot() + geom_density(aes(wtSkew,fill=rep(c('first','additional'),c(20,30))), alpha=.5)+ ggtitle('Mean of wt images by sample') + labs(fill='Sample') + xlab('Mean')


ggplot() + geom_density(aes(c(mSkew[1:22],wtSkew[1:20]),fill=rep(c('mutant','wildtype'),c(22,20))), alpha=.5) + ggtitle('Skewness in first sample images by class')
ggplot() + geom_density(aes(c(mSkew[23:50],wtSkew[21:50]),fill=rep(c('mutant','wildtype'),c(28,30))), alpha=.5) + ggtitle('Skewness in second sample images by class')


ggplot() + geom_boxplot(aes(c(mSkew[1:22],wtSkew[1:20]),fill=rep(c('mutant','wildtype'),c(22,20))), alpha=.5) + ggtitle('Skewness in first sample images by class')
ggplot() + geom_boxplot(aes(c(mSkew[23:50],wtSkew[21:50]),fill=rep(c('mutant','wildtype'),c(28,30))), alpha=.5) + ggtitle('Skewness in second sample images by class')


ggplot() + geom_boxplot(aes(mSkew,fill=rep(c('first','additional'),c(22,28))), alpha=.5) + ggtitle('Skewness of mutant images by sample') 
ggplot() + geom_boxplot(aes(wtSkew,fill=rep(c('first','additional'),c(20,30))), alpha=.5)+ ggtitle('Skewness of wt images by sample') + labs(fill='Sample') + xlab('Skewness')





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



mSkew <- sapply(rawData_mats[51:100],function(x) skewness(x$value))
wtSkew <- sapply(rawData_mats[1:50],function(x) skewness(x$value))


ggplot() + geom_density(aes(mSkew,fill=rep(c('first','additional'),c(22,28))), alpha=.5)
ggplot() + geom_density(aes(wtSkew,fill=rep(c('first','additional'),c(20,30))), alpha=.5)




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


plot(mutants[[10]])
plot(mutants_threshold[[10]]>.05)



















library(imager)


BiocManager::install("EBImage")



library(EBImage)


r <- wildTypes[[10]]

im <- EBImage::as.Image(as.matrix(r))

im <-medianFilter(im,5)

display(im)


f <- matrix(c(1,0,-1,2,0,-2,1,0,1),nrow=3,byrow = t)
filt <- as.Image(f)
filt_t <- as.Image(t(f))
filtered <- filter2(im,filt)
filtered_t <- filter2(im,filt_t)
fullFilter <- sqrt(filtered^2 + filtered_t^2)
display(fullFilter)



fullFilterMat <- as.matrix(fullFilter@.Data)
fullFilterMat <- (fullFilterMat-min(fullFilterMat))/(max(fullFilterMat)-min(fullFilterMat))

newr <- raster(fullFilterMat)
newr@extent <- extent(r)
plot(r)
plot(newr)


plot()


edges <- deriche(im,10,order=2,axis="y") 
plot(edges)



la = matrix(1, nc=3, nr=3)
la[2,2] = -6
y = filter2(im, la)
display(y, title='Filtered image')
display(medianFilter(y,5))



rcim <- as.cimg(as.matrix(r))
plot(rcim)


rcimHough <- hough_line(rcim,data.frame = F)
plot(rcimHough)






plot(rasterToContour(as.raster(fullFilterMat)))



px <- px.square(30,80,80) %>% boundary
plot(px)
#Hough transform
hough_line(px,ntheta=200) %>% plot

df <- hough_line(px,ntheta=800,data.frame=TRUE)
#Plot lines with the highest score
plot(px)
with(subset(df,score > quantile(score,.95)),nfline(theta,rho,col="red"))

plot(boats)
df <- hough_line(boats,ntheta=800,data=TRUE)
# }

with(subset(df,score > quantile(score,.999)),nfline(theta,rho,col="red"))














library(seriation)

absDiff <- function(matrix1,matrix2)
{
  r <- nrow(matrix1)
  c <- ncol(matrix1)
  destMatrix <- matrix1
  for(r in 0:r-1)
  {
    for(c in 0:c-1)
    {
      destMatrix[r,c] <- abs(matrix1[r,c]-matrix1[r,c])
    }
  }
  return(destMatrix)
}

countNonZero <- function(inputMatrix)
{
  return(length(inputMatrix[inputMatrix > 0]))
}

thinningIteration <- function(imageMatrix, iter)
{
  imageInput <- imageMatrix
  r <- nrow(imageInput) - 1
  c <- ncol(imageInput) - 1
  for(i in 2:r)
  {
    for(j in 2:c)
    {
      p2 <- imageInput[i-1, j]
      p3 <- imageInput[i-1, j+1]
      p4 <- imageInput[i, j+1]
      p5 <- imageInput[i+1, j+1]
      p6 <- imageInput[i+1, j]
      p7 <- imageInput[i+1, j-1]
      p8 <- imageInput[i, j-1]
      p9 <- imageInput[i-1, j-1]
      A  <- (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) + 
        (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) + 
        (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
        (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1)
      B  <- p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
      if(iter == 0){
        m1 <- (p2 * p4 * p6)
        m2 <- (p4 * p6 * p8)
      }
      else {
        m1 <- (p2 * p4 * p8)
        m2 <- (p2 * p6 * p8)
      }
      if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
      {
        imageInput[i,j] <- 0
      }
    }
  }
  return(imageInput)
}

thinImage <- function(imageMatrix)
{
  
  im <- imageMatrix
  prev <- im
  repeat {
    im <- thinningIteration(im, 0)
    im <- thinningIteration(im, 1)
    diff <- absDiff(im, prev)
    prev <- im
    if(countNonZero(diff) <= 0)
    {
      break
    }
  } 
  
  return(im)
}



r <- wildTypes[[1]]

im <-as.cimg(r)
blurred <- imager::blur_anisotropic(im,amplitude=1e4)
plot(blurred)




singleImageMatrix <- as.matrix(wildTypes_threshold[[1]])



singleImageMatrix  <- singleImageMatrix>quantile(as.numeric(singleImageMatrix),.8)

pimage(singleImageMatrix)

thin <- thinImage(singleImageMatrix)

pimage(thin)
plot(wildTypes[[1]])





