library(TDA)
library(caret)
library(tictoc)
library(doParallel)

setwd("~/Projects/cytoplasmic_streaming/data/original/actin_13092019/")

rawdata = list() # txt file, coordinate of actins
alldata = list() # betta 0 and 1 dim features

DiagLim <- 10
maxdimension <- 1


cl <- makeCluster(detectCores()-14)  
registerDoParallel(cl)  

computeRips <- function(data) TDA::ripsDiag(data, maxdimension,DiagLim)

loadData <- function(i, pattern){
  filename=paste0(pattern,toString(i),".txt")
  return(read.table(filename))
}


rawData <- foreach(i=seq(50)) %dopar% loadData(i, "actins_case_100_11_2.0625_")
rawData <- c(rawData, foreach(i=seq(50)) %dopar% loadData(i, "actins_case_100_11_4.125_"))
rawData <- c(rawData, foreach(i=seq(50)) %dopar% loadData(i, "actins_case_100_11_8.25_"))
sapply(rawData, nrow)

tic()
alldata <- foreach(x=rawData) %dopar% computeRips(x)
toc()

## Ran in 79 min

# save(alldata, file='../../altered/alldata_13092019.Rdata')
load('../../altered/alldata_13092019.Rdata')

# cut the top 20% of topological features
cuttopv2 <- function(pd,p0,p1){
    bds=pd[["diagram"]]
    Dime=attr(pd[["diagram"]],"dimension")
    Maxd=attr(pd[["diagram"]],"maxdimension")
    Scal=attr(pd[["diagram"]],"scale")
    Call=attr(pd[["diagram"]],"call")
    
    bds0=bds[bds[,1]==0,]
    bds1=bds[bds[,1]==1,]
    pers0=bds0[,3]-bds0[,2]
    pers1=bds1[,3]-bds1[,2]
    
    threshold0=quantile(pers0,p0) 
    bds0=bds0[pers0>=threshold0,]
    
    threshold1=quantile(pers1,p1) 
    bds1=bds1[pers1>=threshold1,]
    
    bds=rbind(bds0,bds1)
    pd[["diagram"]]=bds
    
    class(pd[["diagram"]])="diagram"
    attr(pd[["diagram"]],"dimension")=Dime
    attr(pd[["diagram"]],"maxdimension")=Maxd
    attr(pd[["diagram"]],"scale")=Scal
    attr(pd[["diagram"]],"call")=Call
    return(pd)
}

newdata <- lapply(alldata, function(d) cuttopv2(d,0.8,0.8))

#####################################################################
## compute the Wasserstein distance matrix
   
#load("../../altered/wasserstein_distances_top20_second.Rdata")

# ss: each class has 50 samples
ss<-50
actuallabels <- as.factor(rep(c(1,2,3),each=ss))

# Distance from one sample to one class
# Input: item - label of the sample, class - No. of the class, dismatrix - distance matrix
Dis_to_class <- function(item,class,dismatrix,classTrain){
  if(class == 1){
    average_dist=mean(dismatrix[item,classTrain[[1]]])
  }
  if(class == 2){
    average_dist=mean(dismatrix[item,classTrain[[2]]])
  }
  if(class == 3){
    average_dist=mean(dismatrix[item,classTrain[[3]]])
  }
  return(average_dist)
}

# classify the sample
# Input: item, w=weight of 0 feature
classify <- function(item,w,dist_0,dist_1,classTrain){
  Dis <- matrix(0,2,3)
  for(i in 1:3){
    Dis[1,i] <- w*Dis_to_class(item,i,dist_0,classTrain)
    Dis[2,i] <- (1-w)*Dis_to_class(item,i,dist_1,classTrain)
  }
  bsum <- colSums(Dis)
  label <- which(bsum == min(bsum))
  return(label)
}

# cross validation: 9 folds vs  1 fold
# returns vector of accuracy rates
kFoldClassify <- function(classify, dist_0, dist_1, actuallabels, ss, Nf){
  # create 10 folds
  set.seed(1)
  folds <- createFolds(actuallabels,k = Nf)
  
  accrate <- rep(NA, 21)
  for (j in 0:20) {
    dataIndex <- seq(actuallabels)
    testLabels <- rep(NA, length(actuallabels))
    Ts <- ss-(ss/Nf)
    for (i in 1:Nf) {
      Testing  <- dataIndex[folds[[i]]]
      Training <- dataIndex[-folds[[i]]]
      c1Training <- Training[Training<=ss]
      c2Training <- Training[Training>ss & Training<=(ss*2)]
      c3Training <- Training[Training>(ss*2)]
      testLabels[folds[[i]]] <- sapply(Testing,w=0.05*j,dist_0=dist_0, dist_1=dist_1,
                                       classTrain=list(c1Training,c2Training,c3Training),classify)
    }
    testLabels <- factor(testLabels, levels=levels(actuallabels))
    accrate[j+1]=sum(testLabels==actuallabels)/length(actuallabels)
  }
  return(accrate)
}

accrate <- kFoldClassify(classify, dist_0=was0, dist_1=was1, actuallabels, ss=50, Nf=10)
plot(0:20, accrate, type='o', ylab='Accuracy Rate', xlab='Threshold')


######################################################################
# compute the dpc distance matrix
source("../../../src/R/dpc.R")

nd_diagrams <- lapply(newdata, '[[', 'diagram')
ss <- 50

computeRow <- function(i, newdata, beta, p, c, ss){
  return(c(rep(0, i), lapply((i+1):(ss*3), function(j){
    schudist(nd_diagrams[[i]],nd_diagrams[[j]],beta=beta,p,c)
  })))
}

# cl <- makeCluster(6)
# registerDoParallel(cl)
#642
p=2
c=0.8
r0 <- foreach(i=seq(149)) %dopar% computeRow(i, nd_diagrams, 0, p, c, ss)
r1 <- foreach(i=seq(149)) %dopar% computeRow(i, nd_diagrams, 1, p, c, ss)

dpc0 <- as.matrix(data.table::rbindlist(r0))
dpc1 <- as.matrix(data.table::rbindlist(r1))
dpc0 <- rbind(dpc0, rep(0, ss*3))
dpc1 <- rbind(dpc1, rep(0, ss*3))

dpc0 <- t(dpc0)+dpc0
dpc1 <- t(dpc1)+dpc1

#save(dpc0, dpc1, file = "../../altered/dpc_distances_15092019.Rdata")
load("../../altered/dpc_distances_15092019.Rdata")

dpc0.adj <- dpc0*floor((max(dpc1)/max(dpc0))) # Ask Le

accrate <- kFoldClassify(classify, dist_0=dpc0.adj, dist_1=dpc1, actuallabels, 50, 10)

plot(0:20, accrate, type='o', ylab='Accuracy Rate', xlab='Threshold')


#############################################################
# Run this on AWS if you want it done this year.
load('alldata_13092019.Rdata')
newdata <- lapply(alldata, function(d) cuttopv2(d,0.8,0.8))

#Just copy and paste in computeRow_testC for now

library(clue)
library(tictoc)
library(doParallel)


cl <- makeCluster(71)
registerDoParallel(cl)

computeRow_testC <- function(i, diagrams, beta, p, cReg, ss){
  js <- (i+1):(ss*3)
  l <- lapply(js, function(j){
    schudist_testC(diagrams[[i]],diagrams[[j]],beta=beta,p,cReg)
  })
  costmat <- matrix(NA,ncol=(length(js)+i),nrow=length(cReg))
  for(x in seq(length(cReg))){
    costmat[x,] <- c(rep(0, i),as.numeric(lapply(l, '[', x)))
  }
  return(costmat)
}

cReg <- c(seq(.1,.9,by=.1),1:6)
nd_diagrams <- lapply(newdata, '[[', 'diagram')

ss <- 50
p=2
tic()
r0 <- foreach(i=seq(149)) %dopar% computeRow_testC(i, nd_diagrams, 0, p, cReg, ss)
toc()
tic()
r1 <- foreach(i=seq(149)) %dopar% computeRow_testC(i, nd_diagrams, 1, p, cReg, ss)
toc()

# Swtich to json!!!
#save(r0, r1, file = 'dpc_ctest_lists_0.6to6_18092019.Rdata')

################################################################
r0 <- jsonlite::read_json('../../altered/dpc_ctest_lists_0.6to6_18092019_r0.json')
r1 <- jsonlite::read_json('../../altered/dpc_ctest_lists_0.6to6_18092019_r1.json')


dpcByC <- lapply(1:15, function(c){
  beta0 <- rbind(as.matrix(data.table::rbindlist(lapply(1:149, function(i) r0[[i]][[c]]))), rep(0, ss*3))
  beta1 <- rbind(as.matrix(data.table::rbindlist(lapply(1:149, function(i) r1[[i]][[c]]))), rep(0, ss*3))
  
  beta0 <- t(beta0)+beta0
  beta1 <- t(beta1)+beta1
  
  return(list(beta0=beta0,beta1=beta1))
})

names(dpcByC) <- paste('c', cReg,sep = '_')


# Test out on scaled dpc distence over c from 0.1 to 1

# create 10 folds
# set.seed(20)
Nf=10
folds = createFolds(actuallabels,k = Nf)

accrate=matrix(0,1,41)
overall.accrate=matrix(0,1,15)


# overall.accrate <- data.frame(c=c(),threshold=c(),accrate=c())
# for (i in seq(cReg)) {
#   c <- cReg[i]
#   dpc0 <- dpcByC[[paste('c', c, sep='_')]]$beta0
#   dpc1 <- dpcByC[[paste('c', c, sep='_')]]$beta1
#   dpc0.adj=dpc0*floor((max(dpc1)/max(dpc0)))
# 
#   overall.accrate <- rbind(overall.accrate,data.frame(
#     c=rep(c,21),
#     threshold=seq(21),
#     accrate=kFoldClassify(classify, dist_0=dpc0.adj, dist_1=dpc1, actuallabels, 50, 10)
#   ))
# }

overall.accrate <- matrix(NA, nrow=60,ncol=21)
for (i in seq(cReg)) {
  c <- cReg[i]
  dpc0 <- dpcByC[[paste('c', c, sep='_')]]$beta0
  dpc1 <- dpcByC[[paste('c', c, sep='_')]]$beta1
  dpc0.adj=dpc0*floor((max(dpc1)/max(dpc0)))
  overall.accrate[61 - (c*10),]=kFoldClassify(classify, dist_0=dpc0.adj, dist_1=dpc1, actuallabels, 50, 10)
}


accRas <- raster::raster(overall.accrate)
extent(accRas) <- extent(c(1,21,0,6))

par(mfrow=c(2,2))
plot(accRas, xlab='w', ylab='Regularization - C', main='No Cut',xaxt='n')
axis(side=1, labels=seq(0,20, by=5)/20, at=seq(1,21, by=5))
plot(accRas>.85, xlab='w', ylab='Regularization - C',legend=FALSE, main='Cut at .85',xaxt='n')
axis(side=1, labels=seq(0,20, by=5)/20, at=seq(1,21, by=5))
plot(accRas>.86, xlab='w', ylab='Regularization - C',legend=FALSE,  main='Cut at .86',xaxt='n')
axis(side=1, labels=seq(0,20, by=5)/20, at=seq(1,21, by=5))
plot(accRas>.87, xlab='w', ylab='Regularization - C',  main='Cut at .87',xaxt='n')
axis(side=1, labels=seq(0,20, by=5)/20, at=seq(1,21, by=5))



# overall.accrate

# Test out on scaled dpc distence over c from 1 to 6
# Nf=10
# folds = createFolds(actuallabels,k = Nf)

accrate=matrix(0,1,41)
# overall.accrate=matrix(0,1,6)

for (k in 2:6) {
  c=k
  filename=paste0("dpcp",toString(2),"c",toString(c),"_top20percent.RData")
  load(filename)
  
  dpc0.adj=dpc0*floor((max(dpc1)/max(dpc0)))
  
  for (j in 0:40) {
    datalabels = as.integer(1:150)
    testlabels = integer(3*ss)
    Ts=ss-ss/Nf
    for (i in 1:10) {
      Testing  <- datalabels[folds[[i]]]
      Training = datalabels[-folds[[i]]]
      c1Training = Training[1:Ts]
      c2Training = Training[(Ts+1):(2*Ts)]
      c3Training = Training[(2*Ts+1):(3*Ts)]
      testlabels[folds[[i]]] <- lapply(Testing,w=0.025*j,Classify)
    }
    
    testlabels = as.factor(unlist(testlabels))
    confusion <- confusionMatrix(data = testlabels,actuallabels)
    accrate[j+1]=confusion$overall[1]
  }
  
  overall.accrate[k+9]=max(accrate)
}
overall.accrate
# [,1]      [,2] [,3] [,4]      [,5]      [,6]      [,7]      [,8]      [,9] [,10]     [,11]
# [1,] 0.8666667 0.8933333 0.88 0.88 0.8666667 0.8733333 0.8733333 0.8733333 0.8733333  0.88 0.8733333
# [,12]     [,13] [,14] [,15]
# [1,] 0.8733333 0.8666667  0.86  0.86

# [,1]      [,2]      [,3] [,4]      [,5]      [,6]      [,7] [,8]      [,9]     [,10] [,11] [,12]
# [1,] 0.86 0.8866667 0.8666667 0.88 0.8666667 0.8666667 0.8666667 0.86 0.8733333 0.8733333  0.88  0.88
# [,13]     [,14]     [,15]
# [1,] 0.8666667 0.8666667 0.8666667



