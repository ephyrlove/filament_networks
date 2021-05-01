library(TDA)
library(caret)
library(tictoc)
library(doParallel)
library(class)

setwd("~/Projects/cytoplasmic_streaming/data/original/actin_13092019/")

rawdata = list() # txt file, coordinate of actins
alldata = list() # betta 0 and 1 dim features

DiagLim <- 10
maxdimension <- 1


cl <- makeCluster(40)  
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

set.seed(1)
rawData_sampled <- lapply(rawData, function(dat) dat[sample(seq(1100), 220, replace = F),])
sapply(rawData_sampled, nrow)

alldata <- foreach(x=rawData_sampled) %dopar% computeRips(x)

loopOver <- function(dat, i, dim){
  c(rep(0, i),
    sapply((i+1):150, function(j){
      TDA::wasserstein(dat[[i]],dat[[j]],p=2,dimension =dim)^(0.5)
    }))
}

alldata <- lapply(alldata, '[[', 'diagram')

tic()
a <- foreach(i=seq(149)) %dopar% loopOver(alldata, i, 0)
b <- foreach(i=seq(149)) %dopar% loopOver(alldata, i, 1 )
toc()

# This 10% sample took <10min to run, which is pretty nice
# This 20% sample took <2.6hrs to run, which is not so nice

was0 <- as.matrix(rbind(do.call(rbind, a), rep(0, 150)))
was1 <- as.matrix(rbind(do.call(rbind, b), rep(0, 150)))
was0=t(was0)+was0
was1=t(was1)+was1

save(was0, was1, file='../../altered/wasserstein_distances_sample20.Rdata')
#load('../../altered/wasserstein_distances_sample20.Rdata')

ss<-50
actualClassSize <- 3

metadat = data.frame(row=seq(actualClassSize*ss), label=as.factor(rep(c(1,2,3),each=ss)), fold=NA)

# Creates a distance_metadata class object
createDistanceMetadat <- function(labels){
  metadat <- data.frame(row=seq(labels), label=labels, fold=NA)
  setClass('distance_metadata', slots=list(data='data.frame'))
  metadat <- new('distance_metadata', data=metadat)
  return(metadat)
}

metadat <- createDistanceMetadat(as.factor(rep(1:actualClassSize, each=ss)))


#Set folds for 10-fold crossval
set.seed(20)
Nf = 10
folds <- createFolds(metadat@data$label,k = Nf)

for(i in 1:length(folds)){
  metadat@data$fold[folds[[i]]] <- i
}


dist_decomp <- function(distmat){
  # Graham decomp
  M <- matrix(nrow=nrow(distmat),ncol=ncol(distmat))
  for(i in 1:nrow(M)){
    for(j in 1:ncol(M))
      M[i, j] <- 0.5*(distmat[1, j]^2 + distmat[i, 1]^2 - distmat[i, j]^2)
  }
  
  # M = USUt; X=U*sqrt(S)
  # Find unitary and diagonal matrix of eigenvalues
  
  e <- eigen(M)
  U = e$vectors
  s = e$values
  
  U <- U[,which(s>0)]
  s <- s[which(s>0)]
  
  
  smat <- matrix(0, nrow = length(s), ncol=length(s))
  diag(smat) <- s
  return(U%*%sqrt(smat))
}

X1 <- dist_decomp(was1)
X0 <- dist_decomp(was0)
w=seq(0.05, 1, by=.05)
wasAccs_knn <- lapply(w, function(w){
  
  X1_X0 <- dist_decomp((was0*w)+(was0*(w-1)))

  lapply(list(X0,X1,X1_X0), function(data_comb){
    s <- sapply(unique(metadat@data$fold), function(fold){
      pr <- knn(data_comb[metadat@data$fold!=fold,],data_comb[metadat@data$fold==fold,],
                cl=metadat@data$label[metadat@data$fold!=fold],k=3)
      sum(pr==metadat@data[metadat@data$fold==fold,'label'])/length(pr)
    })
    sum(s)/length(s)
  })
})
wasAccs_knn <- data.table::rbindlist(wasAccs_knn)






kfoldCrossval <- function(dist_info, distances){
  
  distances[[1]]->dist0
  distances[[2]]->dist1
  
  # Check that all distqance matrices are square and of same dim
  dims <- unlist(lapply(distances, dim))
  if(any(dims!=dims[1])) stop('All Distance Matrices must be square and of equal dimension')
  
  w=seq(0.05, 1, by=.05)
  dt <- data.table::rbindlist(lapply(w, function(w){
    lapply(seq(Nf),w=w, function(i, w){
      train <- dist_info[dist_info$fold!=i,]
      test <- dist_info[dist_info$fold==i,]
      
      dist0_train <- dist0[train$row,]
      dist1_train <- dist1[train$row,]
      
      pred <- sapply(seq(nrow(test)), function(i){
        item <- test[i,]
        sums <- colSums(rbind(w*sapply(split(dist0_train[,item[['row']]],train$label), mean),
                              (1-w)*sapply(split(dist1_train[,item[['row']]],train$label), mean)))
        as.numeric(which(sums==min(sums)))
      })
      
      sum(pred==test$label)/length(pred)
    })
  }))
  
  apply(dt, MARGIN = 1, mean)
}


##########################################################
# Wassersein
w=seq(0.05, 1, by=.05)
wasAcc <- kfoldCrossval(metadat@data, list(was0, was1))
plot(w, wasAcc, type='o', ylab='Accuracy')
lines(w,unlist(wasAccs_knn[,3]))





