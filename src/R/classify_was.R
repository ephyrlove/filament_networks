library(TDA)
library(caret)
library(tictoc)
library(doParallel)
library(class)

setwd("C:/Users/ephy/Projects/cytoplasmic_streaming/data/original/actin_13092019/")


#save(was0, was1, file = "../../altered/wasserstein_distances_top20.Rdata")
load("../../altered/wasserstein_distances_top20.Rdata")

# ss: each class has 50 samples
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

#########################################################
# Dpc
load("../../altered/dpc_distances_15092019.Rdata")
dpc0.adj <- dpc0*floor((max(dpc1)/max(dpc0))) # Ask Le

dpcAcc <- kfoldCrossval(metadat@data, list(dpc0.adj, dpc1))
plot(w, dpcAcc, type='o', ylab='Accuracy')
max(dpcAcc)



#TODO: Make function takes as many distances dynamically



# Look at summary stats for dpc0 and dpc1

l <- lapply(1:3, function(class){
  m <-  dpc1[metadat@data$row[metadat@data$label==class],metadat@data$row[metadat@data$label==class]]
  diag(m) <- NA
  apply(m, MARGIN = 1, function(r) mean(r, na.rm=T))
})
hist(log(l[[1]]))
hist(log(l[[2]]))
hist(log(l[[3]]))

# Anova
cl <- do.call('c', l)
labs <- rep(1:3, each=length(l[[1]]))
boxplot(log(cl)~labs)


l <- lapply(1:3, function(class){
  m <-  dpc0[metadat@data$row[metadat@data$label==class],metadat@data$row[metadat@data$label==class]]
  diag(m) <- NA
  apply(m, MARGIN = 1, function(r) mean(r, na.rm=T))
})
hist(log(l[[1]]))
hist(log(l[[2]]))
hist(log(l[[3]]))

# Anova
cl <- do.call('c', l)
labs <- rep(1:3, each=length(l[[1]]))
boxplot(log(cl)~labs)



# The groups seem to to very diferent within groups, so maybe we should be using a more advanced classifier than ML on distance...





# Try clustering
k=6


#Function to get anchor-rows
rowDists <- function(rows){
  sqrt((dpc0.adj[rows,]^2)+(dpc1[rows,]^2))
}
  

groupSize=nrow(metadat@data)/k


#Random Instantiation for k-folds

#TODO: Actually need to create new nodes!!

randomRows <- sample(metadat@data$row, k, replace = F)

anchors <- rowDists(randomRows)

for(i in 1:100){
  distances <- do.call(rbind, lapply(seq(nrow(anchors)), function(a){
    data.frame(anchor=a, dist=anchors[a,], col=seq(ncol(anchors)))
  }))
  
  distances <- distances[order(distances$dist),]
  
  # Assign anchors
  allocations <- matrix(ncol = k, nrow = groupSize)
  for(row in split(distances, seq(distances$anchor))){
    if(!row$col %in% allocations){
      ind <- (sum(!is.na(allocations[,row$anchor]))+1)
      if(ind<51){
        allocations[ind,row$anchor] <- row$col
      }
    }
  }
  
  # Calculate new centroids
  new_anchors <- sapply(seq(k), function(class){
    limited_dists <- dpc0.adj[allocations[,class],allocations[,class]]
    d <- sapply(seq(groupSize), function(i){
      mean(limited_dists[i,])
    })
    which(d==min(d))
  })
  
  anchors <- rowDists(new_anchors)
}


metadat@data$final_assigned <- sapply(metadat@data$row, function(r) which(colSums(allocations==r)>0))

table(metadat@data$label, metadat@data$final_assigned)


#################################################


# Can I decompose back to an approximate euclidean system to cluster?

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

k=3

X1 <- dist_decomp(dpc1)
X0 <- dist_decomp(dpc0.adj)
X1_X0 <- dist_decomp((dpc0.adj*.8)+(dpc1*.2))

ws <- seq(0,1, by=.05)
acc <- sapply(ws, function(w){
  kmodel <- kmeans(cbind(X0*w,X1*(w-1)), k)
  t <- table(metadat@data$label, kmodel$cluster)
  colnames(t)[1] <- which(t[,1]==max(t[,1]))
  colnames(t)[2] <- which(t[,2]==max(t[,2]))
  colnames(t)[3] <- as.character(which(!(as.character(seq(1:k)) %in% colnames(t)[1:2])))
  
  correct=0
  for(c in 1:k){
    correct <- correct+t[as.character(c),as.character(c)]
  }
  correct/150
})

plot(ws, acc, type='o')



################################################

# Maybe k-NN makes more sense, since we need to use nodes that we have already?

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

X1 <- dist_decomp(dpc1)
X0 <- dist_decomp(dpc0.adj)
dpcAccs_knn <- lapply(w, function(w){
  X1_X0 <- dist_decomp((dpc0.adj*w)+(dpc1*(1-w)))
  
  lapply(list(X0,X1,X1_X0), function(data_comb){
    s <- sapply(unique(metadat@data$fold), function(fold){
      pr <- knn(data_comb[metadat@data$fold!=fold,],data_comb[metadat@data$fold==fold,],
                cl=metadat@data$label[metadat@data$fold!=fold],k=3)
      sum(pr==metadat@data[metadat@data$fold==fold,'label'])/length(pr)
    })
    sum(s)/length(s)
  })
})
dpcAccs_knn <- data.table::rbindlist(dpcAccs_knn)


kmod <- kmeans(X1_X0, centers = 3)
table(pred=kmod$cluster,truth=metadat@data$label)

X1 <- dist_decomp(was1)
X0 <- dist_decomp(was0)
wasAccs_knn <- lapply(w, function(w){
  
  X1_X0 <- dist_decomp((was0*w)+(was1*(1-w))) # WHAT??? These are perfect reflections???
  #print(paste(w, dim(X1_X0)))
  
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

X1 <- dist_decomp(was1)
X0 <- dist_decomp(was0)
wasAccs_svm <- lapply(w, function(w){
  
  X1_X0 <- dist_decomp((was0*w)+(was1*(1-w))) # WHAT??? These are perfect reflections???
  #print(paste(w, dim(X1_X0)))
  
  lapply(list(X0,X1,X1_X0), function(data_comb){
    s <- sapply(unique(metadat@data$fold), function(fold){
      mod <- svm(metadat@data$label[metadat@data$fold!=fold]~data_comb[metadat@data$fold!=fold,])
      preds <- predict(mod, data_comb[metadat@data$fold==fold,])
      sum(pr==metadat@data[metadat@data$fold==fold,'label'])/length(pr)
    })
    sum(s)/length(s)
  })
})
wasAccs_svm <- data.table::rbindlist(wasAccs_svm)

X1 <- dist_decomp(dpc1)
X0 <- dist_decomp(dpc0)
dpcAccs_svm <- lapply(w, function(w){
  
  X1_X0 <- dist_decomp((dpc0*w)+(dpc1*(1-w))) # WHAT??? These are perfect reflections???
  #print(paste(w, dim(X1_X0)))
  
  lapply(list(X0,X1,X1_X0), function(data_comb){
    s <- sapply(unique(metadat@data$fold), function(fold){
      mod <- svm(data=data_comb[metadat@data$fold!=fold,], metadat@data$label[metadat@data$fold!=fold]~.)
      preds <- predict(mod, data_comb[metadat@data$fold==fold,])
      sum(pr==metadat@data[metadat@data$fold==fold,'label'])/length(pr)
    })
    sum(s)/length(s)
  })
})
dpcAccs_svm <- data.table::rbindlist(dpcAccs_svm)


kmod <- kmeans(X1_X0, centers = 3)
table(pred=kmod$cluster,truth=metadat@data$label)

# Maybe that dpc retains more eigenvectors and values than was for some reason? Question for VM/Cassy


w=seq(0.05, 1, by=.05)

Accuracy <- c(dpcAcc,wasAcc,as.numeric(unlist(dpcAccs_knn[,3])),as.numeric(unlist(wasAccs_knn[,3])))
Distance <- rep(c('Dpc', 'Wasserstein', 'Dpc', 'Wasserstein'), each=length(w))
ClusterMethod <- rep(c('Traditional', 'Traditional', 'k-NN', 'k-NN'), each=length(w))


ggplot() + 
  geom_line(aes(x=rep(w, times=4) , y=Accuracy, color=Distance, linetype=ClusterMethod), size=2) +
  geom_abline(aes(intercept=.95, slope=0), linetype=2)+
  xlab('w') + ylim(c(.6,1))




