#Returns the dpc distance for persistence disgrams d1 and d2 for Betti number beta

library('clue')

schudist <- function(d1,d2, beta,p,c){
  d1 = as.matrix(d1)
  d2 = as.matrix(d2)
  
  wd1 = d1[which(d1[,1]==beta),2:3, drop = FALSE]
  wd2 = d2[which(d2[,1]==beta),2:3, drop = FALSE]
  n = dim(wd1)[1]
  m = dim(wd2)[1]
  
  if(n == 0 || m == 0){
    if(n == 0 && m == 0){
      return(0)
    }
    return(c)
  }
  
  #always make wd1 smaller
  if(m < n){
    dum1 = n
    dum2 = wd1
    n = m
    m = dum1
    wd1 = wd2
    wd2 = dum2
  }
  
  #Create distance matrix
  diffs <- matrix(abs(wd1[rep(seq(n),each=m),] - wd2[rep(seq(m),n),]), ncol=2) # matrix wrapper because we have to make sure this does not reduce down to vector in stupid R world
  distmat <- matrix(apply(diffs, MARGIN = 1, function(r) min(c,max(r))^p),nrow=n,ncol=m, byrow = T) # Takes infnorm
  
  solution <- as.numeric(clue::solve_LSAP(distmat)) #Uses HUngarian Algorithm
  cost <- sum(diag(distmat[seq(n),solution]))
  cost <- ((cost + (c^p)*(m-n))*(1/m))^(1/p)
  return(cost)
}



# Make a version for testing c's without rerunning if not neccessary...
# takes vector of c's and stops testing when cost does not increase!
schudist_testC <- function(d1,d2, beta,p,c){
  d1 = as.matrix(d1)
  d2 = as.matrix(d2)
  
  wd1 = d1[which(d1[,1]==beta),2:3, drop = FALSE]
  wd2 = d2[which(d2[,1]==beta),2:3, drop = FALSE]
  n = dim(wd1)[1]
  m = dim(wd2)[1]
  
  if(n == 0 || m == 0){
    if(n == 0 && m == 0){
      return(0)
    }
    return(c)
  }
  
  #always make wd1 smaller
  if(m < n){
    dum1 = n
    dum2 = wd1
    n = m
    m = dum1
    wd1 = wd2
    wd2 = dum2
  }
  
  #Create distance matrix
  diffs <- abs(wd1[rep(seq(n),each=m),] - wd2[rep(seq(m),n),])
  
  # We don't need to keep calculating costs if max(r) is greater than c for all cases. We could hash also hash solutions,
  # but probably not worth the memory and io
  costs <- rep(NA, length(c))
  for(i in seq(c)){
    distmat <- matrix(apply(diffs, MARGIN = 1, function(r) min(c[i],max(r))^p),nrow=n,ncol=m, byrow = T) # Takes infnorm
    solution <- as.numeric(clue::solve_LSAP(distmat)) #Uses HUngarian Algorithm
    cost <- sum(diag(distmat[seq(n),solution]))
    costs[i] <- ((cost + (c[i]^p)*(m-n))*(1/m))^(1/p)
    if(i>1 && i<length(c)){
      if(costs[i]==costs[i-1]){
        costs[(i+1):length(c)] <- costs[i]
        break()
      }
    }
  }
  
  return(costs)
}











