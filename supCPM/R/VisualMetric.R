#' Visualization Quality Metrics
#'
#'
#' @param X matrix; The coordinates of the high-dimensional data.
#' @param Y matrix; The coordinates of the embedding data.
#' @param label vector; The label vector indicating which cluster that each cell belongs to.
#' @param KncNeighbor integer; The size of nearest neighbor for KNC metric
#' @param KnnNeighbor integer; The size of nearest neighbor for KNN metric.
#'
#' @return Returns a list contains the following elements:
#' \item{CorV}{Spearman correlation between within-cluster total variation in the high-dimensional space and those in the low.}
#' \item{KNC}{The fraction of k-nearest neighbors of cluster means in the high dimensional data which are still preserved as the k-nearest cluster means in the embedding space.}
#' \item{CPD}{Spearman correlation between pairwise distances in the high and low dimensions.}
#' \item{CSM}{Cluster Separation Measure.}
#' \item{KNN}{The fraction of k-nearest neighbor pairs in the high dimensional data that are still$k-nearest neighbor pairs in the embedding dimensions.}
#' @export
#'
VisualMetric <- function(Y,X,label,KncNeighbor,KnnNeighbor){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Cov_high <- vector()
  Cov_low <- vector()
  CMean_high <- matrix(nrow=length(unique(label)),ncol=ncol(X))
  CMean_low <- matrix(nrow=length(unique(label)),ncol=ncol(Y))

  # CorV
  k = 1
  for( i in unique(label) ){
      cl <- which(label==i)
      Cov_high[k] <- sum(diag(cov(Y[cl,])))
      Cov_low[k] = sum(diag(cov(X[cl,])))
      k = k+1
  }
  CorV = cor(Cov_high,Cov_low, method='spearman')

  # KNC
  k = 1
  for ( i in unique(label) ){
    cl <- which(label==i)
    CMean_high[k,] <- colMeans(X[cl,])
    CMean_low[k,] <- colMeans(Y[cl,])
    k<-k+1
  }
  CNeighbor_high <- find.nn(as.matrix(dist(CMean_high)),KncNeighbor)$neighbors[,-1]
  CNeighbor_low <- find.nn(as.matrix(dist(CMean_low)),KncNeighbor)$neighbors[,-1]
  Cpreserved <- 0
  for ( i in 1:length(unique(label)) ){
    Cpreserved <- Cpreserved + length(intersect(CNeighbor_low[i,],CNeighbor_high[i,]))
  }
  KNC <- Cpreserved / (KncNeighbor*length(unique(label)))

  # CPD
  Dist_high <- as.vector(dist(X))
  Dist_low <- as.vector(dist(Y))
  CPD <- cor(Dist_high,Dist_low,method='spearman')

  # KNN
  Neighbor_high <- find.nn(as.matrix(dist(X)),KnnNeighbor)$neighbors[,-1]
  Neighbor_low <- find.nn(as.matrix(dist(Y)),KnnNeighbor)$neighbors[,-1]
  preserved <- 0
  n <- nrow(X)
  for ( i in 1:n ){
    preserved <- preserved + length(intersect(Neighbor_low[i,],Neighbor_high[i,]))
  }
  KNN <- preserved / (KnnNeighbor*n)

  # logFisher
#  L_tilde <- matrix(rep(0,n^2),n)
#  for (i in label){
#    cl = which(label==i)
#    L_tilde[cl,cl] = 1/length(cl)
#  }
# I = diag(1,n)
#  one = matrix(1,n,1)
#  Sb = (L_tilde-one%*%t(one)/n)
#  Sw = (I - L_tilde)
#  logFisher = log(1+sum(diag(t(Y)%*%Sw%*%Y))/sum(diag(t(Y)%*%Sb%*%Y)))

  # Cluster Separation
  CS = metric(Y,label,q=1,p=2)

  return(list(CorV=CorV,KNC=KNC,CPD=CPD,KNN=KNN,CS=CS))
}


metric = function(y,label,q=1,p=2){
  cl = length(unique(label))
  center = matrix(0,cl,dim(y)[2])
  R = matrix(0,cl,cl)
  M = matrix(0,cl,cl)
  S = rep(0,cl)
  cluster = unique(label)
  for(i in 1:cl){
    center[i,] = colMeans(y[label==cluster[i],])
  }
  for(i in 1:cl){
    #S[i] = mean(abs(y[label==cluster[i],]-center[i,])^q)^(1/q)
    S[i] = mean(pracma::pdist2(y[label==cluster[i],],center[i,])^q)^(1/q)
  }
  for(i in 1:cl){
    for(j in 1:cl){
      M[i,j] = sum(abs(center[i,]-center[j,])^p)^(1/p)
    }
  }
  for(i in 1:cl){
    for(j in 1:cl){
      R[i,j] = (S[i]+S[j])/M[i,j]
    }
  }
  diag(R) = 0
  Rmean = mean(apply(R,1,max))
  return(Rmean)
}
