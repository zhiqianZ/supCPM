compute.dist <- function(X,geodesic,k){
  
  n <- nrow(X)
  p <- ncol(X)
  initial_dims <- min(500,p)
  X <- X - min(X)
  X <- X / max(X)
  X <- X - matrix(rep(colMeans(X),n),n,byrow=T)
  if (p<n){
    C <- t(X) %*% X
  }else{
    C <- (1/n)*(X%*%t(X))
  }
  eigen.sys <- eigen(C)
  lambda <- eigen.sys$values
  lambda.ord <- order(lambda,decreasing=T)
  lambda <- sort(lambda, decreasing=T)
  M <- eigen.sys$vectors
  M <- M[,lambda.ord[1:initial_dims]]
  lambda <- lambda[1:initial_dims]
  if (p>=n){
    temp1 <- t(X)%*%M
    temp2 <- 1/sqrt(n*lambda)
    temp2 <- matrix(rep(temp2,nrow(temp1)),nrow(temp1),byrow=T)
    M <- temp1 * temp2
  }
  X <- X %*% M;
  
  # distance matrix
  if (geodesic !=1){
    sum_X <- rowSums(X^2)
    D <- sum_X + t(sum_X-2*X%*%t(X))
    diag(D) <- 0 
    D[D<0] <- 0
  }else{
    D <- geodesic.dis(X,k)
    D <- D^2
  }
  return(D)
}
