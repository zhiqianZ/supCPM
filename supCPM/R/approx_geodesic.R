approx.geo <- function(W,G){
 n <- length(igraph::V(G))
 select <- 50
 randind = sample(n)
 test_set = tail(randind,11)
 test_col = t(igraph::shortest.paths(G,test_set))
 sample_col = vector()
 for( draw in 1:floor(n/select) ){
   sample_col <- cbind(sample_col,t(igraph::shortest.paths(G,randind[(1+(draw-1)*select) : (select*draw)])))
   Y <- nystrom(sample_col^2,randind[1:(select*draw)],test_set)^0.5
   if( norm(Y-test_col,'2')<3e-2*norm(test_col,'2') ){
     break
   }
 }

 D = nystrom(sample_col^2,randind[1:(select*draw)],1:n)^0.5
 D[W!=0]=W[W!=0]

 return(D)
}


nystrom <- function(X,set,test_set){
  n <- length(set)
  svd_sys <- svd(t(X[set,]),nu=nrow(t(X[set,])),nv=ncol(t(X[set,])))
  U <- svd_sys$u
  V <- svd_sys$v
  S <- matrix(0,nrow(t(X[set,])),ncol(t(X[set,])))
  diag(S) <- svd_sys$d
  n0 <- round(n*0.5)
  lambda <- S[n0,n0]
  S[S<=lambda] <- 0
  Y = X %*% pracma::pinv(V%*%S%*%t(U)) %*% t(X[test_set,])
  Y[Y<0] <- 0
  return(Y)
}
