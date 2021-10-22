compute.dim <- function(D,nscales,smallest_scale){

  factor <- (1/smallest_scale)^(1/nscales)
#  val <- sort(D)
  r <- rep(0,nscales)
  n = nrow(D)
#  for ( i in 1:n ){
#    dist <- sort(D[,i])
#    D5 <- rbind(D5,as.matrix(dist[2:10]))
#  }
  D5 = t(apply(D,1,DescTools::Small,10))
  D5 = D5[,2:10]

  r[1] <- max(D5)
  r1 <- median(D5)
  percens = vector()
  for(i in 2:nscales){
    percens[i-1] = factor^i/factor^nscales
  }
  percens = quantile(D,percens,type=5)
  for ( i in 2:nscales){
    r[i] <- percens[i-1]
    r[i] <- max(r[i-1]*1.1,r[i])
  }

  dim <- rep(0,nscales)

  corrdim <- corr.dim(D,r1,r[1])
  dim[1] <- corrdim$dim
  c <- corrdim$c
  s <- rep(0,nscales)

  for( i in 1:n ){
    dist = D[,i]
    for( j in 1:nscales ){
      s[j] <- s[j] + length(which(dist<r[j]))
    }
  }
  Cr <- vector()
  for( j in 1:nscales ){
    Cr[j] <- (1 / (n * (n - 1))) * s[j]
  }
  rho <- vector()
  for( j in 1:(nscales-1) ){
    rho[j] <- Cr[j+1]-Cr[j]
  }

  for( j in 1:(nscales-1) ){
    a <- rho[j]
    b <- r[j]
    c <- max(c,1e-6)
    func <- function(x){a-c*x*b^(x-1)}
    dim[j+1] <- pracma::fzero(func,1)$x
  }
  return(list(r=r,dim=dim))
}

corr.dim <- function(D,r1,r2){
  n <- nrow(D)
  s1<-0;s2<-0
  for ( i in 1:n ){
    dist <- D[,i]
    s1 = s1 + length(which(dist<r1))
    s2 = s2 + length(which(dist<r2))
  }
  Cr1 <- (1 / (n * (n - 1))) * s1
  Cr2 <- (1 / (n * (n - 1))) * s2
  no_dims = (log(Cr2) - log(Cr1)) / (log(r2) - log(r1))
  temp1 = rbind(r1^no_dims, r2^no_dims)
  temp2 = rbind(Cr1,Cr2)
  c = pracma::pinv(temp1)%*%temp2
  return(list(dim=no_dims,c=c))
}
