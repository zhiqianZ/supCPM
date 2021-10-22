cdist <- function(D,no_dims,compel_force,epsilon=1){
  n <- nrow(D)
  nscale <- 7
  dim <- compute.dim(sqrt(D),nscale,0.01)
  r <- dim$r
  dim1 <- dim$dim
  if (compel_force == 1){
    dim1[2] <- 1.5 * dim1[2]
  }
  max_dim <- 10*(no_dims-1)
  ind <- (dim1 < max_dim)
  dim <- rep(1,length(dim1)) * max_dim
  dim[ind] <- dim1[ind]
  dim[dim1<no_dims] <- no_dims

  r1 = quantile(sqrt(D),0.01,type=5)*epsilon
  I = order(D)
#  Dst <- matrix(rep(0,n^2),n)
#  for(i in 1:n){
#    for(j in 1:n){
#        curr_dim = find.dim(r,dim,sqrt(D[i,j])) # find the estimated dimension at the scale r
#    }
#  }
  curr_dim = matrix(0,n,n)
  r0 = c(-1,r)
  for (k in 2:(length(r)+1)){
    curr_dim[sqrt(D)<=r0[k]&sqrt(D)>r0[k-1]] = dim[k-1]
  }
  Dst = (r1^2+D)^(curr_dim/no_dims)
  Dst[I] = make.mono(sqrt(Dst[I]),sqrt(D[I]+r1^2))
  Dst = sqrt(Dst)
  return(Dst)
}


find.dim <- function(r,dim,a){
  i = 1
  while (i <= length(r)){
    if (a <= r[i]){
      d <- dim[i]
      break
    }else{
      i <- i+1
    }
  }
  return(d)
}

make.mono <- function(x,y){
  d <- max(x[1]-y[1],0)
  for (i in 2:length(x)){
    x[i] <- max(x[i],x[i-1],y[i]+d)
    d <- max(d,max(x[i]-y[i],0))
  }
  x = x^2
  return(x)
}


