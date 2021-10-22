find.nn <- function(D,k){
  order.mat <- t(apply(D,1,order))
  NN <- matrix(rep(0,nrow(order.mat)^2),nrow(order.mat))
  neighbors <-order.mat[,1:(k+1)]
  for(i in 1:nrow(order.mat)){
    NN[i,order.mat[i,1:(k+1)]] <- D[i,order.mat[i,1:(k+1)]]
    NN[order.mat[i,1:(k+1)],i] <- D[order.mat[i,1:(k+1)],i]
  }
  NN.graph <- igraph::graph_from_adjacency_matrix(NN, mode = c("undirected"), weighted = TRUE)
  return(list(matrix=NN,graph=NN.graph,neighbors=neighbors))
}

geodesic.dis <- function(X,k){
  n <- dim(X)[1]
  # calculate distance matrix
  sum_X <- apply(X^2,1,sum)
  DD <- sum_X + t(sum_X-2*X%*%t(X))
  diag(DD) <- 0
  DD  <- sqrt(DD)
  D <- find.nn(DD,k)
  D_graph <- D$graph
  D <- D$matrix
  blocks <- igraph::components(D_graph)
  if(blocks$no>1){
    for (i in 1:(blocks$no-1)){
      for (j in (i+1): blocks$no){
        c1 <- which(blocks$membership==i)
        c2 <- which(blocks$membership==j)
        Ds <- DD[c1,c2]
        ind <- which(Ds == min(Ds), arr.ind = T)
        a <- ind[1]; b <- ind[2]
        D[c1[a],c2[b]] <- 2*DD[c1[a],c2[b]]
        D[c2[b],c1[a]] <- 2*DD[c2[b],c1[a]]
      }
    }
  }
  D_graph <- igraph::graph_from_adjacency_matrix(D, mode = c("undirected"), weighted = TRUE)
  # shortest path
  if (n < 500){
    D = igraph::shortest.paths(D_graph,algorithm = 'dijkstra')
  }else{
    # if n is too large, use nystrom extension to speed up
    D = approx.geo(D,D_graph)
  }
  return(D)
}
