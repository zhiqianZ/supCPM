#' Run Supervised Capacity Preserved Mapping on Big Data
#'
#' Run supCPM visualization on the data matrix. supCPM is a clustering guided visualization method designed for single-cell RNA data.
#' supCPM requires both high dimensional data and label information as input. Labels could be obtained by clustering algorithm or experiments.
#' With the adjustment of parameters, supCPM could also output the result of CPM (Capacity Preserving Mapping).
#'
#'
#' @param data matrix; Input data matrix with rows representing cells and columns representing genes.
#' @param label vector; The label vector indicating which cluster that each cell belongs to.
#' @param alpha numeric; The quantile of distance pairs, which is used to determine the threshold for downsampling.
#' @param ratio numeric between 0 and 1; This is a trade-off between the KL-divergence and trace ratio term.
#' @param no_dims integer; The dimensionality of the embedding space (default is 2).
#' @param compel_force 0 or 1; Compel_force equals to 0 if user wants to pull clusters a bit apart, and equals to 1 if user wants the best preservation of the geometry (default is 0).
#' @param dist 'euclidean' or 'geodesic'; Choice of the distance used in the high dimensions (default is 'euclidean').
#' @param degree numeric; The degree of freedom of the t-distribution in the high dimensions (default 2).
#' @param k integer; This controls the size of nearest neighbors when calculating the geodesic distance (default is 7).
#' @param niter1 integer; Then iteration number of the first phase (default is 500).
#' @param niter2 integer; The iteration of the second phase (default is 700).
#' @param seed integer; Random seed for supCPM. (default is 40)
#' @param factor numeric; The factor multiples on the high-dimensional distance matrix between different clusters.
#' @param verbose logical; Whether the progress of the objective function should be printed (default is TRUE).
#' @param init logical; Whether to ues the random initialization or not. If the init=T,
#'              then the initialization will be the MDS of the adjusted capacity distance matrix. (default is FALSE)
#' @param epsilon numeric; The multiple of the default epsilon, which is used in converting distance to probability; Default value is 1.
#' @param lr numeric; The learning rate for optimization. Default value is 500.
#'
#' @return Returns an matrix containing supCPM coordinates in the embedding dimenisons.
#' @references Zhiqian Zhai, Yu L. Lei, Rongrong Wang, Yuying Xie, Supervised Capacity Preserving Mapping: A Clustering Guided Visualization Method for scRNAseq data,
#' bioRxiv 2021.06.18.448900; doi: https://doi.org/10.1101/2021.06.18.448900
#' @export
supCPM_downsample = function(data,label,alpha =0.02,no_dims=2,compel_force=0,dist='euclidean',degree=2,ratio,k=7,niter1=700,niter2=700,seed=42,factor=1.3,verbose=T,init=T){
  n = nrow(data)
  D = sqDist(data)
  I = order(D[,])
  base = ceiling(length(I)*alpha)
  base_id = arrayInd(base,c(n,n))
  threshold = D[ I[base] ]
  merge = rep(0,n)
  indtable = matrix(0,n,2)
  cl_size =  table(label)
  large_cluster = which(cl_size>n*0.05)
  large_cluster = sort(unique(label))[large_cluster]
  k=1
  for (i in large_cluster){
    cl = which(label==i)
    for( j in 1:length(cl) ){
      if( merge[cl[j]] == 0){ # whether the point is reomved
        smallind = which(D[cl[j],cl]< threshold )
        smallind = intersect(smallind, which(merge[cl]==0))
        if ( length(smallind) > 1){
          merge[cl[smallind]] = k
          indtable[k,] = c(cl[j],k)
          k=k+1
        }
      }
    }
  }
  indtable = indtable[-which(indtable[,2]==0),]
  ind = indtable[,1]; ind = c(ind,which(merge==0))
  data_downsample = data[ind,]
  label_downsample = label[ind]
  print(dim(data_downsample))
  embed = supCPM(data_downsample,label_downsample,no_dims,compel_force,dist,degree,ratio,k,niter1,niter2,seed,factor,verbose,init)
  embed_recover = matrix(0,n,no_dims)
  embed_recover[ind,] = embed
  rsmall = dist(embed_recover[base_id,])
  D = sqDist(embed)
  diag(D)=0
  rsmall = sqrt(min(D[D!=0]))
  for ( i in 1:nrow(indtable)){
    i = 1
    cl = which(merge==i)
    center = indtable[i,1]
    random = matrix(rnorm(length(cl)*no_dims),ncol=no_dims)
    #id = which(abs(random)>1)
    #random[id] = sign(random[id])
    embed_recover[cl,] = matrix(rep(embed_recover[center,],length(cl)),ncol=no_dims,byrow = T) + random * sqrt(rsmall/2)
  }
  return(embed_recover)
}








































