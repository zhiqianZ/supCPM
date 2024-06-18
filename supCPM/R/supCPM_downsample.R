#' Run Supervised Capacity Preserved Mapping on Big Data
#'
#' Run supCPM visualization on the data matrix. supCPM is a clustering guided visualization method designed for single-cell RNA data.
#' supCPM requires both high dimensional data and label information as input. Labels could be obtained by clustering algorithm or experiments.
#' With the adjustment of parameters, supCPM could also output the result of CPM (Capacity Preserving Mapping).
#'
#' supCPM downsample use Geometric sketching to downsample a subset of cells to run the supCPM algoirthm, and the remaining unselected cell's visualization coordinates will be projected the center of its kNN among the downsampled cells plus a random perturbation.
#'
#' Prior to run this function, you need to install reticulate and python package geoskech first.
#'
#' @param data matrix; Input data matrix with rows representing cells and columns representing genes.
#' @param label vector; The label vector indicating which cluster that each cell belongs to.
#' @param alpha numeric; The proportion of cells to downsample (default is 0.1).
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
#' @param init logical; Whether to use the random initialization or not. If the init=T,
#'              then the initialization will be the MDS of the adjusted capacity distance matrix. (default is FALSE)
#' @param return.idx logical; Whether to return the indices of cells that are downsampled to actually run the supCPM algorithm (default is FALSE).
#' @param epsilon numeric; The multiple of the default epsilon, which is used in converting distance to probability; Default value is 1.
#' @param lr numeric; The learning rate for optimization. Default value is 500.
#' @return Returns an matrix containing supCPM coordinates in the embedding dimenisons.
#' @references Zhiqian Zhai, Yu L. Lei, Rongrong Wang, Yuying Xie, Supervised Capacity Preserving Mapping: A Clustering Guided Visualization Method for scRNAseq data,
#' bioRxiv 2021.06.18.448900; doi: https://doi.org/10.1101/2021.06.18.448900
#' @references geometric sketching algorithm described by Brian Hie, Hyunghoon Cho, Benjamin DeMeo, Bryan Bryson, and Bonnie Berger in "Geometric sketching compactly summarizes the single-cell transcriptomic landscape", Cell Systems (2019).
#' @export
#' @import RcppAnnoy
#' @import reticulate
supCPM_downsample = function(data,label,alpha =0.1,
                             no_dims=2,compel_force=0, dist='euclidean',
                             degree=2,ratio,k=7,niter1=700,niter2=700,seed=42,factor=1.3,epsilon=1,lr=500,
                             verbose=T,init=T,return.idx=F){
  library(RcppAnnoy)
  n = length(label)
  set.seed(seed)
  geosketch <-  reticulate::import('geosketch')
  sketch.size <- as.integer(n*alpha)
  sketch.indices <- geosketch$gs(data, sketch.size)
  sketch.indices = unlist(sketch.indices)
  if(sum(table(label[sketch.indices])<10)>0){
    #if there is any cell types contains fewer than 10 cells after downsampling, use all
    celltype = names(which(table(label[sketch.indices]))<10)
    sketch.indices = union(sketch.indices,which(label%in%celltype))
  }
  data_downsample = data[sketch.indices,]
  label_downsample = as.vector(label[sketch.indices])
  embed = supCPM(data_downsample,label_downsample,no_dims=no_dims,compel_force=compel_force,
                 dist=dist,degree=degree,ratio=ratio,k=k,niter1=niter1,niter2=niter2,
                 seed=seed,factor=factor,verbose=verbose,init=init,
                 epsilon=epsilon,lr=lr
                 )
  embed_recover = matrix(0,n,no_dims)
  embed_recover[sketch.indices,] = embed
  select = rep(0, length(label))
  select[sketch.indices] = "sketched cell"
  select[-sketch.indices] = "projected cell"
  for( i in unique(label) ){
    idx = which(label==i)
    idx_subsampled = idx[idx%in%sketch.indices]
    idx_new = idx[!idx%in%sketch.indices]
    d = as.matrix(distances::distances(embed[label_downsample==i,]))
    min_dist = min(d[d!=0])
    Annoy <- new(AnnoyEuclidean, ncol(data[idx_subsampled,]))
    Annoy$setSeed(seed)
    for (i in 1:length(idx_subsampled)) Annoy$addItem(i-1, data[idx_subsampled[i],])
    Annoy$build(50)
    neighbor.size = min(10, length(idx_subsampled))
    NN = t(sapply(1:length(idx_new), function(i) idx_subsampled[Annoy$getNNsByVector(as.vector(data[idx_new[i],]), neighbor.size+1)+1]))
    impute = t(apply(NN,1,function(idx) colMeans(embed_recover[idx,])))
    embed_recover[idx_new, ] = impute + matrix(rnorm(nrow(impute)*ncol(impute),0, mean = sqrt(min_dist/no_dims)),nrow=nrow(impute))
  }
  if(return.idx==T){
    return(list(supCPM=embed_recover,idx=select))
  }
  return(embed_recover)
}








































