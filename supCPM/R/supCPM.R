#' Run Supervised Capacity Preserved Mapping (supCPM)
#'
#' Run supCPM visualization on the data matrix. supCPM is a clustering guided visualization method designed for single-cell RNA data.
#' supCPM requires both high dimensional data and label information as input. Labels could be obtained by clustering algorithm or experiments.
#' With the adjustment of parameters, supCPM could also output the result of CPM (Capacity Preserving Mapping).
#'
#'
#' @param data matrix; Input data matrix with rows representing cells and columns representing genes.
#' @param label vector; The label vector indicating which cluster that each cell belongs to.
#' @param ratio numeric between 0 and 1; This is a trade-off between the KL-divergence and trace ratio term. Default is 0.7.
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
#' @param epsilon numeric; The multiple of the default epsilon, which is used in converting distance to probability; Default value is 1.
#' @param lr numeric; The learning rate for optimization. Default value is 500.
#' @param intermediate logic; Whether to output coordinates of the embedding points before switching objective function or not. Default is False.
#' @return If intermediate=False, returns an matrix containing supCPM coordinates in the embedding dimenisons. If intermediate=True, return a list with the final embedding and the intermediate embedding points.
#' @references Zhiqian Zhai, Yu L. Lei, Rongrong Wang, Yuying Xie, Supervised Capacity Preserving Mapping: A Clustering Guided Visualization Method for scRNAseq data,
#' bioRxiv 2021.06.18.448900; doi: https://doi.org/10.1101/2021.06.18.448900
#' @export
#' @import igraph
#' @importFrom DescTools Small
#' @import pracma
supCPM <- function(data,label,ratio=0.7,no_dims=2,compel_force=1,dist='euclidean',degree=2,k=7,niter1=500,niter2=700,seed=40,factor=1.3,verbose=T,init=T,epsilon=1,lr=500,intermediate=F){
  if(any(is.na(data))){stop('Data contains NAs!')}
  data = as.matrix(data)
  if(any(is.na(data))){stop('Data contains NAs by coercion!')}
  if((no_dims%%1!=0)|no_dims<=0){stop('The embedding dimension should be an integer larger than 0!')}
  if(compel_force!=0&compel_force!=1){stop('Parameter \'compel_force\' should be 0 or 1!')}
  if(degree<=0){stop('Parameter \'degree\' should be larger than 0!')}
  if(is.numeric(dist)){stop('Distance metric used should be indicated by string!')}
  if(tolower(dist)!='euclidean'&tolower(dist)!='geodesic'){stop('Distance could only be either Euclidean or geodesic!')}
  if(ratio>1|ratio<0){stop('Parameter \'ratio\' should be between 0 and 1!')}
  if((k%%1!=0)|k<=0){stop('Parameter \'k\' should be an integer larger than 0!')}
  if((niter1%%1!=0)|niter2<0){stop('The first number of interation should be an integer equal or larger than 0!')}
  if((niter2%%1!=0)|niter2<0){stop('The second number of iteration should be an integer equal or larger than 0!')}
  if(factor<0){stop('Parameter \'factor\' should be larger than 0!')}
  if(epsilon<0){stop('Parameter \'epsilon\' should be larger than 0!')}
  if(lr<0){stop('Learning rate should be larger than 0!')}
  if(lr<0){stop('Learning rate should be larger than 0!')}
  if(!is.logical(verbose)){
    if(verbose!=0&verbose!=1){
      warning('Parameter \'verbose\' should be a logical value! \'verbose\' is coerced into logic!')
      verbose = as.logical(verbose)
    }
    verbose = as.logical(verbose)
  }
  if(!is.logical(init)){
    if(init!=0&init!=1){
      warning('Parameter \'init\' should be a logical value! \'init\' is coerced into logic!')
      init = as.logical(init)
    }
    init = as.logical(init)
  }
  if(!is.vector(label)){
    if(!any(dim(label)==1)){stop('label must be a vector!')}
    if(is.data.frame(label)){label=as.matrix(label)}
  }
  if(any(is.na(label))){stop('labels shouldn\'t contain any NAs!')}
  if(is.character(label[1])){
    label = as.numeric(as.factor(label))
    if(any(is.na(label))){stop('labels contain NAs by coercion!')}}

  # sorting labels
  if(min(label)==0){label=label+1}
  label_sort = sort(label)
  recover_ind = rank(label,ties.method = 'first')
  perm_ind = order(label)
  data_sort = data[perm_ind,]
  cl_size = table(label)
  cl_num = length(cl_size)
  if (dist == 'euclidean'){
    D = compute.dist(data_sort,0)
  }else if (dist == 'geodesic'){
    cat('computing the geodesic distance\n')
    D = compute.dist(data_sort,1,k)
  }else{
    stop('Input dist should be either \'euclidean\' or \'geodesic\' ')
  }
  diag(D) =  0
  D = D/max(D)
  # compute the Capacity adjusted distance
  D1 = cdist(D,no_dims,compel_force,epsilon)
  # could be improved
  # construction of label matrix
  n = nrow(data)
  L_tilde <- matrix(rep(0,n^2),n)
  for (i in label_sort){
    cl = which(label_sort==i)
    L_tilde[cl,cl] = 1/length(cl)
  }
  D1[L_tilde==0] = D1[L_tilde==0] * factor
  n = nrow(D)
  D <- D / max(D)
  D[D<.Machine$double.xmin] <- .Machine$double.xmin
  # change the distances to probabilities
  P_D <- 1/((D1^2)^((1+degree)/2))
  P <- P_D/sum(P_D)
  diag(P) <- 0
  P <- P/sum(P)
  P[P<.Machine$double.xmin] <- .Machine$double.xmin
  P <- 0.5*(P+t(P))
  const <- Re(sum((P*log(P))))

  P <- 4*P

  set.seed(seed)
  ydata <- 0.0001*matrix(rnorm(n*no_dims),nrow=n)
  #initialization
  if(init==T){  ydata = cmdscale(D) }
  momentum <- 0.5
  final_momentum <- 0.8
  mom_switch_iter <- 500
  stop_lying_iter <- 100
  epsilon <- lr
  min_gain <- 0.01
  y_incs  <- matrix(0,nrow(ydata),ncol(ydata))
  gains <- matrix(1,nrow(ydata),ncol(ydata))

  if(niter1>0){
  for (iter in 1:niter1){
     D_ydata = sqDist(ydata)
     num = 1/(1+D_ydata)
     diag(num) = 0
     s = sum(num)
     Q = num/s
     Q[Q<.Machine$double.xmin] = .Machine$double.xmin
     L = (P-Q) * num
     y_grads = 4 * colSums(L) * ydata
     y_grads = y_grads - 4*L%*%ydata
     gains <- (gains + 0.2) * (sign(y_grads) != sign(y_incs))+(gains*0.8)*(sign(y_grads) == sign(y_incs))
     gains[gains < min_gain] <- min_gain
     y_incs <- momentum * y_incs - epsilon * (gains * y_grads)
     ydata <- ydata + y_incs
     ydata <- t(apply(ydata, 1, function(x,y){x+y}, -colMeans(ydata)))
     if (iter == mom_switch_iter){momentum <- final_momentum}
     if (iter == stop_lying_iter){P = P/4}
     if (iter%%10==0){
       cost <- const - sum(P*log(Q))
       cat('Iteration ', iter, ': Objective Function Value ', cost, '\n')
     }
  }
  }
  ydata_temp = ydata
  epsilon <- epsilon/5e2
  #P = P*4
  if(niter2>0){
  for( iter in 1:niter2){
    D_ydata = sqDist(ydata)
    num = 1/(1+D_ydata)
    diag(num) = 0
    s = sum(num)
    Q = num/s
    Q[Q<.Machine$double.xmin] = .Machine$double.xmin
    L = (P-Q) * num
    y_grads = 4 * colSums(L) * ydata
    unsup_grads = y_grads - 4*L%*%ydata

    sup_grads = supGradient(ydata,label_sort,cl_num,cl_size)
    TraceSw = sup_grads$TraceSw
    TraceSb = sup_grads$TraceSb
    sup_grads = sup_grads$sup_grads
    y_grads = (1-ratio) * unsup_grads + ratio * sup_grads
    gains <- (gains + 0.2) * (sign(y_grads) != sign(y_incs))+(gains*0.8)*(sign(y_grads) == sign(y_incs))
    gains[gains < min_gain] <- min_gain
    y_incs <- momentum * y_incs - epsilon * (gains * y_grads)
    ydata <- ydata + y_incs
    ydata <- t(apply(ydata, 1, function(x,y){x+y}, -colMeans(ydata)))
    if (iter%%10==0){
      cost <-  (1-ratio) * const - (1-ratio) * sum(P*log(Q)) + ratio * TraceSb/TraceSw
      cat('Iteration ', iter+niter1, ': Objective Function Value ', cost, '\n')
    }
   # if (iter == stop_lying_iter){P = P/4}
  }
  }
  if(intermediate==T){
    return(list(embedding=ydata[recover_ind,],intermeidate=ydata_temp[recover_ind,]))}
  if(intermediate==F){
    return(embedding=ydata[recover_ind,])
  }
}






