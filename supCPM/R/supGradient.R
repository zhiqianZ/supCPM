supGradient = function(ydata,label,cl_num,cl_size){
  n = dim(ydata)[1]
  p = dim(ydata)[2]
  TraceSw = 0
  end = cumsum(cl_size)
  start = end - cl_size + 1
  for (i in 1:cl_num){
    ypart = ydata[start[i]:end[i],]
    TraceSw = TraceSw + sum(colSums(ypart)^2)/cl_size[i]
  }
  temp = TraceSw
  TraceSw = TraceSw - sum(colSums(ydata)^2)/n

  TraceSb = sum(diag(t(ydata)%*%ydata)) - temp
  y_mean = colMeans(ydata)
  colsum = matrix(0,cl_num,p)
  for (i in 1:p){
    colsum[,i] =  pracma::accumarray(label,ydata[,i])
  }
  for (i in 1:cl_num){
    colsum[i,] = colsum[i,]/cl_size[i];
  }
  temp = colsum[rep(1:nrow(colsum), times = cl_size), ]
  yb = ydata - temp
  yw = t(apply(temp,1,function(x,y){x-y},y_mean))
  grad = (TraceSw * yb- TraceSb * yw)/ TraceSw^2
  return(list(sup_grads = grad,TraceSb = TraceSb,TraceSw =TraceSw))
}
