sqDist <- function(data){
  D = outer(rowSums(data^2), rowSums(data^2), '+') - tcrossprod(data, 2 * data)
  return(D)
}
