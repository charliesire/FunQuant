#' @title Associates a fold to each element of a dataset
#'
#' @param size The size of the dataset
#' @param k  The number of folds

#'
#' @return A vector with the fold number associated to each integer in 1:size.
#' @export
#' @examples
#' kfold(16,3)


kfold = function(size, k){
  res = rep(1:k, each = size%/%k)
  if(size%%k > 0){res = c(res, 1:(size%%k))}
  res = sample(res)
  return(res)
}
