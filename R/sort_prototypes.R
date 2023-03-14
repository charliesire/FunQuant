#' @title Sorting the prototypes by increasing sum of their elements (absolute value)
#'
#' @param prototypes A set of prototypes
#'
#' @return The same set of prototypes but sorted by increasing sum of their elements (absolute value)
#' @export
#'
#' @examples
#'
#' sort_prototypes(prototypes = list(array(10, dim = c(3,3)),
#' array(5, dim = c(3,3)), array(6, dim = c(3,3))))
sort_prototypes = function(prototypes){
  prototypes_sorted = prototypes
  sums = Vectorize(function(k){sum(prototypes[[k]])})(1:length(prototypes))
  for(k in 1:length(prototypes)){prototypes_sorted[[rank(sums)[k]]] = prototypes[[k]]}
  return(prototypes_sorted)
}
