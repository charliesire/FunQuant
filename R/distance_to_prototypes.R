#' @title Compute the distance between a point and its nearest prototype, returning this distance and the associated cell number
#'
#' @param x A point in the space of the data elements
#' @param prototypes A set of prototypes
#' @param distance_func A function computing a distance between two data elements
#'
#' @return The distance between a point and its nearest centroid
#' @export
#'
#' @examples
#' distance_to_prototypes(array(1:9, dim = c(3,3)), list(array(10, dim = c(3,3)),
#' array(5, dim = c(3,3)), array(6, dim = c(3,3))),
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))})
distance_to_prototypes = function(x, prototypes, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}){
  distance = Vectorize(function(k){distance_func(x, prototypes[[k]])})(1:length(prototypes))
  return(list(cellule = which.min(distance), dist = min(distance)))
}
