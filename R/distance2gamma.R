#' @title Compute the distance between a point and its nearest centroid, returning this distance and the associated cell number
#'
#' @param x A point in the output space
#' @param gamma A set of prototypes
#' @param distance_func A function computng a distance between two elements in the output spaces
#'
#' @return
#' @export
#'
#' @examples
distance_to_gamma = function(x, gamma, distance_func){
  distance = Vectorize(function(k){distance_func(x, gamma[[k]])})(1:length(gamma))
  return(list(cellule = which.min(distance), dist = min(distance)))
}
