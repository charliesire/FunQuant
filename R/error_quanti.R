#' @title Compututation of the empirical quantization error
#'

#' @param data The data that needs to be quantized. Useful only if cell_numbers == NULL.
#' @param prototypes A set of prototypes. Useful only if cell_numbers == NULL.
#' @param density_ratio density_ratio indicates the weight fX/g of each data element.
#' @param distance_func A function computing a distance between two data elements. Useful only if cell_numbers == NULL.
#' @param batch A boolean indicating whether the computations must be performed by batch or not. If TRUE, data, cell_numbers and density_ratio must be lists. Default is False.
#'
#' @return An estimation of the quantization error
#' @export
#' @import abind
#' @examples
#' prototypes = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3)))
#' data = array(runif(9*20)*20, dim = c(3,3,20))
#' density_ratio = rep(1,20)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' quanti_error(data = data, prototypes = prototypes, density_ratio = density_ratio,
#' distance_func = distance_func)
quanti_error = function(data, prototypes, density_ratio, batch = FALSE, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}){
  drop = "selected"
  if(sum(dim(prototypes[[1]]) == 1) == length(dim(prototypes[[1]]))){drop = FALSE}
  if(!batch){
    distances = Vectorize(function(it){distance_to_prototypes(x = asub(x = data, dims = length(dim(data)), idx = it,drop = drop), prototypes = prototypes, distance_func = distance_func)$dist})(1:dim(data)[length(dim(data))])
    res = sqrt(mean(distances^2*density_ratio))
  }
  else{
    distances = as.numeric(sapply(1:length(data), function(b){Vectorize(function(it){distance_to_prototypes(x = asub(x = data[[b]], dims = length(dim(data[[b]])), idx = it,drop = drop), prototypes = prototypes, distance_func = distance_func)$dist})(1:dim(data[[b]])[length(dim(data[[b]]))])}))
    res = sqrt(mean(distances^2*unlist(density_ratio)))
  }
  return(res)
}
