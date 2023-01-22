#' @title Compututation of the empirical quantization error
#'

#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma A set of prototypes. Useful only if cell_numbers == NULL.
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param distance_func A function computng a distance between two elements in the output spaces. Useful only if cell_numbers == NULL.

#' @return An estimation of the quantization error
#' @export
#'
#' @examples
#' gamma = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3)))
#' outputs = array(runif(9*20)*20, dim = c(3,3,20))
#' density_ratio = rep(1,20)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' quanti_error(outputs = outputs, gamma = gamma, density_ratio = density_ratio, distance_func = distance_func)
quanti_error = function(outputs, gamma, density_ratio, distance_func){
  drop = "selected"
  if(sum(dim(gamma[[1]]) == 1) == length(dim(gamma[[1]]))){drop = F}
  distances = Vectorize(function(it){distance_to_gamma(x = asub(x = outputs, dims = length(dim(outputs)), idx = it,drop = drop), gamma = gamma, distance_func = distance_func)$dist})(1:dim(outputs)[length(dim(outputs))])
  return(sqrt(mean(distances^2*density_ratio)))
}
