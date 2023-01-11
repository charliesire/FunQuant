#' @title Compututation of the empirical quantization error
#'

#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma A set of prototypes. Useful only if cell_numbers == NULL.
#' @param density_ratio density_ratio indicates the weight fX/g of each output

#' @return An estimation of the quantization error
#' @export
#'
#' @examples
quanti_error = function(outputs, gamma, density_ratio){
  drop = "selected"
  if(sum(dim(gamma[[1]]) == 1) == length(dim(gamma[[1]]))){drop = F}
  distances = Vectorize(function(it){distance_to_gamma(x = Subset(x = outputs, along = length(dim(outputs)), indices = it,drop = drop), gamma = gamma, distance_func = distance_func)$dist})(1:dim(outputs)[length(dim(outputs))])
  return(sqrt(mean(distances^2*density_ratio)))
}
