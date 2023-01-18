
#' @title Provide the Voronoï cell number associated to each sample
#'
#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma A set of prototypes. Useful only if cell_numbers == NULL.
#' @param distance_func A function computng a distance between two elements in the output spaces. Useful only if cell_numbers == NULL.
#'
#' @return A vector providing the Voronoï cell number associated to each sample.
#' @import ClimProjDiags
#' @export
#'
#' @examples
#' gamma = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3)))
#' outputs = array(runif(9*20)*20, dim = c(3,3,20))
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' get_cell_numbers(outputs = outputs, gamma = gamma, distance_func = distance_func)
get_cell_numbers = function(outputs, gamma, distance_func){
  drop = "selected"
  if(sum(dim(gamma[[1]]) == 1) == length(dim(gamma[[1]]))){drop = F}
  return(Vectorize(function(it){distance_to_gamma(x = Subset(x = outputs, along = length(dim(outputs)), indices = it,drop = drop), gamma = gamma, distance_func = distance_func)$cellule})(1:dim(outputs)[length(dim(outputs))]))
}



