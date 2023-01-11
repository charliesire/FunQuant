
#' @title Provide the Voronoï cell number associated to each sample
#'
#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma A set of prototypes. Useful only if cell_numbers == NULL.
#' @param distance_func A function computng a distance between two elements in the output spaces. Useful only if cell_numbers == NULL.
#'
#' @return A vector providing the Voronoï cell number associated to each sample.
#' @export
#'
#' @examples
get_cell_numbers = function(outputs, gamma, distance_func){
  drop = "selected"
  if(sum(dim(gamma[[1]]) == 1) == length(dim(gamma[[1]]))){drop = F}
  return(Vectorize(function(it){distance_to_gamma(x = Subset(x = outputs, along = length(dim(outputs)), indices = it,drop = drop), gamma = gamma, distance_func = distance_func)$cellule})(1:dim(outputs)[length(dim(outputs))]))
}



