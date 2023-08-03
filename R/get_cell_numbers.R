
#' @title Provide the Voronoï cell number associated to each sample
#'
#' @param data The data that needs to be quantized. An array of any dimension is expected, the kth element must be selected with data[,..,k]. Useful only if cell_numbers == NULL.
#' @param prototypes A set of prototypes. Useful only if cell_numbers == NULL.
#' @param distance_func A function computing a distance between two data elements. Useful only if cell_numbers == NULL.
#'
#' @return A vector providing the Voronoï cell number associated to each sample.
#' @import abind
#' @export
#'
#' @examples
#' prototypes = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3)))
#' data = array(runif(9*20)*20, dim = c(3,3,20))
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' get_cell_numbers(data = data, prototypes = prototypes, distance_func = distance_func)
get_cell_numbers = function(data, prototypes, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}){
  drop = "selected"
  if(sum(dim(prototypes[[1]]) == 1) == length(dim(prototypes[[1]]))){drop = FALSE}
  return(Vectorize(function(it){distance_to_prototypes(x = asub(x = data, dims = length(dim(data)), idx = it,drop = drop), prototypes = prototypes, distance_func = distance_func)$cellule})(1:dim(data)[length(dim(data))]))
}



