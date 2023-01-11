#' @title Computing the probability masses of each voronoi cells
#'
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param cell_numbers The voronoi cell number of every output. If NULL, then outputs, gamma and distance_func must be provided.
#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma A set of prototypes. Useful only if cell_numbers == NULL.
#' @param distance_func A function computng a distance between two elements in the output spaces. Useful only if cell_numbers == NULL.
#'
#' @return A vector providing the probability masses of each Voronoï cell.
#' @export
#'
#' @examples
get_probas = function(density_ratio, method_IS = "unique", cell_numbers = NULL, outputs = NULL, gamma = NULL, distance_func, cells = NULL, bias = NULL){
  if(is.null(cell_numbers)){
    if(method_IS == "unique"){
      cell_numbers = get_cell_numbers(outputs, gamma, distance_func)
    }
    else if(method_IS == "percell"){
      cell_numbers = lapply(cells, function(k){get_cell_numbers(outputs[[k]], gamma, distance_func)})
    }
  }

  probas = c()
  for(cell in cells){
    if(method_IS == "unique"){
      probas = c(probas, sum(density_ratio[cell_numbers == cell])/length(density_ratio)-bias[cell])
    }
    else if(method_IS == "percell"){
      probas = c(probas, sum(density_ratio[cell_numbers[[cell]] == cell])/length(density_ratio[[cell]])-bias[cell])
    }

  }

  return(probas)
}
