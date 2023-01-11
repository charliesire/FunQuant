#' @title Compute the centroid and the probability mass of the Voronoï cells
#'
#' @param outputs The output samples that need to be quantized
#' @param cell_numbers The voronoi cell number of every output
#' @param method_IS The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell.
#'
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#'
#' @return The centroid and the probability mass of each probability cell
#' @export
#'
#' @examples
compute_centroids_and_proba = function(outputs, cell_numbers, method_IS = "unique", density_ratio = NULL, bias = rep(0,length(unique(unlist(cell_numbers))))){
  n = length(cell_numbers)#nb of outputs
  nb_cells = length(unique(unlist(cell_numbers)))
  centroids = list()
  probas = c()
  for(j in 1:nb_cells){
    if(method_IS == "unique"){
      outputs_j = outputs
      cell_numbers_j = cell_numbers
      density_j = density_ratio
    }
    else if(method_IS == "percell"){
      outputs_j = outputs[[j]]
      cell_numbers_j = cell_numbers[[j]]
      density_j = density_ratio[[j]]
    }
    numerator =  estim_1(outputs = outputs_j, cell_numbers = cell_numbers_j, density_ratio = density_j, cell = j) ## Sum the Y(X)f/nu of the cell
    denominator = estim_2(density_ratio = density_j, cell_numbers = cell_numbers_j, cell = j, bias = bias[j])
    centroids[[j]] = numerator/denominator
    probas = c(probas, denominator/n)
  }
  return(list(centroids = centroids, probas = probas))
}
