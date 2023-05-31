#' @title Compute the centroid and the probability mass of the Voronoï cells
#'
#' @param data The data that needs to be quantized. If method = "percell", a list of data samples must be provided, of length equal to the number of Voronoï cells.
#' @param cell_numbers The voronoi cell number of every data element
#' @param method_IS The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell.
#' @param sampling_cells If method == "percell" and data is provided, sampling_cells is a vector indicating for each cell, the index element of data associated to this cell.
#' @param density_ratio A vector indicating the weight fX/g of each data element. Default is a vector of 1. If method = "percell", a list of density_ratio must be provided, of length equal to the number of Voronoï cells.
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoi cell. Default is 0 for all Voronoi cells.
#' @param batch A boolean indicating whether the computations must be performed by batch or not. If TRUE, data, cell_numbers and density_ratio must be lists. Default is False.

#'
#' @return The centroid and the probability mass of each probability cell
#' @import abind

compute_centroids_and_proba = function(data, cell_numbers, method_IS = "unique", sampling_cells = NULL, density_ratio = rep(1, dim(data)[length(dim(data))]), bias = rep(0,length(unique(unlist(cell_numbers)))), batch = FALSE){
  if(method_IS == "unique"){n = length(unlist(cell_numbers))}
  else if(method_IS == "percell"){n = length(unlist(cell_numbers))/length(cell_numbers)}
  nb_cells = length(unique(unlist(cell_numbers)))
  centroids = list()
  probas = c()
  for(j in 1:nb_cells){
    if(method_IS == "unique"){
      data_j = data
      cell_numbers_j = cell_numbers
      density_j = density_ratio
    }
    else if(method_IS == "percell"){
      data_j = data[[sampling_cells[j]]]
      cell_numbers_j = cell_numbers[[sampling_cells[j]]]
      density_j = density_ratio[[sampling_cells[j]]]
    }
    numerator =  estim_num_centroid(data = data_j, cell_numbers = cell_numbers_j, density_ratio = density_j, cell = j, batch = batch) ## Sum the Y(X)f/nu of the cell
    if(batch){denominator = sum(Vectorize(function(p){estim_denom_centroid(density_ratio = density_j[[p]], cell_numbers = cell_numbers_j[[p]], cell = j, bias = bias[j])})(1:length(density_j)))}
    else{denominator = estim_denom_centroid(density_ratio = density_j, cell_numbers = cell_numbers_j, cell = j, bias = bias[j])}
    centroids[[j]] = numerator/denominator
    probas = c(probas, denominator/n)
  }
  return(list(centroids = centroids, probas = probas))
}
