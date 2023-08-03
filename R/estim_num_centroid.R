#' @title Compute the estimator which is the numerator of the centroid estimation
#'
#' @param data The data that needs to be quantized. An array of any dimension is expected, the kth element must be selected with data[,..,k].
#' @param cell_numbers The voronoi cell number of every data element
#' @param density_ratio density_ratio indicates the weight fX/g of each data element
#' @param cell The cell number of the computed centroid
#' @param batch A boolean indicating whether the computations must be performed by batch or not. If TRUE, data, cell_numbers and density_ratio must be lists. Default is False.
#' @return An array having the same dimension as an data element, which is the numerator of the centroid estimator
#' @export
#' @import abind
#' @examples
#' data = array(runif(9*20)*15, dim = c(3,3,20))
#' cell_numbers = c(1,3,2,1,2,1,1,2,3,3,2,2,2,2,2,3,1,1,3,3)
#' density_ratio = rep(1,20)
#' cell = 3
#' estim_num_centroid(data = data,cell_numbers = cell_numbers,
#' density_ratio = density_ratio, cell = cell)
estim_num_centroid = function(data, cell_numbers, density_ratio, cell, batch = FALSE){
  if(batch){
    res = 0
    for(batch_i in 1:length(data)){
      data_cell = asub(x = data[[batch_i]], dims = length(dim(data[[batch_i]])), idx = which(cell_numbers[[batch_i]] == cell))
      data_cell = matrix(data_cell, nrow = prod(dim(data[[batch_i]])[1:(length(dim(data[[batch_i]]))-1)]))
      data_cell = t(data_cell)
      data_cell = data_cell*density_ratio[[batch_i]][cell_numbers[[batch_i]] == cell]
      res = res + apply(data_cell,2,sum)
    }
    res = array(res, dim = dim(data[[1]])[1:(length(dim(data[[1]]))-1)])

  }
  else{
  data_cell = asub(x = data, dims = length(dim(data)), idx = which(cell_numbers == cell))
  data_cell = matrix(data_cell, nrow = prod(dim(data)[1:(length(dim(data))-1)]))
  data_cell = t(data_cell)
  data_cell = data_cell*density_ratio[cell_numbers == cell]
  res = apply(data_cell,2,sum)
  res = array(res, dim = dim(data)[1:(length(dim(data))-1)])
  }
  return(res)
}




