#' @title Compute the estimator which is the numerator of the centroid estimation
#'
#' @param outputs The output samples that need to be quantized
#' @param cell_numbers The voronoi cell number of every output
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param cell The cell number of the computed centroid
#'
#' @return An array having the same dimension as an output, which is the numerator of the centroid estimator
#' @export
#'
#' @examples
estim_1 = function(outputs, cell_numbers, density_ratio, cell){
  outputs_cell = Subset(x = outputs, along = length(dim(outputs)), indices = which(cell_numbers == cell))
  outputs_cell = matrix(outputs_cell, nrow = prod(dim(outputs)[1:(length(dim(outputs))-1)]))
  outputs_cell = t(outputs_cell)
  res = outputs_cell*density_ratio[cell_numbers == cell]
  res = apply(res,2,sum)
  res = array(res, dim = dim(outputs)[1:(length(dim(outputs))-1)])
  return(res)
}




