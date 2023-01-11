#' Title Compute the estimator which is the denominator of the centroid estimation
#'
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param cell_numbers The output samples that need to be quantized
#' @param cell The cell number of the computed centroid
#'
#' @return A real number which is the denominator of the centroid estimation
#' @export
#'
#' @examples
estim_2 = function(density_ratio, cell_numbers, cell, bias){
  return(sum(density_ratio[as.numeric(cell_numbers) == cell])-bias*length(density_ratio))
}
