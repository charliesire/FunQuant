#' @title Compute the weights fX/g of a sample
#'
#' @param f The real density function of the inputs
#' @param g The biased density function that helped sampling the inputs
#' @param inputs The value of the sampled inputs
#'
#' @return A vector with the weights fX/g of the inputs
#' @export
#'
#' @examples
compute_density_ratio = function(f, g, inputs){
  f_vec = apply(inputs, 1, f)
  g_vec = apply(inputs, 1, g)
  return(f_vec/g_vec)
}

