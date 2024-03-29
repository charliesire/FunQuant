#' @title Compute the weights fX/g of a sample
#'
#' @param f The real density function of the inputs
#' @param g The biased density function that helped sampling the inputs
#' @param inputs The value of the sampled inputs
#'
#' @return A vector with the weights fX/g of the inputs
#' @import foreach
#' @export
#'
#' @examples
#' g = function(x){prod(sapply(x, function(y){dnorm(y)}))}
#' f = function(x){1}
#' inputs = array(rnorm(30), dim = c(10,3))
#' compute_density_ratio(f,g, inputs)
compute_density_ratio = function(f, g, inputs){
  res = Vectorize(function(it){
    return(as.numeric(f(inputs[it,])/g(inputs[it,])))
    })(1:nrow(inputs))
  return(res)
}


