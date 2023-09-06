
#' @title Get the FPCA basis
#'
#' @param fpca an object of class \code{Fpca2d}.
#'
#' @return The FPCA basis maps
#'
#' @export
#'

find_basis = function(fpca){
  scores = diag(ncol(fpca$x))

  return(inverse_Fpca2d(scores, fpca))
}
