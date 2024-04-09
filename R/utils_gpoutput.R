#' @title Inverse transform of wavelet or orthonormal B-splines basis.
#'
#' @description embeds the coefficients of "Wavelet2D" or "OrthonormalBsplines2D" basis
#' onto two-dimensional domain.
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @param object an object of class \code{\link{Wavelet2D}}.
#'
#' @return a three dimensional array.
#' The first two dimensions correspond to maps dimensions.
#' The third dimension corresponds to the size of the data set.
#'
#' @seealso \code{\link{Inverse2D.Wavelet2D}}
#' @export
#' 
Inverse2D <- function(object){
  UseMethod("Inverse2D",object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Mean proportion of energy
#'
#' @description Wrapper function to compute mean proportion of energy for given coefficents.
#' The valid objects are "Wavelet2D" and "CoefOrthonormalBsplines2D".
#'
#' @param object an object of class \code{\link{Wavelet2D}}
#' which correspond to coefficients from wavelets.
#'
#' @return a matrix K×L which contains the mean proportion of energy of each coefficient,
#'         with K×L is the dimensions of the functional basis.
#'
#'
#'
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @keywords internal
#' @export
#' 
MeanPoe<-function(object){
  UseMethod("MeanPoe",object)
}


