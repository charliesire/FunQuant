#' @title Computation of GP models in the PCA space for different ncoeff values
#'
#' @param outputs The training output samples on which the metamodel will be trained
#' @param ncoeff_vec a vector providing the different values of ncoeff that will be tested. ncoeff fixes the number of coefficients used for PCA.
#' @param npc Maximal number of principal components to be tested
#' @param formula  an object of class "formula"
#' (or a list of "formula" which the length is equal to the number of modeled principal components)
#' specifying the linear trend of the kriging model (see \code{\link{lm}}) on each principal component.
#'  This formula should concern only the input variables (\code{design}), and not the output (\code{response}).
#'  The default is ~1, which defines a constant trend on each principal component.
#' @param design a data frame representing the design of experiments.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param covtype optional character string or vector of character strings
#' specifying the covariance structure to be used on each modeled principal component
#' (see \code{\link{km}} for possible inputs of \code{covtype}).
#' If a vector, the length should be equal to the number of modeled principal components.
#' @param boundary a character string which specifies how boundaries are treated. Only "periodic" is currently implemented (see \code{\link{dwt.2d}}).
#' @param J depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2). Default is 1.
#' @param coef.trend,coef.cov,coef.var optional vectors or matrices containing
#' the values for the trend, covariance and variance parameters.
#' If matrices, the number of rows should be equal to the number of modeled principal components.
#' For details, see \code{\link{km}}).
#' @param nugget an optional variance value or vector standing for the homogeneous nugget effect.
#' If vector, the length should be equal to the number of modeled principal components.
#' @param noise.var an optional vector or matrix containing the noise variance
#' at each observation on each modeled principal component.
#' @param lower,upper optional vectors or matrices containing the bounds of the correlation parameters
#' of each principal component for optimization. For details, see \code{\link{km}}).
#' @param parinit an optional vector or matrix containing the initial values for the variables to be optimized over.
#' For details, see \code{\link{km}}).
#' @param multistart an optional integer indicating the number of initial points from which running the BFGS optimizer.
#'  (see \code{\link{km}}).
#' @param kernel an optional function or list of functions containing a new covariance structure
#' for each principal component. At this stage, the parameters must be provided as well, and are not estimated.
#' @param control an optional list of control parameters for optimization. For details, see \code{\link{km}}).
#' @param ... other parameters of \code{\link{km}} function from \code{DiceKriging}.
#'
#' @return A list of lists of objects of class \code{\linkS4class{km}}. The length of the output list is the length of ncoeff_vec. For each ncoeff of ncoeff_vec, a list of npc objects of class \code{\linkS4class{km}} is computed.
#' @export
#'
#' @examples
create_models_tuning = function(outputs, ncoeff_vec, npc, formula = ~1,design, covtype="matern5_2", boundary = "periodic",J=1,
                                coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                                nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                                parinit = NULL, multistart=1,
                                kernel=NULL, control = NULL,...){
  model_tuning = list()
  for(i in 1:length(ncoeff_vec)){
    print(i)
    set.seed(1)
    ncoeff = ncoeff_vec[i]
    fp = Fpca2d.Wavelets(outputs, wf = "d4", boundary = boundary, J = J, ncoeff = ncoeff, rank = npc)
    model_tuning[[i]] = km_Fpca2d(formula = formula, design = design, response = fp, covtype = covtype, coef.trend = coef.trend, coef.var = coef.var, coef.cov = coef.cov, control = control, nugget = nugget, multistart=multistart, lower = lower)
  }
  return(model_tuning)
}
