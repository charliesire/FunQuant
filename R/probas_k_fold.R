#' @title Computing the relatives errors when predicting the membership probabilities kfold cross validation, for different values of ncoeff and npc
#'
#' @param outputs The output samples on which the metamodel will be trained and tested by kfold cross validation
#' @param nb_folds Number of folds
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param gamma A set of l prototypes defining the Voronoï cells
#' @param distance_func  A function computng a distance between two elements in the output spaces.
#' @param ncoeff_vec A vector providing the different values of ncoeff to be tested. ncoeff fixes the number of coefficients used for PCA.
#' @param npc_vec A vector providing the different numbers of principal components to be tested.
#' @param return_pred A boolean indicating whether the predicted outputs should be returned or not
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
#' @param type A character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK"). Default is "UK.
#' @param seed An optional random seed
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoï cell.
#' @param ... other parameters of \code{\link{km}} function from \code{DiceKriging}.
#'
#'
#' @return A list containing several outputs :
#' - probas_pred_df a dataframe indicating for each pair (npc, ncoeff) the obtained predicted membership probabilities
#' - relative_error_df a dataframe indicating for each pair (npc, ncoeff) the relative error when predicting the membership probabilities
#' - outputs_pred an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.
#' @export
#'
#' @examples
probas_k_fold = function(outputs, nb_folds, density_ratio, gamma, distance_func, ncoeff_vec,npc_vec, return_pred = FALSE,formula = ~1,design, covtype="matern5_2",boundary = "periodic",J=1,
                         coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                         nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                         parinit = NULL, multistart=1,
                         kernel=NULL,control = NULL,type = "UK",seed = NULL,bias = NULL,...){
  if(is.null(seed) == FALSE){set.seed(seed)}
  length_dim = length(dim(outputs))
  dimnames(outputs) = NULL
  grid_cv = expand.grid(ncoeff_vec, npc_vec)
  folds = kfold(length(density_ratio), nb_folds)
  probas_true = get_probas(density_ratio = density_ratio, outputs = outputs, gamma = gamma, distance_func = distance_func, cells = 1:length(gamma), bias = bias)
  probas_pred_df = data.frame()
  outputs_pred_list = list()
  outputs_pred = list()
  relative_error_df = data.frame()
  for(k in 1:nb_folds){
    outputs_pred_list[[k]] = probas_training_test(outputs_train = Subset(outputs,along = length(dim(outputs)), indices = which(folds != k)),outputs_test = Subset(outputs,along = length(dim(outputs)), indices = which(folds == k)), density_ratio = density_ratio, gamma = gamma, distance_func = distance_func, ncoeff_vec = ncoeff_vec,npc_vec = npc_vec, return_pred = TRUE,formula = formula,design_train = design[folds !=k,], design_test = design[folds == k,], covtype=covtype,boundary = boundary,J=J,
                                                  coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                                                  nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                                                  parinit = parinit, multistart=multistart,
                                                  kernel=kernel,control = control,type = type,bias = bias,...)$outputs_pred
  }
  for(i in 1:nrow(grid_cv)){
    outputs_pred[[i]] = array(NA, dim = dim(outputs))
    for(k in 1:nb_folds){
      outputs_pred[[i]][,,folds == k] = outputs_pred_list[[k]][[i]]
    }
    probas_pred_cv = get_probas(density_ratio = density_ratio, outputs = outputs_pred[[i]], gamma = gamma, distance_func = distance_func, cells = 1:length(gamma), bias = bias)
    probas_pred_df = rbind(probas_pred_df, c(as.numeric(grid_cv[i,]), probas_pred_cv))
    relative_error_df = rbind(relative_error_df, c(as.numeric(grid_cv[i,]), abs(probas_pred_cv - probas_true)/probas_true))
  }
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(probas_pred = probas_pred_df, error = relative_error_df, outputs_pred = outputs_pred))
}


