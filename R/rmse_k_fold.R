#' @title Computation of the Root Mean Square Error at each pixel of the outputs by kfold Cross Validation, for different values of ncoeff and npc
#'
#' @param outputs The output samples on which the metamodel will be trained and tested by kfold cross validation
#' @param nb_folds Number of folds
#' @param ncoeff_vec A vector providing the different values of ncoeff to be tested. ncoeff fixes the number of coefficients used for PCA.
#' @param npc_vec A vector providing the different numbers of principal components to be tested.
#' @param return_pred A boolean indicating whether the predicted outputs should be returned or not
#' @param outputs_pred A list of the predicted outputs already obtained with the same parameters. Default is NULL.
#' @param only_positive A boolean indicating whether the predicted outputs should only contained positive values or not. Default is FALSE.
#' @param design a data frame representing the design of experiments.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param wf name of the wavelet filter to use in the decomposition
#' @param boundary a character string which specifies how boundaries are treated. Only "periodic" is currently implemented .
#' @param J depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2). Default is 1.
#' @param kernel Character defining the covariance model: "exp", "gauss", "matern3_2", "matern5_2".
#' @param regmodel Universal Kriging linear trend: "constant", "linear", "interactive".
#' @param normalize Logical. If TRUE both the input matrix X and the response y in normalized to take values in the interval [0, 1].
#' @param optim Character giving the Optimization method used to fit hyper-parameters. Possible values are: "BFGS", "Newton" and "none", the later simply keeping the values given in parameters. The method "BFGS" uses the gradient of the objective. The method "Newton" uses both the gradient and the Hessian of the objective.
#' @param objective  Character giving the objective function to optimize. Possible values are: "LL" for the Log-Likelihood, "LOO" for the Leave-One-Out sum of squares and "LMP" for the Log-Marginal Posterior.
#' @param parameters Initial values for the hyper-parameters. When provided this must be named list with elements "sigma2" and "theta" containing the initial value(s) for the variance and for the range parameters. If theta is a matrix with more than one row, each row is used as a starting point for optimization.
#' @param seed An optional random seed
#' @param noise Boolean specifying whether to execute NoiseKriging or not
#' @param nugget Boolean specifying whether to execute NuggetKriging or not
#'
#' @return A list containing several outputs :
#' - grid_cv is the grid made with the pairs (npc, ncoeff) that are tested
#' - outputs_rmse is a list of objects that have the same dimension as an output, obtained for each pair (npc, ncoeff). Each element (called pixel here) of the objects is the RMSE computed between the predicted values of the pixel and the true value of the pixel.
#' - outputs_pred is an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.
#' @export
#' @import waveslim
#' @import foreach
#' @import GpOutput2D
#' @import rlibkriging
#' @import abind
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2
#' -10*X[i,])**2)/(60*X[i,]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 20))
#' outputs = func2D(design)

#' list_rmse_k_fold = rmse_k_fold(outputs = outputs, nb_folds = 5,
#' ncoeff_vec = c(50,100,200,400), npc_vec = 2:4, design = design)

rmse_k_fold = function(outputs, nb_folds, ncoeff_vec,npc_vec, outputs_pred = NULL, return_pred = FALSE, only_positive = FALSE,design, kernel="matern5_2", wf = "d4", boundary = "periodic",J=1,
                       regmodel = "constant", normalize = FALSE, optim = "BFGS", objective = "LL", parameters = NULL,noise=FALSE, nugget = FALSE,seed = NULL){
  if(is.null(seed)==FALSE){set.seed(seed)}
  bool_outputs_pred = outputs_pred
  length_dim = length(dim(outputs))
  grid_cv = expand.grid(ncoeff = ncoeff_vec, npc = npc_vec)
  folds = kfold(dim(outputs)[length(dim(outputs))], nb_folds)
  outputs_rmse = list()
  relative_error_df = data.frame()
  if(is.null(bool_outputs_pred)){
    outputs_pred_list = list()
    outputs_pred = list()
    for(k in 1:nb_folds){
      outputs_pred_list[[k]] = rmse_training_test(outputs_train = asub(outputs,dims = length(dim(outputs)), idx = which(folds != k)),outputs_test = asub(outputs,dims = length(dim(outputs)), idx = which(folds == k)), ncoeff_vec = ncoeff_vec,npc_vec = npc_vec, return_pred = TRUE, only_positive = only_positive,design_train = as.data.frame(design[folds !=k,]), design_test = as.data.frame(design[folds == k,]),
                                                  kernel=kernel, wf = wf, boundary = boundary,J=J,
                                                  regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters,noise=noise, nugget = nugget)$outputs_pred
    }
  }
  for(i in 1:nrow(grid_cv)){
    if(is.null(bool_outputs_pred)){
      outputs_pred[[i]] = array(0, dim = dim(outputs))
      dimnames(outputs_pred[[i]]) = lapply(dim(outputs_pred[[i]]), function(i){1:i})
      for(k in 1:nb_folds){
        dimnames(outputs_pred_list[[k]][[i]]) = c(lapply(dim(outputs_pred_list[[k]][[i]])[-length(dim(outputs_pred_list[[k]][[i]]))], function(i){1:i}), list(which(folds == k)))
        afill(outputs_pred[[i]]) = outputs_pred_list[[k]][[i]]
      }
    }
    err = (outputs_pred[[i]] - outputs)^2
    outputs_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
  }
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(grid_cv = grid_cv, outputs_rmse = outputs_rmse, outputs_pred = outputs_pred))
}


