#' @title Computation of the Root Mean Square Error at each pixel of the outputs by leave one out, for different values of ncoeff and npc
#'
#' @param outputs The training output samples on which the metamodel will be trained
#' @param model_tuning An optional list of models created for each ncoeff values.
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
#' @param noise Boolean specifying whether to execute NoiseKriging or not
#' @param nugget Boolean specifying whether to execute NuggetKriging or not
#' @param seed An optional random seed
#' @return A list containing several outputs :
#' - grid_cv is the grid made with the pairs (npc, ncoeff) that are tested
#' - outputs_rmse is a list of objects that have the same dimension as an output, obtained for each pair (npc, ncoeff). Each element (called pixel here) of the objects is the RMSE computed between the predicted values of the pixel and the true value of the pixel.
#' - outputs_pred is an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.

#' @export
#' @import waveslim
#' @import foreach
#' @rawNamespace import(GpOutput2D, except = Fpca2d.Wavelets)
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
#' design = data.frame(X = seq(-1,1,l= 8))
#' outputs = func2D(design)

#' list_rmse_loo = rmse_loo(outputs = outputs, ncoeff_vec =
#' c(50,100,200,400), npc_vec = 2:4, design = design)
rmse_loo = function(outputs, model_tuning = NULL, ncoeff_vec,npc_vec, return_pred = FALSE, outputs_pred = NULL, only_positive = FALSE,design,
                    kernel="matern5_2", wf = "d4", boundary = "periodic",J=1,
                    regmodel = "constant", normalize = FALSE, optim = "BFGS", objective = "LL", parameters = NULL,noise=FALSE, nugget = FALSE,seed = NULL){

  grid_cv = expand.grid(ncoeff = ncoeff_vec, npc = npc_vec)
  outputs_rmse = list()
  outputs_loo_list = list()
  if(is.null(outputs_pred)){
    if(is.null(model_tuning)){model_tuning = create_models_tuning(outputs = outputs, ncoeff_vec = ncoeff_vec, npc = max(npc_vec), design = design, kernel=kernel, wf = wf, boundary = boundary,J=J,
                                                                  regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters,noise=noise, nugget = nugget)}
  }
  for(i in 1:nrow(grid_cv)){
    if(is.null(outputs_pred)){
    ncoeff = grid_cv[i,1]
    npc = grid_cv[i,2]
    indice_coeff = which(ncoeff_vec == ncoeff)
    fp = Fpca2d.Wavelets(outputs, wf = wf, boundary = boundary, J = J, ncoeff = ncoeff, rank = npc)
    model = lapply(1:npc, function(k){model_tuning[[indice_coeff]][[k]]})
    pred = sapply(1:npc, function(i){model[[i]]$leaveOneOutVec(theta = model[[i]]$theta())$mean})
    outputs_loo = inverse_Fpca2d(pred,fp)
    if(only_positive){outputs_loo = (outputs_loo > 0)*outputs_loo}
    }
    else{outputs_loo = outputs_pred[[i]]}
    if(return_pred){outputs_loo_list[[i]] = outputs_loo}
    err = (outputs_loo - outputs)^2
    outputs_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
  }
  if(return_pred == FALSE){outputs_loo = NULL}
  return(list(grid_cv = grid_cv, outputs_rmse = outputs_rmse, outputs_pred = outputs_loo_list))
}

