#' @title Computing the relatives errors when predicting the membership probabilities by leave one out, for different values of ncoeff and npc
#'
#' @param outputs The training output samples on which the metamodel will be trained
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param gamma A set of l prototypes defining the Voronoï cells
#' @param distance_func  A function computing a distance between two elements in the output spaces.
#' @param model_tuning An optional list of models created for each ncoeff values.
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
#' @param wf name of the wavelet filter to use in the decomposition

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
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoï cell.
#' @param ... other parameters of \code{\link{km}} function from \code{DiceKriging}.
#'
#' @return A list containing several outputs :
#' - probas_pred_df a dataframe indicating for each pair (npc, ncoeff) the obtained predicted membership probabilities
#' - relative_error_df a dataframe indicating for each pair (npc, ncoeff) the relative error when predicting the membership probabilities
#' - outputs_pred an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.
#' @export
#' @import waveslim
#' @import foreach
#' @import DiceKriging
#' @import abind
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*
#' Zgrid$z2-10*X[i,])**2)/(60*X[i,]**2))*
#' (Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 20))
#' outputs = func2D(design)
#' gamma = lapply(c(1,3,6,8,10,14,16,18), function(i){outputs[,,i]})
#' density_ratio = rep(1, 20)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}

#' list_probas_loo = probas_loo(outputs = outputs, density_ratio = density_ratio,
#'  gamma = gamma, distance_func = distance_func, ncoeff_vec = c(50,100,200,400),
#'   npc_vec = 2:4, design = design, control = list(trace = FALSE))
probas_loo = function(outputs, density_ratio, gamma, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},model_tuning = NULL, ncoeff_vec,npc_vec,return_pred = FALSE, formula = ~1,design, covtype="matern5_2", wf = "d4", boundary = "periodic",J=1,
                      coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                      nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                      parinit = NULL, multistart=1,
                      kernel=NULL,control = NULL,type = "UK",bias = NULL,...){

  if(is.null(bias)){bias = rep(0,length(gamma))}
  if(is.null(model_tuning)){model_tuning = create_models_tuning(outputs = outputs, ncoeff_vec = ncoeff_vec, npc = max(npc_vec), formula = formula,design = design, covtype=covtype,
                                                                coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                                                                nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                                                                parinit = parinit, multistart=multistart,
                                                                kernel=kernel,control = control,...)}
  grid_cv = expand.grid(ncoeff = ncoeff_vec, npc = npc_vec)
  probas_true = get_probas(density_ratio = density_ratio, outputs = outputs, gamma = gamma, distance_func = distance_func, cells = 1:length(gamma), bias = bias)
  probas_pred_df = data.frame()
  relative_error_df = data.frame()
  outputs_loo_list = list()
  for(i in 1:nrow(grid_cv)){
    ncoeff = grid_cv[i,1]
    npc = grid_cv[i,2]
    indice_coeff = which(ncoeff_vec == ncoeff)
    fp = Fpca2d.Wavelets(outputs, wf = wf, boundary = boundary, J = J, ncoeff = ncoeff, rank = npc)
    model = lapply(1:npc, function(k){model_tuning[[indice_coeff]][[k]]})
    pred = sapply(1:npc, function(i){leaveOneOut.km(model[[i]], type=type)$mean})
    outputs_loo = inverse_Fpca2d(pred,fp)
    if(return_pred == TRUE){outputs_loo_list[[i]] = outputs_loo}
    probas_pred_cv = get_probas(density_ratio = density_ratio, outputs = outputs_loo, gamma = gamma, distance_func = distance_func, cells = 1:length(gamma), bias = bias)
    probas_pred_df = rbind(probas_pred_df, c(as.numeric(grid_cv[i,]), probas_pred_cv))
    relative_error_df = rbind(relative_error_df, c(as.numeric(grid_cv[i,]), abs(probas_pred_cv - probas_true)/probas_true))
  }
  return(list(probas_pred = probas_pred_df, error = relative_error_df, outputs_pred = outputs_loo_list))
}
