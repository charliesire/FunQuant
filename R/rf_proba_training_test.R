#' @title Computing the relative errors when predicting the membership probabilities of a validation output sample with a classification step, for different values of hyperparameters
#'
#' @param design_train a data frame representing the design of experiments of the training part.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param design_test a data frame representing the design of experiments of the validation part.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param outputs_train The training output samples on which the metamodel will be trained
#' @param outputs_test  The validation output samples on which the metamodel performance will be evaluated
#' @param threshold_classification The threshold that creates the two classes of maps for the classification
#' @param threshold_fpca The threshold used for the training of the FPCA. Only the maps for which the sum of the pixel is above this threshold are used for the training. If NULL, this threshold takes the value of threshold_classification.
#' @param list_search A list containing for each hyperparameters to be tested a list of the tested values.
#' @param seed An optional random seed
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param prototypes A set of l prototypes defining the Voronoï cells
#' @param distance_func  A function computing a distance between two elements in the output spaces.
#' @param return_pred A boolean indicating whether the predicted outputs should be returned or not
#' @param outputs_pred A list of the predicted outputs already obtained with the same parameters. Default is NULL.
#' @param only_positive A boolean indicating whether the predicted outputs should only contained positive values or not. Default is FALSE.
#' @param ncoeff The number of coefficients used for PCA
#' @param npc The number of principal components
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
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoï cell.
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.
#'
#' @return A list containing several outputs :
#' - list_search the list containing for each hyperparameters to be tested a list of the tested values.
#' - probas_pred_df a dataframe indicating for each combination of hyperparameters values the obtained predicted membership probabilities
#' - relative_error_df a dataframe indicating for each combination of hyperparameters values the relative error when predicting the membership probabilities
#' - outputs_pred an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.#' @export
#' @export
#' @import waveslim
#' @import foreach
#' @import GpOutput2D
#' @import rlibkriging
#' @import abind
#' @importFrom randomForest randomForest
#' @examples
#'  set.seed(5)
#'  func2D <- function(X){
#'  Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#'  n<-nrow(X)
#'  Y <- lapply(1:n, function(i){(X[i,2] > 0)*X[i,2]*X[i,1]*exp(-((0.8*Zgrid$z1+0.2*
#'  Zgrid$z2-10*X[i,1])**2)/(60*X[i,1]**2))*(
#'  Zgrid$z1-Zgrid$z2)*cos(X[i,1]*4)^2*sin(X[i,2]*4)^2})
#'  Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' library(randtoolbox)
#' design = as.data.frame(sobol(300,2))*2-1
#' outputs = func2D(design)
#' design_train = design[1:250,]
#' design_test = design[251:300,]
#' outputs_train = outputs[,,1:250]
#' outputs_test = outputs[,,251:300]
#' df_search = expand.grid(seq(0.1,1,0.3), c(1,5,9,13,17))
#' list_search = list("nodesize" = as.list(df_search[,2]), "classwt" = lapply(1:nrow(df_search),
#' function(i){c(df_search[i,1], 1-df_search[i,1])}))
#' density_ratio = rep(1,50)
#' prototypes = lapply(c(2,3,51,7), function(i){outputs[,,i]})
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}

#' list_rf_prob_train_test = rf_probas_training_test(design_train =
#' design_train,
#' design_test = design_test, outputs_train = outputs_train,
#' outputs_test = outputs_test, threshold_classification = 2,
#'  threshold_fpca = 0, list_search = list_search,
#' density_ratio = density_ratio, prototypes = prototypes,
#' distance_func= distance_func,
#' ncoeff = 400, npc = 6)

rf_probas_training_test = function(design_train, design_test, outputs_train, outputs_test,threshold_classification, threshold_fpca = NULL, list_search, density_ratio, prototypes, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},return_pred = FALSE, outputs_pred = NULL, only_positive = FALSE, ncoeff,npc, kernel="matern5_2", wf = "d4", boundary = "periodic",J=1,
                                   regmodel = "constant", normalize = FALSE, optim = "BFGS", objective = "LL", parameters = NULL,noise=FALSE, nugget = FALSE,seed = NULL,bias = NULL,...){
  if(is.null(seed)==FALSE){set.seed(seed)}
  if(is.null(threshold_fpca)){threshold_fpca = threshold_classification}
  bool_outputs_pred = outputs_pred
  probas_true = get_probas(density_ratio = density_ratio, data = outputs_test, prototypes = prototypes, distance_func = distance_func, cells = 1:length(prototypes), bias = bias)
  probas_pred_df = data.frame()
  relative_error_df = data.frame()

  sum_depth = Vectorize(function(it){sum(asub(x = outputs_train, idx = it, dims = length(dim(outputs_train)), drop = "selected"))})(1:dim(outputs_train)[length(dim(outputs_train))])
  y = as.factor(sum_depth > threshold_classification)
  indexes_train_fpca = which(sum_depth > threshold_fpca)
  if(is.null(bool_outputs_pred)){
    outputs_pred = list()
    fp = Fpca2d.Wavelets(asub(x = outputs_train, dims = length(dim(outputs_train)), idx = indexes_train_fpca,drop = FALSE), wf = wf, boundary = boundary, J = J, ncoeff = ncoeff, rank = npc) #We apply FPCA on the maps with water in the training group
    model = km_Fpca2d_lib(X = design_train[indexes_train_fpca,], response = fp,  kernel=kernel,
                      regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters,noise=noise,nugget=nugget)
  }
  for(i in 1:length(list_search[[1]])){
    if(is.null(bool_outputs_pred)){

      list_cv = list()
      for(v in 1:length(list_search)){
        list_cv[[v]] = list_search[[names(list_search)[v]]][[i]]
      }
      names(list_cv) = names(list_search)
      list_cv = c(list_cv, list("x" = design_train, "y" = y, "xtest" = design_test,...))
      rf = do.call(randomForest, list_cv)
      rf_pred = as.numeric(rf$test$predicted) - 1
      outputs_pred[[i]] = array(0,dim = dim(outputs_test))
      if(sum(rf_pred == 1)>0){
        pred =  matrix(sapply(1:npc, function(g){predict(object = model[[g]], x = design_test[rf_pred == 1,])$mean}), ncol = npc)
        outputs_pred_draft = inverse_Fpca2d(pred,fp)
        dimnames(outputs_pred[[i]]) = lapply(dim(outputs_pred[[i]]), function(i){1:i})
        dimnames(outputs_pred_draft) = c(lapply(dim(outputs_pred_draft)[-length(dim(outputs_pred_draft))], function(i){1:i}), list(which(rf_pred == 1)))
        afill(outputs_pred[[i]]) = outputs_pred_draft
      }
      if(only_positive){outputs_pred[[i]] = (outputs_pred[[i]] > 0)*outputs_pred[[i]]}
    }
    probas_pred_cv = get_probas(density_ratio = density_ratio, data = outputs_pred[[i]], prototypes = prototypes, distance_func = distance_func, cells = 1:length(prototypes), bias = bias)
    probas_pred_df = rbind(probas_pred_df,probas_pred_cv)
    relative_error = abs(probas_pred_cv - probas_true)/probas_true
    relative_error[probas_pred_cv == 0 & probas_true == 0] = 0
    relative_error[probas_pred_cv != 0 & probas_true == 0] = Inf
    relative_error_df = rbind(relative_error_df, relative_error)
  }
  colnames(relative_error_df) = c("ncoeff", "npc", 1:(ncol(relative_error_df)-2))
  colnames(probas_pred_df) = c("ncoeff", "npc", 1:(ncol(probas_pred_df)-2))
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(list_search = list_search,probas_pred = probas_pred_df, error = relative_error_df, outputs_pred = outputs_pred))
}
