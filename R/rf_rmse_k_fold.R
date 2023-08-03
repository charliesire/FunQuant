#' @title Computation of the Root Mean Square Error at each pixel of the outputs by kfold Cross Validation with a classification step, for different values of hyperparameters
#'
#' @param design A dataframe of inputs
#' @param outputs The output samples on which the metamodel will be trained and tested by kfold cross validation
#' @param threshold_classification The threshold that creates the two classes of maps for the classification
#' @param threshold_fpca The threshold used for the training of the FPCA. Only the maps for which the sum of the pixel is above this threshold are used for the training. If NULL, this threshold takes the value of threshold_classification.
#' @param list_search A list containing for each hyperparameters to be tested a list of the tested values.
#' @param nb_folds Number of folds
#' @param seed An optional random seed
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
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.

#'
#' @return A list containing several outputs :
#' - list_search the list containing for each hyperparameters to be tested a list of the tested values.
#' - outputs_rmse is a list of objects that have the same dimension as an output, obtained for each combination of hyperparameters values. Each element (called pixel here) of the objects is the RMSE computed between the predicted values of the pixel and the true value of the pixel.
#' - outputs_pred is an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.
#' @export
#' @import waveslim
#' @import foreach
#' @import GpOutput2D
#' @import rlibkriging
#' @import abind
#' @importFrom randomForest randomForest
#' @importFrom dismo kfold
#' @examples
#' set.seed(5)
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){(X[i,2] > 0)*X[i,2]*X[i,1]*
#' exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-10*X[i,1])**2)/(60*X[i,1]**2))*
#' (Zgrid$z1-Zgrid$z2)*cos(X[i,1]*4)^2*sin(X[i,2]*4)^2})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' library(randtoolbox)
#' design = as.data.frame(sobol(100,2))*2-1
#' outputs = func2D(design)
#' list_search = list("nodesize" = as.list(c(1,3,5,7,9,11)))
#' list_rf_rmse_k_fold = rf_rmse_k_fold(design = design,outputs = outputs,
#'  threshold_classification = 2, threshold_fpca = 0,
#'  list_search = list_search, nb_folds = 10, ncoeff = 400,
#'   npc = 6)

rf_rmse_k_fold = function(design, outputs, threshold_classification, threshold_fpca = NULL, list_search, nb_folds, return_pred = FALSE, outputs_pred = NULL, only_positive = FALSE, ncoeff,npc, kernel="matern5_2", wf = "d4", boundary = "periodic",J=1,
                          regmodel = "constant", normalize = FALSE, optim = "BFGS", objective = "LL", parameters = NULL,noise=FALSE, nugget = FALSE,seed = NULL,...){
  if(is.null(seed)==FALSE){set.seed(seed)}
  if(is.null(threshold_fpca)){threshold_fpca = threshold_classification}
  bool_outputs_pred = outputs_pred
  folds = kfold(dim(outputs)[length(dim(outputs))], nb_folds)
  sum_depth = Vectorize(function(it){sum(asub(x = outputs, dims = length(dim(outputs)), idx = it,drop = "selected"))})(1:dim(outputs)[length(dim(outputs))])
  y = as.factor(sum_depth > threshold_classification)
  model = list()
  fp = list()
  if(is.null(bool_outputs_pred)){
    for(k in 1:nb_folds){
      outputs_pred = list()
      indexes_train = which(folds !=k)
      indexes_test = which(folds == k)
      indexes_train_fpca = which(sum_depth[indexes_train] > 0)
      fp[[k]] = Fpca2d.Wavelets(asub(x = outputs, dims = length(dim(outputs)), idx = indexes_train[indexes_train_fpca],drop = FALSE), wf = wf, boundary = boundary, J = J, ncoeff = ncoeff, rank = npc) #We apply FPCA on the maps with water in the training group
      model[[k]] = km_Fpca2d_lib(X = design[indexes_train[indexes_train_fpca],], response = fp[[k]],   kernel=kernel,
                             regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters,noise=noise,nugget=nugget)
    }
  }
  outputs_rmse = list()
  for(i in 1:length(list_search[[1]])){
    if(is.null(bool_outputs_pred)){
      outputs_pred[[i]] = array(0,dim = dim(outputs))
      dimnames(outputs_pred[[i]]) = lapply(dim(outputs_pred[[i]]), function(i){1:i})
      list_cv = list()
      for(v in 1:length(list_search)){
        list_cv[[v]] = list_search[[names(list_search)[v]]][[i]]
      }
      names(list_cv) = names(list_search)
      outputs_pred_draft = list()
      for(k in 1:nb_folds){
        indexes_train = which(folds !=k)
        indexes_test = which(folds == k)
        indexes_train_fpca = which(sum_depth[indexes_train] > threshold_fpca)
        list_cv_fold = c(list_cv, list("x" = design[folds!=k, ], "y" = y[folds!=k], "xtest" = design[folds==k,],...))
        rf = do.call(randomForest, list_cv_fold)
        rf_pred = as.numeric(rf$test$predicted) - 1
        if(sum(rf_pred == 1)>0){
          pred =  matrix(sapply(1:npc, function(g){predict(object = model[[k]][[g]], x = design[indexes_test[rf_pred == 1],])$mean}), ncol = npc)
          outputs_pred_draft = inverse_Fpca2d(pred,fp[[k]])
          dimnames(outputs_pred_draft) = c(lapply(dim(outputs_pred_draft)[-length(dim(outputs_pred_draft))], function(i){1:i}), list(indexes_test[rf_pred == 1]))
          afill(outputs_pred[[i]]) = outputs_pred_draft
        }
      }

      if(only_positive){outputs_pred[[i]] = (outputs_pred[[i]] > 0)*outputs_pred[[i]]}
    }
    err = (outputs_pred[[i]] - outputs)^2
    outputs_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
  }
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(list_search = list_search, outputs_rmse = outputs_rmse, outputs_pred = outputs_pred))
}

