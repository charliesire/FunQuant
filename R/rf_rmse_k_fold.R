#' @title Computation of the Root Mean Square Error at each pixel of the outputs by kfold Cross Validation with a classification step, for different values of hyperparameters
#'
#' @param design A dataframe of inputs
#' @param outputs The output samples on which the metamodel will be trained and tested by kfold cross validation
#' @param threshold The threshold that creates the two classes of maps for the classification
#' @param list_search A list containing for each hyperparameters to be tested a list of the tested values.
#' @param nb_folds Number of folds
#' @param seed An optional random seed
#' @param return_pred A boolean indicating whether the predicted outputs should be returned or not
#' @param only_positive A boolean indicating whether the predicted outputs should only contained positive values or not. Default is FALSE.
#' @param ncoeff The number of coefficients used for PCA
#' @param npc The number of principal components
#' @param formula  an object of class "formula"
#' (or a list of "formula" which the length is equal to the number of modeled principal components)
#' specifying the linear trend of the kriging model (see \code{\link{lm}}) on each principal component.
#'  This formula should concern only the input variables (\code{design}), and not the output (\code{response}).
#'  The default is ~1, which defines a constant trend on each principal component.#'
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
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.

#'
#' @return A list containing several outputs :
#' - list_search the list containing for each hyperparameters to be tested a list of the tested values.
#' - outputs_rmse is a list of objects that have the same dimension as an output, obtained for each combination of hyperparameters values. Each element (called pixel here) of the objects is the RMSE computed between the predicted values of the pixel and the true value of the pixel.
#' - outputs_pred is an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.
#'#' @export
#'
#' @examples
rf_rmse_k_fold = function(design, outputs, threshold, list_search, nb_folds, return_pred = FALSE, seed = NULL, ncoeff,npc, formula = ~1, covtype="matern5_2",boundary = "periodic",J=1,
                          coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                          nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                          parinit = NULL, multistart=1,
                          kernel=NULL,control = NULL,type = "UK",...){
  if(is.null(seed)==FALSE){set.seed(seed)}
  folds = kfold(length(density_ratio), nb_folds)
  probas_true = get_probas(density_ratio = density_ratio, outputs = outputs, gamma = gamma, distance_func = distance_func, cells = 1:length(gamma), bias = bias)
  probas_pred_df = data.frame()
  relative_error_df = data.frame()
  folds = kfold(length(y), nb_folds)
  sum_depth = Vectorize(function(it){sum(Subset(x = outputs_train, along = length(dim(outputs)), indices = it,drop = "selected"))})(1:dim(outputs)[length(dim(outputs))])
  y = as.factor(sum_depth > threshold)
  list_indexes = lapply(1:length(list_search[[1]]), function(x){data.frame()})
  model = list()
  fp = list()
  for(k in 1:nb_folds){
    indexes_train = which(folds !=k)
    indexes_test = which(folds == k)
    indexes_train_fpca = which(sum_depth[indexes_train] > 0)
    fp[[k]] = Fpca2d.Wavelets(Subset(x = outputs, along = length(dim(outputs)), indices = indexes_train[indexes_train_fpca],drop = FALSE), wf = "d4", boundary = boundary, J = J, ncoeff = ncoeff, rank = npc) #We apply FPCA on the maps with water in the training group
    model[[k]] = km_Fpca2d(formula = formula, design = design[indexes_train[indexes_train_fpca],], response = fp[[k]],  covtype=covtype,
                           coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                           nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                           parinit = parinit, multistart=multistart,
                           kernel=kernel,control = control)
  }
  outputs_pred = list()
  for(i in 1:length(list_search[[1]])){
    list_cv = list()
    for(v in 1:length(list_search)){
      list_cv[[v]] = list_search[[names(list_search)[v]]][[i]]
    }
    names(list_cv) = names(list_search)
    outputs_pred_draft = list()
    for(k in 1:nb_folds){
      indexes_train = which(folds !=k)
      indexes_test = which(folds == k)
      indexes_train_fpca = which(sum_depth[indexes_train] > 0)
      list_cv_fold = c(list_cv, list("x" = design[folds!=k, ], "y" = y[folds!=k], "xtest" = design[folds==k,],...))
      rf = do.call(randomForest, list_cv_fold)
      rf_pred = as.numeric(rf$test$predicted) - 1
      list_indexes[[i]] = rbind(list_indexes[[i]], cbind(indexes_test[rf_pred == 1],k))
      pred =  sapply(1:npc, function(g){predict(object = model[[k]][[g]], newdata = design[indexes_test[rf_pred == 1],], type = type, compute = FALSE)$mean})

      outputs_pred_draft[[k]] = inverse_Fpca2d(pred,fp[[k]])
    }

    outputs_pred[[i]] = array(NA, c(dim(outputs)[-(length(dim(outputs)))],0))
    for(j in 1:dim(outputs)[(length(dim(outputs)))]){
      if(j %in% list_indexes[[i]][,1]){

        related_fold = list_indexes[[i]][list_indexes[[i]][,1] == j,2]
        outputs_pred[[i]] = abind(outputs_pred[[i]], Subset(x = outputs_pred_draft[[related_fold]], along = length(dim(outputs)), indices = which(list_indexes[[i]][list_indexes[[i]][,2] == related_fold,1] == j), drop = "selected"), along = length(dim(outputs)))
      }
      else{
        outputs_pred[[i]] = abind(outputs_pred[[i]], array(0, dim = dim(outputs)[-length(dim(outputs))]), along = length(dim(outputs)))
      }
    }
    dimnames(outputs_pred[[i]]) = NULL
    if(only_positive){outputs_pred[[i]] = (outputs_pred[[i]] > 0)*outputs_pred[[i]]}
    err = (outputs_pred[[i]] - outputs)^2
    outputs_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
    if(return_pred == FALSE){outputs_pred = list()}
  }
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(list_search = list_search, outputs_rmse = outputs_rmse, outputs_pred = outputs_pred))
}

