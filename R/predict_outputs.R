#' @title Predict outputs at new inputs, based on a training database.
#'
#' @param design_train a data frame representing the design of experiments of the training part.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param design_test a data frame representing the design of experiments on which we want to predict outputs.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param outputs_train The training output samples on which the metamodel will be trained
#' @param only_positive A boolean indicating whether the predicted outputs should only contained positive values or not. Default is FALSE.
#' @param seed An optional random seed
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
#' @param classification A boolean indicating whether a classification step should be integrated in the metamodel or not. Default is FALSE.
#' @param control_classification The list of hyperparameters of the classification function. Required only if classfification is TRUE.
#' @param threshold The threshold that creates the two classes of maps for the classification. Required only if classification is TRUE.
#' @param ... other parameters of \code{\link{km}} function from \code{DiceKriging}.

#'
#' @return An array containing the predicted outputs.
#' @export
#'
#' @examples
predict_outputs = function(design_train, design_test, outputs_train, only_positive = FALSE, seed = NULL, ncoeff,npc, formula = ~1, covtype="matern5_2",boundary = "periodic",J=1,
                           coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                           nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                           parinit = NULL, multistart=1,
                           kernel=NULL,control = NULL,type = "UK", classification = FALSE, control_classification = NULL,threshold = NULL,...){
  if(classification){
    sum_depth = Vectorize(function(it){sum(Subset(x = outputs_train, along = length(dim(outputs)), indices = it,drop = "selected"))})(1:dim(outputs)[length(dim(outputs))])
    y = as.factor(sum_depth > threshold)
    indexes_train_fpca = which(sum_depth > 0)
    control_classification = c(control_classification, list("x" = design_train, "y" = outputs_train, "xtest" = design_test))
    rf = do.call(randomForest, list_cv)
    outputs_train_fpca = Subset(x = outputs_train, along = length(dim(outputs)), indices = indexes_train_fpca,drop = FALSE)
    design_train_fpca = design_train[indexes_train_fpca,]
    rf_pred = as.numeric(rf$test$predicted) - 1
    design_test_fpca = design_test[rf_pred == 1,]
  }
  else{
    outputs_train_fpca = outputs_train
    design_train_fpca = design_train
    design_test_fpca = design_test
  }
  fp = Fpca2d.Wavelets(outputs_train_fpca, wf = "d4", boundary = boundary, J = J, ncoeff = ncoeff, rank = npc) #We apply FPCA on the maps with water in the training group
  model = km_Fpca2d(formula = formula, design = design_train_fpca, response = fp,  covtype=covtype,
                    coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                    nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                    parinit = parinit, multistart=multistart,
                    kernel=kernel,control = control)
  pred =  sapply(1:npc, function(g){predict(object = model, newdata = design_test_fpca, type = type, compute = FALSE)$mean})
  outputs_pred = inverse_Fpca2d(pred,fp)
  if(only_positive){outputs_pred = (outputs_pred > 0)*outputs_pred}
  return(outputs_pred)
}
