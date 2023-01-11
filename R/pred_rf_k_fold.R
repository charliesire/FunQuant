#' @title Predict the output class with random forest classifier by kfold cross validation for different hyperparameters values
#'
#' @param x A dataframe of inputs
#' @param y A vector of categorical output
#' @param list_search A list containing for each hyperparameters to be tested a list of the tested values.
#' @param nb_folds Number of folds
#' @param seed An optional random seed
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.
#'
#' @return A list containing a vector of predicted classes for each combination of tested hyperparameters.
#' @export
#'
#' @examples
pred_rf_k_fold = function(x, y, list_search, nb_folds, seed = NULL,...){
  if(is.null(seed)==FALSE){set.seed(seed)}
  folds = kfold(length(y), nb_folds)
  pred = list()
  for(i in 1:length(list_search[[1]])){
    liste_cv = list()
    pred[[i]] = rep(0, length(y))
    for(var in names(list_search)){
      liste_cv = c(liste_cv, var = list_search[[var]][[i]])
    }
      for(k in 1:nb_folds){
        liste_cv_fold = c(liste_cv, list("x" = x[folds!=k, ], "y" = y[folds!=k], "xtest" = x[folds==k,], "ytest" = y[folds == k],...))
        rf = do.call(randomForest, liste_cv_fold)
        pred[[i]][which(folds == k)] = rf$test$predicted
      }
  }
  return(pred)
}
