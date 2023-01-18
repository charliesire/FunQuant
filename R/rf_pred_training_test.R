#' @title Predict the output class of a validation sample with random forest classifier for different hyperparameters values
#'
#' @param xtrain A dataframe of training inputs
#' @param xtest A dataframe of test inputs
#' @param ytrain A vector of training outputs for the random forests
#' @param list_search A list containing for each hyperparameters to be tested a list of the tested values.
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.
#'
#' @return A list containing a vector of predicted classes for each combination of tested hyperparameters.
#' @export
#'
#' @examples
rf_pred_training_test = function(xtrain, xtest, ytrain, list_search,...){
  pred = list()
  for(i in 1:length(list_search[[1]])){
    liste_cv = list()
    for(var in names(list_search)){
      liste_cv = c(liste_cv, var = list_search[[var]][[i]])
    }
    liste_cv = c(liste_cv,list("x" = xtrain, "y" = ytrain, "xtest" = xtest,...))
    rf = do.call(randomForest, liste_cv)
    pred[[i]] = rf$test$predicted
  }
  return(pred)
}


