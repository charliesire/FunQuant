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
rf_proba_k_fold = function(x, outputs, threshold, list_search, nb_folds, density_ratio, gamma, distance_func,return_pred = FALSE, seed = NULL, ncoeff,npc, formula = ~1, covtype="matern5_2",boundary = "periodic",J=1,
                          coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                          nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                          parinit = NULL, multistart=1,
                          kernel=NULL,control = NULL,type = "UK",bias = NULL,...){
  if(is.null(seed)==FALSE){set.seed(seed)}
  folds = kfold(length(density_ratio), nb_folds)
  probas_true = get_probas(density_ratio = density_ratio, outputs = outputs, gamma = gamma, distance_func = distance_func, cells = 1:length(gamma), bias = bias)
  probas_pred_df = data.frame()
  relative_error_df = data.frame()
  #folds = kfold(length(y), nb_folds)
  folds = df_folds_no_breach_rf[, "folds"]
  drop = "selected"
  if(sum(dim(gamma[[1]]) == 1) == length(dim(gamma[[1]]))){drop = F}
  sum_depth = Vectorize(function(it){sum(Subset(x = outputs, along = length(dim(outputs)), indices = it,drop = drop))})(1:dim(outputs)[length(dim(outputs))])
  y = as.factor(sum_depth > threshold)
  list_indexes = lapply(1:length(list_search[[1]]), function(x){data.frame()})
  model = list()
  fp = list()
  # for(k in 1:nb_folds){
  #   print(k)
  #   indexes_train = which(folds !=k)
  #   indexes_test = which(folds == k)
  #   indexes_train_fpca = which(sum_depth[indexes_train] > 0)
  #   fp[[k]] = Fpca2d.Wavelets(Subset(x = outputs, along = length(dim(outputs)), indices = indexes_train[indexes_train_fpca],drop = FALSE), wf = "d4", boundary = boundary, J = J, ncoeff = ncoeff, rank = npc) #We apply FPCA on the maps with water in the training group
  #   model[[k]] = km_Fpca2d(formula = formula, design = x[indexes_train[indexes_train_fpca],], response = fp[[k]],  covtype=covtype,
  #                          coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
  #                          nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
  #                          parinit = parinit, multistart=multistart,
  #                          kernel=kernel,control = control)
  # }
  #mm <<- model
  model = mm
  fp = ff
  #ff <<- fp
  outputs_pred = list()
  maps_rmse = list()
  for(i in 1:length(list_search[[1]])){

    print(i)
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
      list_cv_fold = c(list_cv, list("x" = x[folds!=k, ], "y" = y[folds!=k], "xtest" = x[folds==k,], "ytest" = y[folds == k],...))
      lll <<- list_cv_fold
      rf = do.call(randomForest, list_cv_fold)
      rf_pred = as.numeric(rf$test$predicted) - 1
      list_indexes[[i]] = rbind(list_indexes[[i]], cbind(indexes_test[rf_pred == 1],k))
      pred =  sapply(1:npc, function(g){predict(object = model[[k]][[g]], newdata = x[indexes_test[rf_pred == 1],], type = type, compute = FALSE)$mean})
      pp <<- pred
      fff <<- fp[[k]]

      xx <<- x[indexes_test[rf_pred == 1],]
      outputs_pred_draft[[k]] = inverse_Fpca2d(pred,fp[[k]])
    }
    abcde <<- outputs_pred_draft
    outputs_pred[[i]] = array(NA, c(dim(outputs)[-(length(dim(outputs)))],0))
    for(j in 1:dim(outputs)[(length(dim(outputs)))]){
      if(j %in% list_indexes[[i]][,1]){
        ll <<- list_indexes[[i]]
        jj <<- j
        related_fold = list_indexes[[i]][list_indexes[[i]][,1] == j,2]
        outputs_pred[[i]] = abind(outputs_pred[[i]], Subset(x = outputs_pred_draft[[related_fold]], along = length(dim(outputs)), indices = which(list_indexes[[i]][list_indexes[[i]][,2] == related_fold,1] == j), drop = "selected"), along = length(dim(outputs)))
      }
      else{
        outputs_pred[[i]] = abind(outputs_pred[[i]], array(0, dim = dim(outputs)[-length(dim(outputs))]), along = length(dim(outputs)))
      }
    }
    dimnames(outputs_pred[[i]]) = NULL
    outputs_pred[[i]] = (outputs_pred[[i]] > 0)*outputs_pred[[i]]
    probas_pred_cv = get_probas(density_ratio = density_ratio, outputs = outputs_pred[[i]], gamma = gamma, distance_func = distance_func, cells = 1:length(gamma), bias = bias)
    probas_pred_df = rbind(probas_pred_df,probas_pred_cv)
    relative_error_df = rbind(relative_error_df, abs(probas_pred_cv - probas_true)/probas_true)
    if(return_pred == FALSE){outputs_pred = list()}
  }
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(list_search = list_search,probas_pred = probas_pred_df, error = relative_error_df, outputs_pred = outputs_pred))
}

