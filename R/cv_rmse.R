rmse_loo = function(outputs, model_tuning = NULL, ncoeff_vec,npc_vec, return_pred = FALSE, formula = ~1,design, covtype="matern5_2",boundary = "periodic",J=1,
                      coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                      nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                      parinit = NULL, multistart=1,
                      kernel=NULL,control = NULL,type = "UK",...){
  dimnames(outputs) = NULL
  if(is.null(model_tuning)){model_tuning = create_models_tuning(outputs = outputs, ncoeff_vec = ncoeff_vec, npc = max(npc_vec), formula = formula,design = design, covtype=covtype,
                                                                coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                                                                nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                                                                parinit = parinit, multistart=multistart,
                                                                kernel=kernel,control = control,...)}
  grid_cv = expand.grid(ncoeff_vec, npc_vec)
  maps_rmse = list()
  outputs_loo_list = list()
  for(i in 1:nrow(grid_cv)){
    print(i)
    ncoeff = grid_cv[i,1]
    npc = grid_cv[i,2]
    indice_coeff = which(ncoeff_vec == ncoeff)
    fp = Fpca2d.Wavelets(outputs, wf = "d4", boundary = boundary, J = J, ncoeff = ncoeff, rank = npc)
    model = lapply(1:npc, function(k){model_tuning[[indice_coeff]][[k]]})
    pred = sapply(1:npc, function(i){leaveOneOut.km(model[[i]], type=type)$mean})
    outputs_loo = inverse_Fpca2d(pred,fp)
    if(return_pred){outputs_loo_list[[i]] = outputs_loo}
    err = (outputs_loo - outputs)^2
    maps_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
  }
  if(return_pred == FALSE){outputs_loo = NULL}
  return(list(grid_cv = grid_cv, maps_rmse = maps_rmse, outputs_pred = outputs_loo_list))
}

rmse_training_test = function(outputs_train,outputs_test, model_tuning = NULL, ncoeff_vec,npc_vec, return_pred = FALSE,formula = ~1,design_train, design_test, covtype="matern5_2",boundary = "periodic",J=1,
                                coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                                nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                                parinit = NULL, multistart=1,
                                kernel=NULL,control = NULL,type = "UK", ...){
  dimnames(outputs) = NULL
  if(is.null(model_tuning)){model_tuning = create_models_tuning(outputs = outputs_train, ncoeff_vec = ncoeff_vec, npc = max(npc_vec), formula = formula,design = design_train, covtype=covtype,
                                                                coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                                                                nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                                                                parinit = parinit, multistart=multistart,
                                                                kernel=kernel,control = control,...)}
  grid_cv = expand.grid(ncoeff_vec, npc_vec)
  relative_error_df = data.frame()
  outputs_pred_list = list()
  maps_rmse = list()
  for(i in 1:nrow(grid_cv)){
    ncoeff = grid_cv[i,1]
    npc = grid_cv[i,2]
    indice_coeff = which(ncoeff_vec == ncoeff)
    fp = Fpca2d.Wavelets(outputs_train, wf = "d4", boundary = boundary, J = J, ncoeff = ncoeff, rank = npc)
    model = lapply(1:npc, function(k){model_tuning[[indice_coeff]][[k]]})
    pred =  sapply(1:npc, function(k){predict(object = model[[k]], newdata = design_test, type = type, compute = FALSE)$mean})
    outputs_pred = inverse_Fpca2d(pred,fp)
    if(return_pred){outputs_pred_list[[i]] = outputs_pred}
    err = (outputs_pred - outputs)^2
    maps_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
  }
  return(list(grid_cv = grid_cv, maps_rmse = maps_rmse, outputs_pred = outputs_pred_list))
}

probas_k_fold = function(outputs, nb_folds, density_ratio, gamma, distance_func, model_tuning = NULL, ncoeff_vec,npc_vec, return_pred = FALSE,formula = ~1,design, covtype="matern5_2",boundary = "periodic",J=1,
                         coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                         nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                         parinit = NULL, multistart=1,
                         kernel=NULL,control = NULL,type = "UK",seed = NULL,bias = bias,...){
  set.seed(seed)
  length_dim = length(dim(outputs))
  dimnames(outputs) = NULL
  grid_cv = expand.grid(ncoeff_vec, npc_vec)
  folds = kfold(length(density_ratio), nb_folds)
  maps_rmse = list()
  outputs_pred_list = list()
  outputs_pred = list()
  relative_error_df = data.frame()
  for(k in 1:nb_folds){
    outputs_pred_list[[k]] = rmse_training_test(outputs_train = Subset(outputs,along = length(dim(outputs)), indices = which(folds != k)),outputs_test = Subset(outputs,along = length(dim(outputs)), indices = which(folds == k)), ncoeff_vec = ncoeff_vec,npc_vec = npc_vec, return_pred = TRUE,formula = formula,design_train = design[folds !=k,], design_test = design[folds == k,], covtype=covtype,boundary = boundary,J=J,
                                                  coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                                                  nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                                                  parinit = parinit, multistart=multistart,
                                                  kernel=kernel,control = control,type = type,bias = bias,...)$outputs_pred
  }
  for(i in 1:nrow(grid_cv)){
    outputs_pred[[i]] = array(NA, dim = dim(outputs))
    for(k in 1:nb_folds){
      outputs_pred[[i]][,,folds == k] = outputs_pred_list[[k]][[i]]
    }
    err = (outputs_pred[[i]] - outputs)^2
    maps_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
  }
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(grid_cv = grid_cv, maps_rmse = maps_rmse, outputs_pred = outputs_pred_list))
}


