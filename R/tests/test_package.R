
library(miceadds)
library(waveslim)
library(abind)
library(DiceKriging)
library(foreach)
library(dismo)
library(randomForest)
library(ggplot2)
library(evd)
library(randtoolbox)

source.all("~/FunQuant/R")
source.all("~/FunQuant/R/GpOutput2D-main/GpOutput2D/R")

func2D <- function(X){
  Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
  n<-nrow(X)
  Y <- lapply(1:n, function(i){(X[i,2] > 0)*X[i,2]*X[i,1]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-10*X[i,1])**2)/(60*X[i,1]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,1]*4)^2*sin(X[i,2]*4)^2})
  Ymaps<- array(unlist(Y),dim=c(20,20,n))
  return(Ymaps)
}
f2 = function(x){
  res = 0
  ptrunc = pgev(-1, loc=-0.4,scale=0.1) + 1 - pgev(1, loc=-0.4,scale=0.1)
  if(x>=-1 & x < 1){res = dgev(x, loc=-0.4,scale=0.1)/(1-ptrunc)}
  return(res)
}
f1 = function(x){1/2}
fX = function(x){f1(x[1])*f2(x[2])}
g = function(x){1/4}
design = as.data.frame(sobol(400,2))*2-1
density_ratio = compute_density_ratio(f = fX, g = g, inputs = design)
distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}

outputs = func2D(design)

## Proto map algo sans m?tamod?le
sum_depth = Vectorize(function(i){sum(outputs[,,i])})(1:dim(outputs)[3])
gamma = lapply(1:5, function(i){outputs[,,which.min(abs(as.numeric(quantile(sum_depth,c(0,0.6,0.7,0.8,0.9))[i]) - sum_depth))[1]]})
res_proto = proto_map_algo(gamma = gamma, outputs, density_ratio = density_ratio, distance_func = distance_func, trace = TRUE)
gamma_star_apriori = res_proto$gamma

## K fold with rmse

list_rmse_k_fold = rmse_k_fold(outputs = outputs[,, sum_depth > 0], design = design[sum_depth > 0,], nb_folds = 10, npc_vec = npc_vec, ncoeff_vec = ncoeff_vec, seed = 10, control = list(trace = FALSE))
quantile_90 = sapply(list_rmse_k_fold$outputs_rmse, function(x){quantile(x, 0.9)})

best_params = list_rmse_k_fold$grid_cv[quantile_90 == min(quantile_90),]

## K fold with probas

npc_vec = 2:6
ncoeff_vec = c(50,150,250,350,400)
list_probas_k_fold = probas_k_fold(outputs = outputs[,, sum_depth > 0], design = design[sum_depth > 0,], nb_folds = 10, density_ratio = density_ratio[sum_depth >0], gamma = gamma_star_apriori, distance_func = distance_func, npc_vec = npc_vec, ncoeff_vec = ncoeff_vec, seed = 10, control = list(trace = FALSE))

## Random forest rmse

df_search = expand.grid(seq(0.1,0.9,0.2), c(1,3,5,7))
list_search = list("nodesize" = as.list(df_search[,2]), "classwt" = lapply(1:nrow(df_search), function(i){c(df_search[i,1], 1-df_search[i,1])}))
list_rf_rmse_k_fold = rf_rmse_k_fold(design = design,outputs = outputs, threshold = 0, list_search = list_search, nb_folds = 10, ncoeff = 250, npc = 4, control = list(trace = FALSE), seed = 10)
quantile_90 = sapply(list_rf_rmse_k_fold$outputs_rmse, function(x){quantile(x, 0.9)})

best_params_rf = df_search[quantile_90 == min(quantile_90),]

## Random forest probas

list_rf_probas_k_fold = rf_probas_k_fold(design = design,outputs = outputs, density_ratio = density_ratio, gamma = gamma_star_apriori, distance_func = distance_func,threshold = 0, list_search = list_search, nb_folds = 10, ncoeff = 250, npc = 4, control = list(trace = FALSE), seed = 10)

## Predict

set.seed(10)
design_test = as.data.frame(matrix(runif(2*10^4), ncol=2)*2-1)
outputs_test = func2D(design_test)
outputs_pred = predict_outputs(design_train = design, design_test = design_test, outputs_train = outputs, seed = 10, ncoeff = 250, npc = 4, control = list(trace = FALSE), classification = TRUE, threshold = 0)
density_ratio_pred = compute_density_ratio(fX,g,design_test)

res_proto_pred = proto_map_algo(gamma = gamma_star_apriori, outputs_pred, density_ratio = density_ratio_pred, distance_func = distance_func, trace = FALSE, print_progress = TRUE)

## std probas

res_std_probas = std_proba(outputs_pred, gamma_list = list(res_proto_pred$gamma), density_ratio = density_ratio_pred, distance_func = distance_func, cells = 1:5, nv = 10^4)

## std_centroid

res_std_centroid = std_centroid(outputs = outputs_pred, gamma_list = list(res_proto_pred$gamma), density_ratio = density_ratio_pred, distance_func = distance_func, cells = 1:5, nv = 10^4)

quantiles_90_std = sapply(res_std_centroid[[1]], function(x){quantile(x, 0.9)})

## Compare probas true and probas pred

probas_true = get_probas(density_ratio = density_ratio_pred, gamma = res_proto_pred$gamma, outputs = outputs_test, distance_func = distance_func)
probas_pred = res_proto_pred$probas

## Compare quantization errors

res_proto_test = proto_map_algo(gamma = gamma_star_apriori, outputs_test, density_ratio = density_ratio_pred, distance_func = distance_func, trace = FALSE, print_progress = TRUE)
quanti_error_test = quanti_error(outputs = outputs_test, gamma = res_proto_test$gamma, density_ratio = density_ratio_pred, distance_func = distance_func)
quanti_error_pred = quanti_error(outputs = outputs_test, gamma = res_proto_pred$gamma, density_ratio = density_ratio_pred, distance_func = distance_func)


