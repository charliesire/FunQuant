library(ClimProjDiags)
library(ggplot2)
library(miceadds)
library(abind)
library(POT)
library(randomForest)
library(dismo)


path = getwd()
source(file = "C:/Users/sire-cha/Downloads/test_package/coastal_utils.R")
source.all(path ="C:/Users/sire-cha/Documents/FunQuant/R")
#
# #load("R/data_example/NewFitting_Charlie_v090821.RData")
# pretty_map_bis = function(gamma_k, Hmax = NULL, Hmin = 0){
#   x = NULL
#   y = NULL
#   he = NULL
#   c = 0
#   for (i in 1:length(lon)){
#     for (j in 1:length(lat)){
#       c<-c+1
#       x[c]<-lon[i]
#       y[c]<-lat[j]
#       he[c]<-gamma_k[i,j]
#     }
#   }
#
#   dat_lonlat<-data.frame(lon=x,lat=y,he=he)
#   dat_lonlat[which(he == 0), "he"] = NA
#   p =  ggplot()  +
#     geom_raster(data = dat_lonlat , aes(x = lon, y = lat, fill = he, alpha = he)) +
#     geom_contour(data = dat_lonlat, aes(x = lon, y = lat,z = he), breaks = c(0.15,10)) +
#     scale_alpha(range = c(0.6,1), guide = "none") + labs(fill = "Water depth") +
#     scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
#   if(is.null(Hmax)){
#     p = p + scale_fill_gradient(low = "white", high = "blue", na.value = "transparent")
#   }
#   else{
#     p = p + scale_fill_gradient(low = "white", high = "blue", na.value = "transparent", limits = c(Hmin,Hmax))
#
#   }
#   return(p)
# }
#
#
# load("~/math4flood/Data_Cas_Maritime/Boucholeurs_NOFailure_Charlie.RData")
# doe_no_breach = doe ## inputs without breach
# HE_no_breach = HE ## maps without breach
# area_no_breach = Area ## water areas without breach
#
# doe_no_breach_bis = data.frame(doe_no_breach)
# doe_no_breach_bis$loc = 0
# doe_no_breach_bis$erosion = 0
#
# load("~/math4flood/Data_Cas_Maritime/Boucholeurs_Failure_Charlie.RData")
# doe_breach = data.frame(doe) ## inputs with breach
# HE_breach = HE ## maps with breach
# area_breach = Area ## water areas with breach
#
# colnames(doe_breach) = colnames(doe_no_breach_bis)
# doe_full = rbind(doe_no_breach_bis, doe_breach) ## concatenate with and without breach
#
#
# HE_full = abind(HE_no_breach, HE_breach, along = 3)
#
# HE_full = HE_full[41:(41+63), 70:(70+63),] ## Select only a zone 64x64 of the maps
#
# sum_depth = Vectorize(function(i){sum(HE_full[,,i])})(1:1300) ## Vector indicating for each map the water volume
#
# doe_full_phy = data.frame(doe_full) ## doe_full_phy are the inputs not normalized
#
# doe_full_phy[1] = 0.95 + (3.7-0.95)*doe_full[1]
# doe_full_phy[2] = 0.65 + (2.5-0.65)*doe_full[2]
# doe_full_phy[3] = -6.1 + (6.1+6.1) *doe_full[3]
# doe_full_phy[4] = -12.4 + (-0.5+12.4)*doe_full[4]
# doe_full_phy[5] = 0.2 + (12.2 - 0.2)*doe_full[5]
#
# outputs = HE_full
# gamma = lapply(1:5, function(i){outputs[,,i]})
#
# distL2 = function(A1,A2){
#   return(sqrt(sum((A1-A2)^2)))
# }
#
#
# # list_loi_exp is a list containing the empirical observation of the variables related to the offshore conditions
# list_loi_exp = list()
# list_loi_exp[["T"]] = T #tide
# list_loi_exp[["S"]] = S #surge
# list_loi_exp[["phi"]] = phi #phi
# list_loi_exp[["t1"]] = t1 #t1
# list_loi_exp[["t2"]] = t2 #t2
# integrate = stats::integrate
#
#
# offset = 0.65 # the offset represent an increase of the Surge
#
#
# proba_for_gpd = mean(list_loi_exp[["S"]] <= 0.55)
#
#
# ## The first computation of the density function is directly related to the empirical density. For the surge, we add an information : the surge above 0.55 + offset follows a gpd density
#
# densite_i = function(x, i){
#   res = 1
#   if (i %in% c(1,3,4,5)){
#     res = res * density(list_loi_exp[[i]], from = x, to = x, n = 1)$y[1]}
#   else if (i == 2){res = res *  (density(list_loi_exp[[i]],  from = x-offset, to = x-offset, n = 1)$y[1]*(x-offset  <= 0.55) + dgpd(x-offset,loc = 0.55, scale = fit$fitted.values[1], shape = fit$fitted.values[2])*(x-offset > 0.55)*(1 - proba_for_gpd)) }
#   return(res)
# }
# ## We then perform a truncation of the density function, so that when all historical observations are contained in the interval use to build the training database,
# ## we consider that the support of the density function is contained in this interval
#
# ## df_bornes contains the intervals for each input used to build the 1300 training maps
# df_bornes = data.frame(lower = c(0.95,0.65,-6.1, -12.4, 0.2), upper = c(3.7,2.5,6.1,-0.5,12.2))
#
# ##intervals_truncation indicates the bounds of the truncation. For the tide, no truncation is performed for the lower bound, as the associated interval for the training database is [0.95, 3.7], although empirical observations are recorded below 0.95
# intervals_truncation = list(c(-20,df_bornes[1,2]), c(offset, df_bornes[2,2]), c(df_bornes[3,1],df_bornes[3,2]),c(df_bornes[4,1],df_bornes[4,2]), c(df_bornes[5,1],df_bornes[5,2]))
#
# coef_truncation = c()
# for(i in 1:5){
#   lower = intervals_truncation[[i]][1]
#   upper = intervals_truncation[[i]][2]
#   coef_truncation = c(coef_truncation,min(integrate(f = Vectorize(function(x){densite_i(x, i)}), lower = lower, upper = upper)$value,1))
# }
#
#
#
#
# #densite_i_trunc = function(x,i){
# #if (i %in% c(1,3,4,5)){res = density(list_loi_exp[[i]], from = x, to = x, n = 1)$y[1]*(x>=intervals_truncation[[i]][1])*(x< intervals_truncation[[i]][2])/coef_truncation[i]
# #}
# #else if (i == 2){
# #res = density(list_loi_exp[[2]],  from = x-offset, to = x-offset, n = 1)$y[1]*(x - offset <= 0.55) + dgpd(x,loc = 0.55, scale = 0.1112987, shape = 0.2670383)*(x - offset > 0.55)*(1 - proba_for_gpd)*(x>=intervals_truncation[[i]][1])*(x<=intervals_truncation[[i]][2])/coef_truncation[i]
# #}
# #return(res)
# #}
# densite_i_trunc = function(x,i){
#   densite_i(x,i)*(x>=intervals_truncation[[i]][1])*(x<= intervals_truncation[[i]][2])/coef_truncation[i]
#
# }
#
#
# ## densite_xu computes the density of the 5 variables related to the offshore conditions, considered independantly
# densite_xu = function(x){
#   res = 1
#   for (i in c(1:5)){
#     res = res * densite_i_trunc(x[i],i)
#   }
#   return(as.numeric(res))
# }
#
# ## Regarding the variables related to the breach, we consider that the probability of a breach depends on the signals of the tide and the surge.
# ## If the maximum of the sum of these signals is above 70% of the average size of the embankment, then the probability of breach is 1/2. Else it is 10^-4
#
# #sum_water_signal sums the signals of tide and surge and takes the maximum
# sum_water_signal = function(x){
#   signal_T = x[1]*cos(2*pi/12.5*seq(-12,12,l=1001))
#   function_S = function(t){
#     t_pic = x[3]
#     if(t > t_pic){return(max(0,x[2]*(1-(t-t_pic)/x[5])))}
#     else{return(return(max(0,x[2]*(1+(t_pic-t)/x[4]))))}
#   }
#   signal_S = Vectorize(function_S)(seq(-12,12,l=1001))
#   return(max(signal_S+signal_T))
# }
#
# ## Densite_ratio_full computes the probability f_{X} / nu
#
# densite_ratio_full = function(x){
#   q_x = 1/((3.7-0.95)*(2.5-0.65)*(6.1+6.1)*(-0.5+12.4)*(12.2-0.2))
#   qte_eau = sum_water_signal(x)
#   if (qte_eau <= 0.7*4.76){
#     p_rup = 10^-4
#     if(x[7]==0){
#       return(densite_xu(x[1:5])/q_x*(1-p_rup)/(5/13))
#     }
#     else{
#       return(densite_xu(x[1:5])/q_x*p_rup/(8/13))
#     }
#   }
#   if (qte_eau > 0.7*4.76){
#
#     p_rup = 1/2
#
#     if(x[7]==0){
#       return(densite_xu(x[1:5])/q_x*(1-p_rup)/(5/13))
#     }
#     else{
#       return(densite_xu(x[1:5])/q_x*p_rup/(8/13))
#     }
#   }
# }
#
#
# densite_ratio_full_vec = function(X){
#   return(Vectorize(function(i){densite_ratio_full(as.numeric(X[i,]))})(1:nrow(X)))
# }
#
# density_ratio = densite_ratio_full_vec(doe_full_phy)
#
# gamma = lapply(c(501,502,503,504,506), function(i){outputs[,,i]})
#
#
# quanti_unif = proto_map_algo(gamma = gamma, outputs = outputs, density_ratio = rep(1,1300),print_progress = T,method_IS = "unique", distance_func = distL2)
#
# cste_bias = -integrate(Vectorize(function(x){density(list_loi_exp[[1]], from = x, to = x, n = 1)$y[1]}), lower = -1, upper = df_bornes[1,1])$value
#
# quanti_proba = proto_map_algo(gamma = gamma, outputs = outputs, density_ratio = density_ratio,print_progress = TRUE,method_IS = "unique", distance_func = distL2, bias = c(cste_bias, rep(0,4)))
#
#
# lat = 1:64
# lon = 1:64
#
# for(i in 1:5){print(pretty_map_bis(quanti_unif$gamma[[i]]))}
#
# set.seed(1)
# aa = sample(cut(1:10, breaks = 3, labels = FALSE))
#
# set.seed(1)
# bb = kfold(1:10,3)

classwt_df = data.frame(weight_false = c(1,0.9,0.1), weight_true = c(1,0.1,0.9))
for(i in 2:8){
  classwt_df = rbind(classwt_df, c(1/(i+1), i/(i+1)))
  classwt_df = rbind(classwt_df, c(i/(i+1), 1/(i+1)))
}

df_params = cbind(classwt_df, nodesize = rep(1, nrow(classwt_df)))
for(i in 2:8){df_params = rbind(df_params, cbind(classwt_df, nodesize = rep(i, nrow(classwt_df))))}

densite_vec = Vectorize(function(it){
  return(densite_ratio_full(as.numeric(doe_full_phy[it,])))})(1:1300)

ll_search = list(classwt = lapply(1:nrow(df_params), function(i){c(df_params[i,1],df_params[i,2])}), nodesize = lapply(1:nrow(df_params), function(i){c(df_params[i,3])}))
ll_search = list(classwt = list(c(0.8,0.2)), nodesize = list(4))
dimnames(HE_full) = NULL
cste_bias = -integrate(Vectorize(function(x){density(list_loi_exp[[1]], from = x, to = x, n = 1)$y[1]}), lower = -1, upper = df_bornes[1,1])$value

distL2 = function(A1,A2){
  return(sqrt(sum((A1-A2)^2)))
}

load(file = "gamma_star_1300maps.RData")
nb_breaks = 10
df_folds_no_breach_rf = data.frame(indice = 1:500) #we create 10 folds
df_folds_no_breach_rf$folds = 0

set.seed(1)
df_folds_no_breach_rf[, "folds"] = sample(cut(1:500, breaks = nb_breaks, labels = FALSE))
pp =rf_proba_k_fold(x = doe_full[1:500,1:5],outputs = HE_full[,,1:500], threshold = 5,list_search = ll_search, nb_folds = 10, density_ratio = densite_vec[1:500], gamma = gamma_star, distance_func = distL2,ncoeff = NULL, npc = 2, nugget = 0.1,multistart=5, lower = rep(0.1,5), control = list(trace = 0), bias = c(cste_bias, rep(0,4)))

