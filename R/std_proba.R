#' @title Computation of the IS coefficients of variation of the membership probability for different set of prototypes.
#'
#' @param data The data that needs to be quantized. An array of any dimension is expected, the kth element must be selected with data[,..,k].
#' Useful only if cell_numbers == NULL and outputs_function == NULL.
#' @param prototypes_list A list of set of prototypes on which we want to evaluate the IS centroid standard deviation. Each element is a list of prototypes.
#' @param density_ratio density_ratio indicates the weight fX/g of each data element.
#' @param distance_func  A function computing a distance between two data elements.
#' @param cells The Voronoï cell numbers that we are investigating
#' @param nv The size of the sample for which we want to estimate the IS coefficient of variation of the membership probability
#' @param cell_numbers An optional list providing for each set of prototypes the voronoi cell number of every data element.
#' @param outputs_function Function to compute the outputs from the inputs.
#' @param inputs If data and cell_numbers are not provided, inputs is a dataframe of the inputs at which the outputs will be computed.
#' @param batch_size If the computation is by batch, the number of elements by batch. Optional.
#' @param return_cell_numbers Boolean indicating whether the cell_numbers must be provided or not.
#' @param bootstrap Integer indicating the number of bootstramp samples to generate. If NULL, no boostrap is performed. Default is NULL.
#' @return A list of IS coefficients of variation of the membership probability obtained for each set of prototypes.
#' @export
#' @import abind
#' @importFrom  stats var
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2
#' -10*X[i,])**2)/(60*X[i,]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 50))
#' data = func2D(design)
#' prototypes_list = list(lapply(c(1,3,10,14,18), function(i){data[,,i]}))
#' density_ratio = rep(1, 50)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' list_std_proba = std_proba(data = data, prototypes_list = prototypes_list,
#'  density_ratio = density_ratio, distance_func = distance_func,
#'  cells = 1:length(prototypes_list[[1]]), nv = 50)
std_proba = function(data = NULL, prototypes_list, density_ratio = rep(1,dim(data)[length(data)]), distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},cells, cell_numbers = NULL, nv, outputs_function = NULL, inputs = NULL, batch_size = NULL, return_cell_numbers = FALSE, bootstrap = NULL){
  std_list = list()
  if(is.null(data) & is.null(cell_numbers)){
    nb_batch = nrow(inputs)%/%batch_size
    cell_numbers = vector(mode='list', length=length(prototypes_list))
    for(batch in 1:nb_batch){
      data = outputs_function(inputs[((batch-1)*batch_size+1):(batch_size*batch),])
      for(it in 1:length(prototypes_list)){
        cell_numbers[[it]] = c(cell_numbers[[it]],get_cell_numbers(data = data, prototypes = prototypes_list[[it]], distance_func = distance_func))
      }
    }
  }
  for(it in 1:length(prototypes_list)){
    if(is.null(cell_numbers)){cell_numbers_it = get_cell_numbers(data = data, prototypes = prototypes_list[[it]], distance_func = distance_func)}
    else{cell_numbers_it = cell_numbers[[it]]}
    if(!is.null(bootstrap)){
      df_probas_boot = data.frame()
      estim_probas = get_probas(density_ratio = density_ratio, method_IS = "unique", cell_numbers = cell_numbers_it, distance_func = distance_func, cells = cells)
      for(boot in 1:bootstrap){
        idxs = sample(1:length(cell_numbers_it), size = length(cell_numbers_it), replace = TRUE)
        cell_numbers_boot = cell_numbers_it[idxs]
        df_probas_boot = rbind(df_probas_boot,get_probas(density_ratio = density_ratio[idxs], method_IS = "unique",cell_numbers = cell_numbers_boot, distance_func = distance_func, cells = cells))
      }
      std_list[[it]] = as.numeric(apply(df_probas_boot, 2, function(x){sqrt(var(x))})*sqrt(length(cell_numbers_it))/sqrt(nv)/estim_probas)
    }
    else{
    df_for_std = Vectorize(function(i){density_ratio*(cell_numbers_it == i)})(cells)
    std_list[[it]] = as.numeric(apply(df_for_std, 2, function(x){sqrt(var(x))})/apply(df_for_std,2,mean)/sqrt(nv)) #for all voronoi cells, we compute the relative standard error
    }
  }
  if(!return_cell_numbers){return(std_list)}
  else{return(list(std_list = std_list, cell_numbers = cell_numbers))}

}
