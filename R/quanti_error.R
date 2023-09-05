#' @title Compututation of the empirical quantization error
#'

#' @param data The data that needs to be quantized. An array of any dimension is expected, the kth element must be selected with data[,..,k]. Useful only if inputs == NULL.
#' @param prototypes A set of prototypes. Useful only if cell_numbers == NULL.
#' @param density_ratio density_ratio indicates the weight fX/g of each data element.
#' @param distance_func A function computing a distance between two data elements. Useful only if cell_numbers == NULL.
#' @param batch_size If the computation is by batch, the number of elements by batch. Optional.
#' @param inputs If data is not provided, inputs is a dataframe of the inputs at which the outputs will be computed.
#' @param outputs_function Function to compute the outputs from the inputs.
#' @param method_IS The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell.
#' @param sampling_cells If method == "percell" and data is provided, sampling_cells is a vector indicating for each cell, the index element of data associated to this cell.
#' @return An estimation of the quantization error
#' @export
#' @import abind
#' @examples
#' prototypes = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3)))
#' data = array(runif(9*20)*20, dim = c(3,3,20))
#' density_ratio = rep(1,20)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' quanti_error(data = data, prototypes = prototypes, density_ratio = density_ratio,
#' distance_func = distance_func)
quanti_error = function(data = NULL, prototypes, density_ratio, batch_size = NULL, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}, inputs = NULL, outputs_function = NULL, method_IS = "unique", sampling_cells = 1:length(data)){
  drop = "selected"
  if(method_IS == "unique"){
    if(sum(dim(prototypes[[1]]) == 1) == length(dim(prototypes[[1]]))){drop = FALSE}

    if(is.null(batch_size)){
      if(is.null(data)){
        data = outputs_function(inputs)
      }
      distances = Vectorize(function(it){distance_to_prototypes(x = asub(x = data, dims = length(dim(data)), idx = it,drop = drop), prototypes = prototypes, distance_func = distance_func)$dist})(1:dim(data)[length(dim(data))])
      res = sqrt(mean(distances^2*density_ratio))
    }
    else{
      if(is.null(data)){nb_batch = nrow(inputs)%/%batch_size}
      else{nb_batch = length(data)}
      distances = c()
      for(batch in 1:nb_batch){
          if(is.null(data)){data_batch = outputs_function(inputs[((batch-1)*batch_size+1):(batch_size*batch),])}
          else{data_batch = data[[batch]]}
          distances = c(distances,Vectorize(function(it){distance_to_prototypes(x = asub(x = data_batch, dims = length(dim(data_batch)), idx = it,drop = drop), prototypes = prototypes, distance_func = distance_func)$dist})(1:dim(data_batch)[length(dim(data_batch))]))
      }
      res = sqrt(mean(distances^2*unlist(density_ratio)))
    }
  }
  else if(method_IS == "percell"){
    res = 0
      for(idx_sampling in unique(sampling_cells)){
        dists_cells =  t(Vectorize(function(it){distance_to_prototypes(x = asub(x = data[[idx_sampling]], dims = length(dim(data[[idx_sampling]])), idx = it,drop = drop), prototypes = prototypes, distance_func = distance_func)})(1:dim(data[[idx_sampling]])[length(dim(data[[idx_sampling]]))]))
        for(cell in which(sampling_cells == idx_sampling)){
          res = res + sum(as.numeric(dists_cells[dists_cells[,1] == cell,2])^2*density_ratio[[idx_sampling]][dists_cells[,1] == cell])/length(density_ratio[[idx_sampling]])
        }
      }
    res = sqrt(res)
  }
  return(res)
}
