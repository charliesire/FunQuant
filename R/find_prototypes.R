
#' @title The algorithm to find prototypes with multistart. Providing prototypes and their probability masses.

#' @param starting_proto Optional. If multistart = 1, starting_proto is a list of initial prototypes. Else, starting_proto is a list of list of prototypes, which length is equal to multistart.
#'
#' @param nb_cells Required only if starting_proto is NULL. Indicates the number of prototypes of the quantization.
#' @param data The data that needs to be quantized. An array of any dimension is expected, the kth element must be selected with data[,..,k]. If method = "percell", a list of data samples must be provided.
#' @param multistart Number of starts of the algorithm
#' @param method_IS The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell. Default is "unique".
#' @param sampling_cells If method == "percell" and data is provided, sampling_cells is a vector indicating for each cell, the index element of data associated to this cell.
#' @param density_ratio A vector indicating the weight fX/g of each data element. Default is a vector of 1. If method = "percell", a list of density_ratio must be provided, of length equal to the number of Voronoï cells.
#' @param budget The maximum number of iterations of the algorithm. Default is 10^3.
#' @param distance_func A function computing a distance between two data elements.
#' @param print_progress A boolean indicating whether to print the progress through the start numbers. Default is FALSE.
#' @param threshold A real positive number. When the distance between the new centroids and the previous ones is lower than this value, then we stop the algorithm.
#' @param trace A boolean. If TRUE, tracing information on the progress of the algorithm is produced. Default is FALSE.
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoi cell. Default is 0 for all Voronoi cells.
#' @param index_sampling_error Required only if method_IS = "percell". Indicates which of the data samples must be used for the computation of the quantization error.
#' @param all_starts A boolean indicating whether the function should return the optimal prototypes obtained for each start.
#' @param seed An optional random seed.
#' @param batch A boolean indicating whether the computations must be performed by batch or not. If TRUE, data, cell_numbers and density_ratio must be lists. Default is False.
#' @param density_function Density function of the inputs. Required only if method = "percell" and if the density biased function is updated at each iteration.
#' @param density_biased_function Density functions of the biased importance sampling. It must be a list of length equal to the number of Voronoi cells. Each element is a function of x and the inputs of the associated Voronoi cell.
#' @param inputs_ref Dataframe of inputs to delimit the Voronoi cells.
#' @param data_ref Outputs associated to inputs ref.
#' @param inputs_function List of functions to generate the inputs in each Voronoi cell with density equal to the elements of density biased function.
#' @param outputs_function Function to compute the outputs from the inputs.
#'
#' @return A list containing :
#' - prototypes : the list of optimal prototypes
#' - probas : a vector indicating the probability mass of the prototypes
#' - cell_numbers : a vector indicating the cell number associated to each data element
#' - iterations : an integer indicating the number of iterations performed
#' - record : a list containing all the centroids computed through the iterations of the best start. Provided only if trace = TRUE.
#' - all_errors : a vector indicating the quantization error of each start
#' - all_starts : a list indicating all the best prototypes obtained for each start. Provided only if all_start = TRUE.
#'
#' @export
#' @import abind
#' @examples
#' set.seed(20)
#' data = array(runif(9*20)*15, dim = c(3,3,20))
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' find_prototypes(nb_cells = 3, data = data,
#' multistart = 2, distance_func = distance_func)

find_prototypes = function(starting_proto = NULL, nb_cells = NULL, data = NULL, multistart = 1, method_IS = "unique", sampling_cells = NULL, density_ratio = rep(1, dim(data)[length(dim(data))]), budget = 10^3, threshold = 0, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},print_progress = FALSE, trace = FALSE, all_starts = FALSE,bias = NULL, index_sampling_error = 1, seed = NULL, batch = FALSE, density_function = NULL, density_biased_function = NULL, inputs_ref = NULL, data_ref = NULL, inputs_function = NULL, outputs_function = NULL){
  data_bool = is.null(data)
  if(is.null(dim(data)) & !is.null(data)){data = t(as.matrix(data))}
  if(!is.null(seed)){set.seed(seed)}
  best_error = +Inf
  best_list = list()
  all_errors = c()
  all_proto = list()
  for(it in 1:multistart){
    if(print_progress){print(paste("start number", it))}
    if(is.null(starting_proto)){
      if(is.null(nb_cells)){stop("nb_cells must me provided if starting_proto is not")}
      are_distinct = FALSE
      while(are_distinct == FALSE){
        prototypes_it = lapply(sample(x = 1:dim(data)[length(dim(data))], size = nb_cells), function(i){asub(data,dims = length(dim(data)), idx = i)})
        are_distinct = distinct_prototypes(prototypes_it)
      }
    }
    else{
      if(multistart == 1){prototypes_it = starting_proto}
      else{prototypes_it = starting_proto[[it]]}
    }
    prototypes_it = oned_to_matrix(prototypes_it)
    prototypes_it = sort_prototypes(prototypes_it)
    if(is.null(bias)){bias = rep(0, length(prototypes_it))}
    record = NULL
    if(trace){record = list(prototypes_it)}
    for (i in 1:budget){
      if(print_progress){print(paste("Iteration number", i))}
      if(method_IS == "unique"){
        if(batch){cell_numbers = lapply(data,function(x){get_cell_numbers(x,prototypes_it,distance_func = distance_func)})}
        else{cell_numbers = get_cell_numbers(data,prototypes_it,distance_func = distance_func)}
        if(length(table(unlist(cell_numbers)))<length(prototypes_it)){stop("One of the initial Voronoi cells is empty")}
      }
      else if(method_IS == "percell"){
        if(data_bool){
          numbers_ref = get_cell_numbers(data_ref,prototypes_it,distance_func = distance_func)
          cells_ref = lapply(1:length(prototypes_it), function(h){as.data.frame(inputs_ref[numbers_ref == h,])})
          inputs_list = lapply(1:length(prototypes_it), function(h){inputs_function[[h]](cells_ref[[h]])})
          data = lapply(1:length(prototypes_it), function(h){outputs_function(inputs_list[[h]])})
          sampling_cells = 1:length(prototypes_it)
          density_ratio = lapply(1:length(prototypes_it), function(h){compute_density_ratio(density_function, function(x){density_biased_function[[h]](x, cells_ref[[h]])}, inputs_list[[h]])})
          }
        cell_numbers = lapply(sampling_cells, function(h){get_cell_numbers(data[[h]],prototypes_it,distance_func = distance_func)})
        if(0 %in% sapply(1:length(cell_numbers), function(y){sum(cell_numbers[[y]] == y)}) ){stop("One of the initial Voronoi cells is empty")}
      }
      centro_probas = compute_centroids_and_proba(data = data, cell_numbers = cell_numbers,method_IS = method_IS, sampling_cells = sampling_cells, density_ratio = density_ratio, bias = bias, batch = batch)
      prototypes_new = oned_to_matrix(centro_probas$centroids)
      probas = centro_probas$probas
      prototypes_evol = sum(Vectorize(function(m){distance_func(prototypes_it[[m]], prototypes_new[[m]])})(1:length(prototypes_it)))
      prototypes_it = prototypes_new
      if(trace){record[[i]] = prototypes_it}
      if (prototypes_evol <= threshold){break}
    }
    if(method_IS == "percell"){
      batch_size = NULL
      error = quanti_error(data = data[[index_sampling_error]], prototypes = prototypes_it, density_ratio = density_ratio[[index_sampling_error]], batch_size = batch_size, distance_func = distance_func)
    }
    else if(method_IS == "unique"){
      if(batch){batch_size = dim(data[[1]])[length(dim(data[[1]]))]}
      else{batch_size = NULL}
      error = quanti_error(data = data, prototypes = prototypes_it, density_ratio = density_ratio, batch_size = batch_size, distance_func = distance_func)

    }
    all_errors = c(all_errors, error)
    if(all_starts){all_proto[[it]] = prototypes_it}
    if(error < best_error){
      best_error = error
      best_list = list(prototypes = prototypes_it, probas = probas, cell_numbers = cell_numbers, iterations = i, record = record)
    }
  }
  return(c(best_list, list(all_errors = all_errors, all_starts = all_proto)))
}

