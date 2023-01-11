#' @title Computation of the IS coefficients of variation of the membership probability for different set of prototypes.
#'
#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma_list A list of gamma on which we want to evaluate the IS coefficient of variation of the membership probability. Each gamma is a set of prototypes
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param distance_func  A function computng a distance between two elements in the output spaces.
#' @param cells The Voronoï cell numbers that we are investigating
#' @param nv The size of the sample for which we want to estimate the IS coefficient of variation of the membership probability
#' @param cell_numbers An optional list providing for each set of prototypes the voronoi cell number of every output.
#'
#' @return A list of IS coefficients of variation of the membership probability obtained for each set of prototypes.
#' @export
#'
#' @examples
std_proba = function(outputs, gamma_list, density_ratio, distance_func, cells, cell_numbers = NULL, nv){
  std_list = list()
  for(it in 1:length(gamma_list)){ #for all Gamma in (Gamma^r)

    if(is.null(cell_numbers)){cell_numbers_it = get_cell_numbers(outputs, gamma_list[[it]], distance_func)}
    else{cell_numbers_it = cell_numbers[[it]]}
    list_for_std = Vectorize(function(i){density_ratio*(cell_numbers_it == i)})(cells)
    std_list[[it]] = apply(list_for_std, 2, std)/apply(list_for_std[[it]],2,mean)/sqrt(nv) #for all voronoi cells of Gamma, we compute the relative standard error
  }
  return(std_list)
}

