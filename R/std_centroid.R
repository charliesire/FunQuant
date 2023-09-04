#' @title Computation of the IS centroid standard deviation for different sets of prototypes.
#'
#' @param data The data samples that need to be quantized. An array of any dimension is expected, the kth element must be selected with data[,..,k]. Useful only if cell_numbers == NULL.
#' @param prototypes_list A list of set of prototypes on which we want to evaluate the IS centroid standard deviation. Each element is a list of prototypes.
#' @param density_ratio density_ratio indicates the weight fX/g of each data element.
#' @param distance_func  A function computing a distance between two data elements.
#' @param cells The Voronoï cell numbers that we are investigating.
#' @param nv The size of the sample for which we want to estimate the IS centroid standard deviation.
#' @param cell_numbers An optional list providing for each set of prototypes the voronoi cell number of every data element.
#' @param outputs_function Function to compute the outputs from the inputs.
#' @param inputs If data and cell_numbers are not provided, inputs is a dataframe of the inputs at which the outputs will be computed.
#' @param batch_size If the computation is by batch, the number of elements by batch. Optional.
#' @param return_cell_numbers Boolean indicating whether the cell_numbers must be provided or not.
#' @param bootstrap Integer indicating the number of bootstramp samples to generate. If NULL, no boostrap is performed. Default is NULL.

#' @return A list providing for each set of prototypes a list the IS centroid standard deviation for each voronoi cell
#' @export
#' @import abind
#' @importFrom stats var cov
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-
#' 10*X[i,])**2)/(60*X[i,]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 50))
#' data = func2D(design)
#' prototypes_list = list(lapply(c(1,3,10,14,18), function(i){data[,,i]}))
#' density_ratio = rep(1, 50)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' list_std_centroid = std_centroid(data = data, prototypes_list =
#' prototypes_list, density_ratio = density_ratio, distance_func = distance_func
#' , cells = 1:length(prototypes_list[[1]]), nv = 50)
std_centroid = function(data = NULL, prototypes_list, density_ratio = rep(1,dim(data)[length(data)]), distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},cells, cell_numbers = NULL, nv,outputs_function = NULL, inputs = NULL, batch_size = nrow(inputs), return_cell_numbers = FALSE, bootstrap = NULL){
  bool_cell_numbers = is.null(cell_numbers)
  std_ratio_list = list()
  if(is.null(data)){
    list_density_voronoi = vector(mode='list', length=length(prototypes_list))
    list_mean_maps = lapply(1:length(prototypes_list),function(it){as.list(rep(0,length(cells)))})
    list_covs = lapply(1:length(prototypes_list),function(it){as.list(rep(0,length(cells)))})
    list_vars= lapply(1:length(prototypes_list),function(it){as.list(rep(0,length(cells)))})
    nb_batch = nrow(inputs)%/%batch_size
    if(bool_cell_numbers){cell_numbers = vector(mode='list', length=length(prototypes_list))}
    for(batch in 1:nb_batch){
      data = outputs_function(inputs[((batch-1)*batch_size+1):(batch_size*batch),])
      weighted_map = t(matrix(data, nrow = prod(dim(data)[-length(dim(data))]),ncol = dim(data)[length(dim(data))]))*density_ratio[((batch-1)*batch_size+1):(batch_size*batch)]
      for(it in 1:length(prototypes_list)){
        if(bool_cell_numbers){cell_numbers[[it]] = c(cell_numbers[[it]],get_cell_numbers(data = data, prototypes = prototypes_list[[it]], distance_func = distance_func))}
        sum_map = list()
        sum_map_square = list()
        sum_map_covs = list()
        for(cell in 1:length(cells)){
          j=cells[cell]
          map_loop = weighted_map
          map_loop[cell_numbers[[it]][((batch-1)*batch_size+1):(batch_size*batch)] != j,] = 0
          sum_map[[cell]] = apply(map_loop, 2, sum)
          sum_map_square[[cell]] = apply(map_loop, 2, function(x){sum(x^2)})
          sum_map_covs[[cell]] = apply(map_loop,2,function(x){sum(x*density_ratio[((batch-1)*batch_size+1):(batch_size*batch)])})
        }
        list_covs[[it]] = mapply("+", list_covs[[it]], sum_map_covs, SIMPLIFY = FALSE)
        list_mean_maps[[it]] = mapply("+", list_mean_maps[[it]], sum_map, SIMPLIFY = FALSE)
        list_vars[[it]] = mapply("+", list_vars[[it]], sum_map_square, SIMPLIFY = FALSE)
      }
    }
    for(it in 1:length(prototypes_list)){
      prototypes = prototypes_list[[it]]
      list_density_voronoi[[it]] = Vectorize(function(i){density_ratio*(cell_numbers[[it]] == i)})(cells)
      std_ratio_list[[it]] = vector(mode='list', length=length(cells))
      for(cell in 1:length(cells)){
        j=cells[cell]
        covariance = 1/(nrow(inputs)-1)*(list_covs[[it]][[cell]]-mean(list_density_voronoi[[it]][,cell])*list_mean_maps[[it]][[cell]])/nv
        moy1 = 1/nrow(inputs)*list_mean_maps[[it]][[cell]]
        moy2 = mean(list_density_voronoi[[it]][,cell]) #This is the empirical expectation of B
        var1 = 1/(nrow(inputs)-1)*(list_vars[[it]][[cell]]-1/nrow(inputs)*list_mean_maps[[it]][[cell]]^2)/nv #This is the empirical variance of Ai for all i
        var2 = var(list_density_voronoi[[it]][,cell])/nv
        std_ratio_list[[it]][[j]] = sqrt(1/moy2^2*(var1+var2/moy2^2*moy1^2-2*covariance*moy1/moy2)) #This is the empirical standard error of the ratio
      }
    }
  }
  else{
    weighted_map = t(matrix(data, nrow = prod(dim(data)[-length(dim(data))]),ncol = dim(data)[length(dim(data))]))*density_ratio
    for(it in 1:length(prototypes_list)){
      prototypes = prototypes_list[[it]]
      std_ratio_list[[it]] = vector(mode='list', length=length(cells))
      if(is.null(cell_numbers)){cell_numbers_it = get_cell_numbers(data = data, prototypes = prototypes, distance_func = distance_func)}
      else{cell_numbers_it = cell_numbers[[it]]}
      if(!is.null(bootstrap)){
        list_centro_boot = list()
        for(boot in 1:bootstrap){
          idxs = sample(1:length(cell_numbers_it), size = length(cell_numbers_it), replace = TRUE)
          cell_numbers_boot = cell_numbers_it[idxs]
          data <<- asub(data,dims = length(dim(data)), idx = idxs)
          density_ratio <<- density_ratio[idxs]
          cell_numbers <<- cell_numbers_boot
          aa <<- compute_centroids_and_proba(data = asub(data,dims = length(dim(data)), idx = idxs), density_ratio = density_ratio[idxs], method_IS = "unique",cell_numbers = cell_numbers_boot, cells = cells)$centroids
          list_centro_boot[[boot]] = compute_centroids_and_proba(data = asub(data,dims = length(dim(data)), idx = idxs), density_ratio = density_ratio[idxs], method_IS = "unique",cell_numbers = cell_numbers_boot, cells = cells)$centroids
        }
        for(cell in 1:length(cells)){
          j=cells[cell]
          std_ratio_list[[it]][[cell]] = array(0, dim = dim(list_centro_boot[[1]][[1]]))
          for(pix in 1:length(std_ratio_list[[it]][[cell]])){
            std_ratio_list[[it]][[cell]][pix] = sqrt(var(sapply(lapply(list_centro_boot, function(x){x[[cell]]}), function(h){h[pix]})))*sqrt(dim(data)[length(dim(data))])/sqrt(nv)
          }
        }
      }
      else{
        for(cell in 1:length(cells)){#for all voronoi cells
          j = cells[cell]
          map_loop = weighted_map #weighted map is the set of maps multiplied by the weights f_{x}/mu
          density_voronoi = density_ratio*(cell_numbers_it == j) #density_voronoi is the vector of the weights multiplied by the characteristic function of the voronoi cell
          map_loop[cell_numbers_it != j,] = 0 #map_loop is now the set of maps multiplied by the weights f_{x}/mu multiplied by the characteristic function of the voronoi cell
          covariance = apply(map_loop,2,function(x){cov(x,density_voronoi)})/nv #this is the covariance between Ai and B for all i
          moy1 = apply(map_loop,2, mean) #This is the empirical expectation of Ai for all i
          moy2 = mean(density_voronoi) #This is the empirical expectation of B
          var1 = apply(map_loop,2,var)/nv #This is the empirical variance of Ai for all i
          var2 = var(density_voronoi)/nv #This is the empirical variance of B
          std_ratio_list[[it]][[cell]] = sqrt(1/moy2^2*(var1+var2/moy2^2*moy1^2-2*covariance*moy1/moy2)) #This is the empirical standard error of the ratio
        }
      }
    }
  }
  if(return_cell_numbers){return(list(std_ratio_list = std_ratio_list, cell_numbers = cell_numbers))}
  else{return(std_ratio_list)}
}



