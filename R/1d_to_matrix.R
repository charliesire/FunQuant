#' @title Transform every 1D prototype into a column matrix
#'
#' @param prototypes a list of prototypes
#'
#' @return A list of matrices prototypes
#' @export
#'
#' @examples
#' oned_to_matrix(list(1:5, runif(5), rep(0,5)))
oned_to_matrix = function(prototypes){
  prototypes = lapply(1:length(prototypes), function(j){
    if(is.null(dim(prototypes[[j]])) | length(dim(prototypes[[j]])) == 1){t(as.matrix(prototypes[[j]]))}
    else{prototypes[[j]]}})
  return(prototypes)
}
