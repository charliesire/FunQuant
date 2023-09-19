#' @title Check that all the elements in a list of prototypes are distinct
#'
#' @param prototypes A list of prototypes
#'
#' @return A boolean indicating whether the elements in prototypes are distinct
#' @export
#'
#' @examples
#' distinct_prototypes(list(1,2,34,1))
distinct_prototypes = function(prototypes){
  if(length(prototypes) == 1){return(TRUE)}
  else{
    for(i in 1:(length(prototypes)-1)){
      dist_prototypes = distance_to_prototypes(prototypes[[i]], lapply((i+1):length(prototypes), function(j){prototypes[[j]]}))$dist
      if(dist_prototypes == 0){return(FALSE)}
    }
    return(TRUE)
  }
}
