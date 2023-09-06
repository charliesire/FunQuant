
#' @title Get the FPCA basis
#'
#' @param fpca an object of class \code{Fpca2d}.
#'
#' @return The FPCA basis maps
#'
#' @export
#'

find_basis = function(fpca){
  rot <- attr(fpca,"pca")$rotation #rotation matrix of the PCA : it is the eigen vectors Omega (dimension:  Ktilde x npc)

  idx_pca = attr(fpca, "idx_pca") # getting alpha coefficients that were not null and then used for the pca

  #res_coeff is the matrix Omega but with 0 on the 64^2 - Ktilde coefficients that were left out of the pca
  d = dim(mm$fp$EigFct)

  coeff <- attr(fpca,"coeff")
  mu <-as.vector(apply(coeff, c(1,2), mean))

  res_coeff <- sapply(1:ncol(rot), function(i){
    ci <- rep(0,prod(d[1:2]))
    ci[idx_pca] <- coeff_pca[,i]
    return(ci+mu)
  })# end res_coeff

  res_coeff <- sapply(1:ncol(rot), function(i){
    ci <- rep(0,prod(d))
    ci[idx_pca] <- rot[,i]
    return(ci)
  })# end res

  res_coeff <- array(res_coeff,dim=c(d[1:2], ncol(rot)))

  attr(res_coeff,"type")<-attr(fpca,"type")
  attr(res_coeff,"J")<-attr(fpca,"J")  # depth of wavelet decomposition
  attr(res_coeff,"wf")<-attr(fpca,"wf") #wavelet type
  attr(res_coeff,"boundary")<-attr(fpca,"boundary")
  attr(res_coeff,"dim")<-c(d,n)

  class(res_coeff)<-"Wavelet2D"
  bases_nb = Inverse2D(res_coeff)
}
