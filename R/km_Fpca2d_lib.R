###########################################################################################

#' @title Gaussian Process Model on principal components of \code{Fpca2d},
#' by using \code{rlibkriging} package.
#'
#' @description  the function \code{Kriging} of the rlibkriging package is use to
#' fit kriging models on each principal component.
#'
#' @author SIRE Charlie,Tran Vi-vi Elodie PERRIN
#'
#' @param X a data frame representing the design of experiments.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param response an object of class \code{Fpca2d} which contains eigen
#' decomposition of the model/function ouput.
#' @param kernel Character defining the covariance model: "exp", "gauss", "matern3_2", "matern5_2".
#' @param regmodel Universal Kriging linear trend: "constant", "linear", "interactive".
#' @param normalize Logical. If TRUE both the input matrix X and the response y in normalized to take values in the interval [0, 1].
#' @param optim Character giving the Optimization method used to fit hyper-parameters. Possible values are: "BFGS", "Newton" and "none", the later simply keeping the values given in parameters. The method "BFGS" uses the gradient of the objective. The method "Newton" uses both the gradient and the Hessian of the objective.
#' @param objective  Character giving the objective function to optimize. Possible values are: "LL" for the Log-Likelihood, "LOO" for the Leave-One-Out sum of squares and "LMP" for the Log-Marginal Posterior.
#' @param parameters Initial values for the hyper-parameters. When provided this must be named list with elements "sigma2" and "theta" containing the initial value(s) for the variance and for the range parameters. If theta is a matrix with more than one row, each row is used as a starting point for optimization.
#' @param noise Boolean specifying whether to execute NoiseKriging or not
#' @param nugget Boolean specifying whether to execute NuggetKriging or not
#' @import foreach
#' @importFrom graphics image
#' @importFrom rlibkriging Kriging NuggetKriging NoiseKriging
#'
#' @return a list of Kriging objects for each modeled principal component.
#'
#' @examples
#'
#'  set.seed(5)
#'  func2D <- function(X){
#'  Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#'  n<-nrow(X)
#'  Y <- lapply(1:n, function(i){(X[i,2] > 0)*X[i,2]*X[i,1]*
#'  exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-10*X[i,1])**2)/(60*X[i,1]**
#'  2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,1]*4)^2*sin(X[i,2]*4)^2})
#'  Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' library(randtoolbox)
#' X = as.data.frame(sobol(300,2))*2-1
#' Y = func2D(X)
#'  fp = Fpca2d.Wavelets(Y, wf = "d4", boundary = "periodic", J = 1, ncoeff = 200, rank = 3)
#' mm <- km_Fpca2d_lib(X=X,response=fp)
#' @export
km_Fpca2d_lib <- function(X,response,kernel="matern5_2",
                          regmodel = "constant", normalize = FALSE, optim = "BFGS", objective = "LL", parameters = NULL,noise=FALSE, nugget = FALSE){
  X=as.matrix(X)
  if(!is(response,"Fpca2d")){
    stop("response must be an object of class 'Fpca2d'.")
  }# end if

  # FPCA outputs & number of principal components
  y <- response$x; nPC <- ncol(y)


  # if only one principal component
  if((nPC==1)|(is.null(nPC))){
    if(noise){
      m = NoiseKriging(X=X, y=y, kernel=kernel,
                  regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters)
    }
    else if(nugget){
      m<-NuggetKriging(X=X, y=y, kernel=kernel,
                 regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters)
    }
    else{
    m<-Kriging(X=X, y=y, kernel=kernel,
               regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters)
    }
    m<-list(m)
  }else{

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #
  #     Repeat function input for each principal component
  #
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   #_________________
   # to vector
   #_________________

  test <- kernel
    if((!is.null(test)) & (is.vector(test))){
      kernel<-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
                   test)
    }
  test <- regmodel
  if((!is.null(test)) & (is.vector(test))){
    regmodel<-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
      test)
  }# end if

  test <- normalize
  if((!is.null(test)) & (is.vector(test))){
    normalize<-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
      test)
  }# end if

  test <- optim
  if((!is.null(test)) & (is.vector(test))){
    optim<-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
      test)
  }# end if

  test <- objective
  if((!is.null(test)) & (is.vector(test))){
    objective<-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
      test)
  }# end if

  test <- parameters
  if((!is.null(test)) & (is.vector(test))){
    parameters<-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
      test)
  }# end if
  #__________
  # to list
  #__________

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  #%%%%%%%%%%%%%
  ### Models ###
  #%%%%%%%%%%%%%
  if(noise){
    m = foreach(i = 1:nPC) %do% {NoiseKriging(X=X, y=as.numeric(y[,i]), kernel=kernel[[i]],
                      regmodel = regmodel[[i]], normalize = normalize[[i]], optim = optim[[i]], objective = objective[[i]], parameters = parameters[[i]])
    }
  }
  else if(nugget){
    m<-foreach(i = 1:nPC) %do% {NuggetKriging( X=X, y=as.numeric(y[,i]),kernel=kernel[[i]],
                regmodel = regmodel[[i]], normalize = normalize[[i]], optim = optim[[i]], objective = objective[[i]], parameters = parameters[[i]])
    }
  }
  else{
    m<-foreach(i = 1:nPC) %do% {Kriging(X=X, y=as.numeric(y[,i]), kernel=kernel[[i]],
               regmodel = regmodel[[i]], normalize = normalize[[i]], optim = optim[[i]], objective = objective[[i]], parameters = parameters[[i]])
    }
  }

  } # end ifelse
  class(m)<-"km_Fpca2d_lib"
  attr(m,"fpca")<-response
  return(m)
} # end km.Fpca2d
