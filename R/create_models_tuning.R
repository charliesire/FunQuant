#' @title Computation of GP models in the PCA space for different ncoeff values
#'
#' @param outputs The training output samples on which the metamodel will be trained
#' @param ncoeff_vec a vector providing the different values of ncoeff that will be tested. ncoeff fixes the number of coefficients used for PCA.
#' @param npc Maximal number of principal components to be tested
#' @param design a data frame representing the design of experiments.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param wf name of the wavelet filter to use in the decomposition
#' @param boundary a character string which specifies how boundaries are treated. Only "periodic" is currently implemented .
#' @param J depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2). Default is 1.
#' @param kernel Character defining the covariance model: "exp", "gauss", "matern3_2", "matern5_2".
#' @param regmodel Universal Kriging linear trend: "constant", "linear", "interactive".
#' @param normalize Logical. If TRUE both the input matrix X and the response y in normalized to take values in the interval [0, 1].
#' @param optim Character giving the Optimization method used to fit hyper-parameters. Possible values are: "BFGS", "Newton" and "none", the later simply keeping the values given in parameters. The method "BFGS" uses the gradient of the objective. The method "Newton" uses both the gradient and the Hessian of the objective.
#' @param objective  Character giving the objective function to optimize. Possible values are: "LL" for the Log-Likelihood, "LOO" for the Leave-One-Out sum of squares and "LMP" for the Log-Marginal Posterior.
#' @param parameters Initial values for the hyper-parameters. When provided this must be named list with elements "sigma2" and "theta" containing the initial value(s) for the variance and for the range parameters. If theta is a matrix with more than one row, each row is used as a starting point for optimization.
#' @param noise Boolean specifying whether to execute NoiseKriging or not
#' @param nugget Boolean specifying whether to execute NuggetKriging or not
#' @return A list of lists of Kriging objects. The length of the output list is the length of ncoeff_vec. For each ncoeff of ncoeff_vec, a list of npc objects of class km is computed.
#' @import waveslim
#' @import foreach
#' @import rlibkriging
#' @export
#'
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*
#' Zgrid$z1+0.2*Zgrid$z2-10*X[i,])**2)/(60*X[i,]**2))*
#' (Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 8))
#' outputs = func2D(design)
#' ncoeff_vec = c(50,100,200,400)
#' npc = 4
#' models = create_models_tuning(outputs = outputs, ncoeff_vec = ncoeff_vec,
#' npc = 4, design = design)
create_models_tuning = function(outputs, ncoeff_vec, npc, design,kernel="matern5_2",
                                    regmodel = "constant", normalize = FALSE,wf = 'd4', boundary = "periodic", J = 1,
                                    optim = "BFGS", objective = "LL",
                                    parameters = NULL,noise=FALSE, nugget = FALSE){
  model_tuning = list()
  for(i in 1:length(ncoeff_vec)){
    ncoeff = ncoeff_vec[i]
    fp = Fpca2d.Wavelets(outputs, wf = wf, boundary = boundary, J = J, ncoeff = ncoeff, rank = npc)
    model_tuning[[i]] = km_Fpca2d_lib(X=design, response =fp, kernel=kernel,
                                      regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters,noise=noise,nugget = nugget)
  }
  return(model_tuning)
}


