#' @title Train the different metamodel steps
#'
#' @param design_train a data frame representing the design of experiments of the training part.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param outputs_train The training output samples on which the metamodel will be trained
#' @param seed An optional random seed
#' @param ncoeff The number of coefficients used for PCA
#' @param npc The number of principal components
#' @param kernel Character defining the covariance model: "exp", "gauss", "matern3_2", "matern5_2".
#' @param regmodel Universal Kriging linear trend: "constant", "linear", "interactive".
#' @param normalize Logical. If TRUE both the input matrix X and the response y in normalized to take values in the interval [0, 1].
#' @param optim Character giving the Optimization method used to fit hyper-parameters. Possible values are: "BFGS", "Newton" and "none", the later simply keeping the values given in parameters. The method "BFGS" uses the gradient of the objective. The method "Newton" uses both the gradient and the Hessian of the objective.
#' @param objective  Character giving the objective function to optimize. Possible values are: "LL" for the Log-Likelihood, "LOO" for the Leave-One-Out sum of squares and "LMP" for the Log-Marginal Posterior.
#' @param parameters Initial values for the hyper-parameters. When provided this must be named list with elements "sigma2" and "theta" containing the initial value(s) for the variance and for the range parameters. If theta is a matrix with more than one row, each row is used as a starting point for optimization.
#' @param classification A boolean indicating whether a classification step should be integrated in the metamodel or not. Default is FALSE.
#' @param control_classification The list of hyperparameters of the classification function. Required only if classfification is TRUE.
#' @param threshold_classification The threshold that creates the two classes of maps for the classification
#' @param threshold_fpca The threshold used for the training of the FPCA. Only the maps for which the sum of the pixel is above this threshold are used for the training. If NULL, this threshold takes the value of threshold_classification.
#' @param wf name of the wavelet filter to use in the decomposition
#' @param boundary a character string which specifies how boundaries are treated. Only "periodic" is currently implemented .
#' @param J depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2). Default is 1.
#' @param noise Boolean specifying whether to execute NoiseKriging or not
#' @param nugget Boolean specifying whether to execute NuggetKriging or not
#'
#' @return An list containing :
#' - a trained randomForest object
#' - a trained Fpca2d object
#' - a list of trained km object
#' @export
#' @import waveslim
#' @import foreach
#' @import rlibkriging
#' @examples
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
#' design = as.data.frame(sobol(250,2))*2-1
#' outputs = func2D(design)
#' fit_metamodel = fit_metamodel(design_train = design, outputs_train = outputs,
#' ncoeff = 400, npc = 6, classification = TRUE,
#' control_classification = list(nodesize = 4), threshold_classification = 2)

fit_metamodel = function(design_train, outputs_train, seed = NULL, ncoeff,npc,kernel="matern5_2", wf = "d4", boundary = "periodic",J=1,
                             regmodel = "constant", normalize = FALSE, optim = "BFGS", objective = "LL", parameters = NULL,noise=FALSE, nugget = FALSE, classification = FALSE, control_classification = NULL,threshold_classification = 0,threshold_fpca = NULL){
  if(is.null(threshold_fpca)){threshold_fpca = threshold_classification}
  pred_fpca = TRUE
  if(classification){
    sum_depth = Vectorize(function(it){sum(asub(x = outputs_train, idx = it, dims = length(dim(outputs_train)), drop = "selected"))})(1:dim(outputs_train)[length(dim(outputs_train))])
    y = as.factor(sum_depth > threshold_classification)
    indexes_train_fpca = which(sum_depth > threshold_fpca)
    control_classification = c(control_classification, list("x" = design_train, "y" = y))
    rf = do.call(randomForest, control_classification)
    outputs_train_fpca = asub(x = outputs_train, dims = length(dim(outputs_train)), idx = indexes_train_fpca,drop = FALSE)
    design_train_fpca = design_train[indexes_train_fpca,]
  }
  else{
    outputs_train_fpca = outputs_train
    design_train_fpca = design_train
    rf = NULL
  }
  fp = Fpca2d.Wavelets(outputs_train_fpca, wf = wf, boundary = boundary, J = J, ncoeff = ncoeff, rank = npc) #We apply FPCA on the maps with water in the training group
  model = km_Fpca2d_lib(X=design_train_fpca, response=fp, kernel=kernel,
                        regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters,noise=noise,nugget=nugget)

  return(list(classifier = rf, fp = fp, model = model))
}

