#' @title Predict outputs at new inputs, based on a training database.
#'
#' @param metamodel_fitted Optional list containing the different metamodel steps fitted on the training data. It contains a classifier, a Fpca2d object, and a list of km objects. Can be obtained with the function fit_metamodel.
#' @param design_train a data frame representing the design of experiments of the training part.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param design_test a data frame representing the design of experiments on which we want to predict outputs.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param outputs_train The training output samples on which the metamodel will be trained
#' @param only_positive A boolean indicating whether the predicted outputs should only contained positive values or not. Default is FALSE.
#' @param seed An optional random seed
#' @param ncoeff The number of coefficients used for PCA
#' @param npc The number of principal components
#' If a vector, the length should be equal to the number of modeled principal components.
#' @param wf name of the wavelet filter to use in the decomposition
#' @param boundary a character string which specifies how boundaries are treated. Only "periodic" is currently implemented .
#' @param J depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2). Default is 1.
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
#' @param noise Boolean specifying whether to execute NoiseKriging or not
#' @param nugget Boolean specifying whether to execute NuggetKriging or not
#' @return An array containing the predicted outputs.
#' @export
#' @import waveslim
#' @import foreach
#' @rawNamespace import(GpOutput2D, except = Fpca2d.Wavelets)
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
#' design = as.data.frame(sobol(300,2))*2-1
#' outputs = func2D(design)
#' design_train = design[1:250,]
#' design_test = design[251:300,]
#' outputs_train = outputs[,,1:250]
#' outputs_test = outputs[,,251:300]
#' outputs_pred = predict_outputs(design_train = design_train,
#' design_test = design_test, outputs_train = outputs_train,
#' ncoeff = 400, npc = 6, classification = TRUE,
#' control_classification = list(nodesize = 4), threshold_classification = 2)

predict_outputs = function(metamodel_fitted = NULL, design_train = NULL, design_test, outputs_train = NULL, only_positive = FALSE, seed = NULL, ncoeff = NULL,npc = NULL,wf = "d4", boundary = "periodic",J=1,kernel="matern5_2",
                               regmodel = "constant", normalize = FALSE, optim = "BFGS", objective = "LL", parameters = NULL,noise=FALSE,
                               nugget = FALSE,classification = FALSE, control_classification = NULL,threshold_classification = NULL,
                               threshold_fpca = NULL){
  colnames(design_test) = colnames(design_train)
  if(is.null(metamodel_fitted)){
    metamodel_fitted = fit_metamodel(design_train = design_train, outputs_train = outputs_train, seed = seed, ncoeff = ncoeff, npc = npc,
                                         kernel=kernel, wf = wf, boundary = boundary,J=J,
                                         regmodel = regmodel, normalize = normalize, optim = optim, objective = objective, parameters = parameters,noise=noise, nugget = nugget, classification = classification, control_classification = control_classification,threshold_classification = threshold_classification,threshold_fpca = threshold_fpca)
  }
  rf = metamodel_fitted$classifier
  fp = metamodel_fitted$fp
  model = metamodel_fitted$model
  pred_fpca = TRUE
  if(classification){
    rf_pred = as.numeric(predict(rf, design_test)) - 1
    design_test_fpca = design_test[rf_pred == 1,]
    if(sum(rf_pred) == 0){pred_fpca = FALSE}
  }
  else{
    design_test_fpca = design_test
  }

  if(pred_fpca){
    pred =  matrix(sapply(1:length(model), function(g){predict(object = model[[g]], x = design_test_fpca)$mean}), ncol = length(model))
    outputs_pred_draft = inverse_Fpca2d(pred,fp)
  }
  outputs_pred = array(0,dim = c(dim(fp$EigFct)[-length(dim(fp$EigFct))], nrow(design_test)))
  if(classification == FALSE){outputs_pred = outputs_pred_draft}
  else if(pred_fpca){
    outputs_pred = array(0,dim = c(dim(fp$EigFct)[-length(dim(fp$EigFct))], nrow(design_test)))
    dimnames(outputs_pred) = lapply(dim(outputs_pred), function(i){1:i})
    dimnames(outputs_pred_draft) = c(lapply(dim(outputs_pred_draft)[-length(dim(outputs_pred_draft))], function(i){1:i}), list(which(rf_pred == 1)))
    afill(outputs_pred) = outputs_pred_draft
  }
  if(only_positive){outputs_pred = (outputs_pred > 0)*outputs_pred}
  return(outputs_pred)
}
