% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/std_proba.R
\name{std_proba}
\alias{std_proba}
\title{Computation of the IS coefficients of variation of the membership probability for different set of prototypes.}
\usage{
std_proba(
  data = NULL,
  prototypes_list,
  density_ratio = rep(1, dim(data)[length(dim(data))]),
  distance_func = function(A1, A2) {
     return(sqrt(sum((A1 - A2)^2)))
 },
  cells,
  cell_numbers = NULL,
  nv = NULL,
  outputs_function = NULL,
  inputs = NULL,
  batch_size = NULL,
  return_cell_numbers = FALSE,
  bootstrap = NULL
)
}
\arguments{
\item{data}{The data that needs to be quantized. An array of any dimension is expected, the kth element must be selected with data[,..,k].
Useful only if cell_numbers == NULL and outputs_function == NULL.}

\item{prototypes_list}{A list of set of prototypes on which we want to evaluate the IS centroid standard deviation. Each element is a list of prototypes.}

\item{density_ratio}{density_ratio indicates the weight fX/g of each data element.}

\item{distance_func}{A function computing a distance between two data elements.}

\item{cells}{The Voronoï cell numbers that we are investigating}

\item{cell_numbers}{An optional list providing for each set of prototypes the voronoi cell number of every data element.}

\item{nv}{The size of the sample for which we want to estimate the IS coefficient of variation of the membership probability}

\item{outputs_function}{Function to compute the outputs from the inputs.}

\item{inputs}{If data and cell_numbers are not provided, inputs is a dataframe of the inputs at which the outputs will be computed.}

\item{batch_size}{If the computation is by batch, the number of elements by batch. Optional.}

\item{return_cell_numbers}{Boolean indicating whether the cell_numbers must be provided or not.}

\item{bootstrap}{Integer indicating the number of bootstramp samples to generate. If NULL, no boostrap is performed. Default is NULL.}
}
\value{
A list of IS coefficients of variation of the membership probability obtained for each set of prototypes.
}
\description{
Computation of the IS coefficients of variation of the membership probability for different set of prototypes.
}
\examples{
func2D <- function(X){
Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
n<-nrow(X)
Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2
-10*X[i,])**2)/(60*X[i,]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
Ymaps<- array(unlist(Y),dim=c(20,20,n))
return(Ymaps)
}
design = data.frame(X = seq(-1,1,l= 50))
data = func2D(design)
prototypes_list = list(lapply(c(1,3,10,14,18), function(i){data[,,i]}))
density_ratio = rep(1, 50)
distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
list_std_proba = std_proba(data = data, prototypes_list = prototypes_list,
 density_ratio = density_ratio, distance_func = distance_func,
 cells = 1:length(prototypes_list[[1]]), nv = 50)
}
