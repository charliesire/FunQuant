% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_probas.R
\name{get_probas}
\alias{get_probas}
\title{Computing the probability masses of each voronoi cells}
\usage{
get_probas(
  density_ratio,
  method_IS = "unique",
  cell_numbers = NULL,
  data = NULL,
  prototypes = NULL,
  distance_func = function(A1, A2) {
     return(sqrt(sum((A1 - A2)^2)))
 },
  cells = NULL,
  bias = NULL
)
}
\arguments{
\item{density_ratio}{density_ratio indicates the weight fX/g of each data element}

\item{method_IS}{The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell. Default is "unique".}

\item{cell_numbers}{The voronoi cell number of every data element. If NULL, then data, prototypes and distance_func must be provided.}

\item{data}{The data that needs to be quantized. Useful only if cell_numbers == NULL.}

\item{prototypes}{A set of prototypes. Useful only if cell_numbers == NULL. If NULL, "cells" must be provided.}

\item{distance_func}{A function computing a distance between two data elements. Useful only if cell_numbers == NULL.}

\item{cells}{The cell numbers that are investigated}

\item{bias}{A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoi cell.}
}
\value{
A vector providing the probability masses of each Voronoï cell.
}
\description{
Computing the probability masses of each voronoi cells
}
\examples{
prototypes = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3)))
data = array(runif(9*20)*20, dim = c(3,3,20))
distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
density_ratio = rep(1,20)
get_probas(density_ratio = density_ratio, data = data,
prototypes = prototypes, distance_func = distance_func)
}
