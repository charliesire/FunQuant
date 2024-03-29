% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_denom_centroid.R
\name{estim_denom_centroid}
\alias{estim_denom_centroid}
\title{Title Compute the estimator which is the denominator of the centroid estimation}
\usage{
estim_denom_centroid(density_ratio, cell_numbers, cell, bias = 0)
}
\arguments{
\item{density_ratio}{density_ratio indicates the weight fX/g of each data element}

\item{cell_numbers}{The voronoi cell number of every data element.}

\item{cell}{The cell number of the computed centroid}

\item{bias}{A number indicating the bias that came out when computing the importance sampling estimators of the membership probabilities of the Voronoi cell. Default is 0.}
}
\value{
A real number which is the denominator of the centroid estimation
}
\description{
Title Compute the estimator which is the denominator of the centroid estimation
}
\examples{
density_ratio = rep(1,20)
cell_numbers = c(1,3,2,1,2,1,1,2,3,3,2,2,2,2,2,3,1,1,3,3)
cell = 3
estim_denom_centroid(density_ratio = density_ratio, cell_numbers = cell_numbers, cell = cell)
}
