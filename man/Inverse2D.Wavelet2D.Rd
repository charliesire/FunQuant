% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Wavelet2D.R
\name{Inverse2D.Wavelet2D}
\alias{Inverse2D.Wavelet2D}
\title{Multiple two-dimensional inverse wavelet transform}
\usage{
\method{Inverse2D}{Wavelet2D}(object, ...)
}
\arguments{
\item{object}{an object of class \code{\link{Wavelet2D}}.}

\item{...}{other arguments}
}
\value{
a three dimensional array.
The first two dimensions correspond to maps dimensions.
The third dimension corresponds to the size of the data set
}
\description{
performs inverse wavelet transform of each set of wavelet coefficients from a data set.
}
\seealso{
\code{Wavelet2D}
}
\author{
Tran Vi-vi Elodie PERRIN
}
