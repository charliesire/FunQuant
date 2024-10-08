% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Wavelet2D.R
\name{as.matrix.dwt.2d}
\alias{as.matrix.dwt.2d}
\title{Wavelet coefficients in matrix form from a \code{dwt.2d} object}
\usage{
\method{as.matrix}{dwt.2d}(x, nrows, ncols, ...)
}
\arguments{
\item{x}{an object of class \code{dwt.2d}.}

\item{nrows, ncols}{the number of rows and columns.}

\item{...}{additional arguments to be passed to or from methods.}
}
\value{
an object of class \code{matrix.dwt.2d} corresponding to a matrix with wavelet coefficients.
}
\description{
organizes the wavelet coefficients from a 2D DWT into a single matrix.
}
\author{
Tran Vi-vi Elodie PERRIN
}
\keyword{internal}
