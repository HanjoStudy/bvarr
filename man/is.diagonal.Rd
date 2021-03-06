% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bvar_conjugate.R
\name{is.diagonal}
\alias{is.diagonal}
\title{check whether matrix is diagonal}
\usage{
is.diagonal(A, epsilon = 0)
}
\arguments{
\item{A}{symmetric matrix}

\item{epsilon}{tolerance, default is 0.}
}
\value{
logical, TRUE/FALSE, TRUE for diagonal matrices
}
\description{
check whether matrix is diagonal
}
\details{
check whether matrix is diagonal
}
\examples{
A <- matrix(c(2, 1, 1, 2), nrow = 2)
sym_inv(A)
}
