% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constrained.R
\name{constrained}
\alias{constrained}
\title{Constrained Precision Matrix}
\usage{
constrained(Sigma, adj)
}
\arguments{
\item{Sigma}{Covariance matrix}

\item{adj}{Matrix with constraints. A zero indicates that element
should be constrained to zero.}
}
\value{
A list containing the inverse covariance matrix and the covariance matrix.
}
\description{
Compute the maximum likelihood estimate, given certain elements are constrained to zero
(e.g., an adjacency matrix). This approach is described in \insertCite{hastie2009elements;textual}{GGMnonreg}.
}
\note{
The algorithm is written in \code{c++}.
}
\examples{
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# constrain to zero
constrained_graph <- constrained(cor(Y), fit$adj)

}
\references{
\insertAllCited{}
}
