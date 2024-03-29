% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictability.R
\name{predictability}
\alias{predictability}
\title{Network Predictability (R2)}
\usage{
predictability(x, ci = 0.95)
}
\arguments{
\item{x}{An object of class \code{ggm_inference}}

\item{ci}{Numeric. The confidence interval to be computed (defaults to \code{0.95}).}
}
\value{
An object of class \code{predictability}, including a matrix of R2
for each node.
}
\description{
Network Predictability (R2)
}
\note{
Predictability is variance explained for each node in the
      network \insertCite{Haslbeck2018}{GGMnonreg}.
}
\examples{
# data
Y <- ptsd

# estimate graph
fit <- ggm_inference(Y, boot = FALSE)

# predictability
r2 <- predictability(fit)

# print
r2

}
\references{
\insertAllCited{}
}
