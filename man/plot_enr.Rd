% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enr.R
\name{plot_enr}
\alias{plot_enr}
\title{Plot \code{enr} Objects}
\usage{
plot_enr(x, iter = 1e+05, fill = "#009E73", alpha = 0.5, ...)
}
\arguments{
\item{x}{An object of class \code{enr}.}

\item{iter}{Integer. How many draws from the
Poisson-binomial distribution (defaults to 1,000)?}

\item{fill}{Which color to fill the density?}

\item{alpha}{Numeric (between 0 and 1). The transparency
for the density.}

\item{...}{Currently ignored}
}
\value{
An object of class \code{ggplot}
}
\description{
Plot the probability mass function for ENR.
}
\examples{
\donttest{
# correlations
cors <- cor(GGMnonreg::ptsd)

# inverse
inv <- solve(cors)

# partials
pcors <-  -cov2cor(inv)

# set values to zero
pcors <- ifelse(abs(pcors) < 0.05, 0, pcors )

est <- enr(net = pcors, n = 500, replications = 2)

# plot
plot_enr(est)
}

}
