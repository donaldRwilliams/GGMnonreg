% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{ptsd_cor4}
\alias{ptsd_cor4}
\title{Cor: Post-Traumatic Stress Disorder  (Sample # 4)}
\format{
A correlation matrix with 16 variables
}
\description{
A correlation matrix that includes 16 variables. The correlation matrix
was estimated from 965 individuals \insertCite{fried2018replicability}{GGMnonreg}.
}
\details{
\itemize{
  \item Intrusive Thoughts
  \item Nightmares
  \item Flashbacks
  \item Physiological/psychological reactivity
  \item Avoidance of thoughts
  \item Avoidance of situations
  \item Amnesia
  \item Disinterest in activities
  \item Feeling detached
  \item Emotional numbing
  \item Foreshortened future
  \item Sleep problems
  \item Irritability
  \item Concentration problems
  \item Hypervigilance
  \item Startle response
}
}
\examples{
data(ptsd_cor4)
Y <- MASS::mvrnorm(n = 965,
                   mu = rep(0, 16),
                   Sigma = ptsd_cor4,
                   empirical = TRUE)

}
\references{
\insertAllCited{}
}
\keyword{datasets}
