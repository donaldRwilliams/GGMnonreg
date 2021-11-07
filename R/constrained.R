#' Precision Matrix with Known Graph
#'
#' @description Compute the maximum likelihood estimate of the precision matrix,
#' given a known graphical structure (i.e., an adjacency matrix).
#' This approach was originally described in
#' "The Elements of Statistical Learning"
#' \insertCite{@see pg. 631, @hastie2009elements}{GGMnonreg}.
#'
#' @param Sigma Covariance matrix
#'
#' @param adj An adjacency matrix that encodes the constraints,
#'            where a zero indicates that element should be zero.
#'
#'
#' @return  A list containing the following:
#'
#' \itemize{
#'
#' \item{\strong{Theta}}: Inverse of the covariance matrix
#' (precision matrix), that encodes the conditional
#' (in)dependence structure.
#'
#' \item{\strong{Sigma}}: Covariance matrix.
#'
#' \item{\strong{wadj}}: Weighted adjacency matrix, corresponding
#' to the partial correlation network.
#'
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @note
#' The algorithm is written in \code{c++}, and should scale to high dimensions.
#'
#' Note there are a variety of algorithms for this purpose. Simulation
#' studies indicated that this approach is both accurate and computationally
#' efficient \insertCite{@HFT therein, @emmert2019constrained}{GGMnonreg}
#'
#'
#'
#' @examples
#' # data
#' Y <- ptsd
#'
#' # estimate graph
#' fit <- ggm_inference(Y, boot = FALSE)
#'
#' # constrain to zero
#' constrained_graph <- constrained(cor(Y), fit$adj)
#'
#' @importFrom GGMncv constrained
#' @export
constrained <- function(Sigma, adj){

  # include diagonal
  # note: an error otherwise
  diag(adj) <- 1

  # call c++
  fit <- GGMncv::constrained(Sigma, adj)

  # precision matrix
  Theta <- round(fit$Theta, 3)

  # covariance matrix
  Sigma <- round(fit$Sigma, 3)

  # weighted adjacency matrix
  # partial correlation network
  wadj <- -(cov2cor(Theta) - diag(ncol(adj)))

  returned_object <- list(Theta = Theta,
                          Sigma = Sigma,
                          wadj = wadj)

  return(returned_object)
}
