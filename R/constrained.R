#' Constrained Precision Matrix
#'
#' @description Compute the maximum likelihood estimate of the precision matrix,
#' given a known graphical structure (i.e., an adjacency matrix).
#' This approach was originally described in
#' "The Elements of Statistical Learning"
#' \insertCite{@see pg. 631, @hastie2009elements}{GGMnonreg}.
#'
#' @param Sigma Covariance matrix
#'
#' @param adj Matrix with constraints. A zero indicates that element
#'            should be constrained to zero.
#'
#'
#' @return  A list containing the inverse covariance matrix
#' (i.e., precision matrix) and the covariance matrix.
#'
#' @references
#' \insertAllCited{}
#'
#' @note
#'
#' The algorithm is written in \code{c++}, and should scale to high dimensions.
#'
#' Note there are a variety of algorithms for this purpose. Simulation
#' studies indicated that this approach is both accurate and computationally
#' efficient \insertCite{@HFT therein, @emmert2019constrained}{GGMnonreg}

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
  # change to zeros
  # adj <- ifelse(adj == 1, 0, 1)

  # include diagonal!
  diag(adj) <- 1

  # call c++
  fit <- GGMncv::constrained(Sigma, adj)

  Theta <- round(fit$Theta, 3)
  Sigma <- round(fit$Sigma, 3)
  wadj <- -(cov2cor(Theta) - diag(ncol(adj)))
  returned_object <- list(Theta = Theta, Sigma = Sigma, wadj = wadj)

  return(returned_object)
}
