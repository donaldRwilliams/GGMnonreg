#' Constrained Precision Matrix
#'
#' @description Compute the maximum likelihood estimate, given certain elements are constrained to zero
#' (e.g., an adjacency matrix). This approach is described in \insertCite{hastie2009elements;textual}{GGMnonreg}.
#'
#' @param Sigma Covariance matrix
#'
#' @param adj Matrix with constraints. A zero indicates that element
#'            should be constrained to zero.
#'
#' @references
#' \insertAllCited{}
#'
#' @return  A list containing the inverse covariance matrix and the covariance matrix.
#'
#' @note The algorithm is written in \code{c++}.
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
