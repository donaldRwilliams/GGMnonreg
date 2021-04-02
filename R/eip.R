#' Edge Inclusion "Probability"
#'
#'  Compute the proportion of bootstrap samples that each relation was selected.
#'
#' @param Y The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes).
#'
#' @param method Character string. Which type of correlation coefficients
#'               to be computed. Options include \code{"pearson"} (default),
#'               \code{"kendall"}, \code{"spearman"}, and \code{"polychoric"}.
#'
#' @param B Integer. Number of bootstrap replicates (defaults to \code{1000}).
#'
#' @param progress. Logical. Should a progress bar be included (defaults to \code{TRUE})?
#'
#' @return An object of class \code{eip}.
#'
#' @details The order is the upper-triangular.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- ptsd
#'
#' # eip
#' fit_eip <- eip(Y, method = "spearman")
#'
#' # print
#' fit_eip
#' }
eip <- function(Y, method = "pearson", B = 1000, progress = TRUE) {

  n <- nrow(Y)
  p <- ncol(Y)

  if(progress){
  pb <- utils::txtProgressBar(max = B, style = 3)
  }
  eips <- t(sapply(1:B, function(x){

  eips_i <-  ggm_inference(Y[sample(1:n, n, replace = T),],
                  boot = FALSE, method = method)$adj[upper.tri(diag(p))]

  if(progress){
  setTxtProgressBar(pb, x)
  }

  eips_i
  }))

  returned_object <- list(eip =  data.frame( eip = colMeans( eips)))
  class(returned_object) <- c("ggmnonreg", "eip")
  return(returned_object)
}

print_eip <- function(x,...){
  print(x$eip,...)
}
