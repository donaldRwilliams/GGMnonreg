#' Compare Gaussian Graphical Models
#'
#' @param Yg1 The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes) for
#'            group one.
#'
#' @param Yg2 The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes) for
#'            group two.
#'
#' @param alpha The desired significance level (defaults to \code{0.05}). Note that
#'              1 - alpha corresponds to specificity.
#'
#' @return An object of class \code{ggm_compare}
#' @export
#'
#' @examples
#' \donttest{
#' # data
#'
#' Yg1 <- na.omit(subset(bfi, gender == 1)[,1:10])
#' Yg2 <- na.omit(subset(bfi, gender == 2)[,1:10])
#'
#' # compare relations
#' fit <- ggm_compare(Yg1, Yg2)
#'}
ggm_compare <- function(Yg1, Yg2, alpha = 0.05){

  fit1 <- ggm_inference(Yg1)
  fit2 <- ggm_inference(Yg2)

  diff <- fit1$boot_samps - fit2$boot_samps

  ci_lower <- alpha / 2
  ci_upper <- 1 - ci_lower

  p <- ncol(Yg1)

  adj <- wadj <- matrix(0, nrow = p, ncol = p)

  cis <- t(apply(diff, 2, quantile, probs = c(ci_lower, ci_upper)))

  adj[upper.tri(adj)] <- ifelse(cis[, 1] < 0 & cis[, 2] > 0, 0, 1)
  adj <- symm_mat(adj)

  wadj[upper.tri(diag(p))] <- colMeans(diff)
  wadj <- symm_mat(wadj) * adj

  returned_object <- list(adj = adj,
                          wadj = wadj,
                          cis = cis)

  class(returned_object) <- c("ggmnonreg",
                              "ggm_compare")

  return(returned_object)
}
