#' Non-parametric bootstrap
#' @param X Data matrix (dimensions n by p)
#' @param sims number of bootstrap samples
#' @param alpha type I error rate (correspondig to approximately 1 - specificity)
#'
#' @return mat_selected adjacency matrix
#' @return mat_mean mean of the bootstrap samples
#' @return boot_results bootstrap samples
#' @export
#'
#' @examples
#'
#' fit <- GGMboot(X, sims = 1000, alpha = 0.01)
#'
#' # adjacency matrix
#' fit_boot$mat_selected
GGMboot <- function(X, sims = 1000, alpha = 0.01){
  # Check data and remove missings:
  if (is.matrix(X)){
    X <- as.data.frame(X)
  }
  if (!is.data.frame(X)){
    stop("'X' is not a data frame")
  }

  # Listwise missingness removal:
  X <- na.omit(X)



  mat_selected <- mat_mean <- matrix(0, ncol(X), ncol(X))
  lw_bound <- alpha / 2
  up_bound <- 1 -   lw_bound
  boot_results <-  t(replicate(sims,  boot_inv_helper(X)$upper_tri, simplify = T))


  quantiles <- t(apply(boot_results, MARGIN = 2,  quantile, probs = c(lw_bound, up_bound)))

  means <-  t(apply(boot_results, MARGIN = 2,  mean))

  mat_selected[upper.tri(mat_selected)] <- mapply(quantile_helper, x = quantiles[,1], y = quantiles[,2])
  mat_selected[lower.tri(mat_selected)] <- t(mat_selected)[lower.tri(mat_selected)]

  mat_mean[upper.tri(mat_mean)] <-  means
  mat_mean[lower.tri(mat_mean)] <- t(mat_mean)[lower.tri(mat_mean)]

  list(mat_selected = mat_selected, mat_mean =  mat_mean, boot_results = boot_results)
}
