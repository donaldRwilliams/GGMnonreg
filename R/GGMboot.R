#' Title
#'
#' @param X Data matrix (n times p)
#' @param sims number of bootstrap samples
#' @param alpha type I error rate (correspondig to 1 - specificity)
#'
#' @return
#' @export
#'
#' @examples
GGMboot <- function(X, sims, alpha){
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

  list(mat_selected = mat_selected, mat_mean =  mat_mean)
}
