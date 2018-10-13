#' Maximum likelihood precision and covariance matrices
#'
#' @param X Data Matrix (dimensions n by p)
#'
#' @return inv_cov_mle MLE precision matrix
#' @return cov_mle MLE covariance matrix
#' @export
#'
#' @examples
mle <- function(X){
      inv_cov_mle <- svd_inv_helper(X)
      cov_mle <- solve(inv_cov_mle)
      list(inv_cov_mle = inv_cov_mle, cov_mle = cov_mle)
}
