#' regression to inverse covariance matrix
#'
#' @param x data matrix (dimensions n by p)
#'
#' @return inverse-cov Inverse of the covariance matrix (i.e., the precision matrix)
#' @return cov_matrix Covariance matrix
#' @export
#'
#' @examples
#' p = 10
#' n = 100
#' x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#'
#' estimate <- regress_2_inverse(x)
#'
#' all.equal(solve(cov(x)), estimate$inverse_cov)
#'
#'
regress_2_inverse <- function(x){
  X <- scale(x, scale = F)
  p <- ncol(x)

  mat <- matrix(0, p, p)

  for(i in 1:p){

    fit <-  lm(X[,i] ~ X[,-i])
    mat[i,-i] <-   (- 1 * coefficients(fit)[-1]) / var(residuals(fit))
    mat[i,i] <- 1 / var(residuals(fit))


    }
list(inverse_cov = mat, cov_matrix = solve(mat))


}



