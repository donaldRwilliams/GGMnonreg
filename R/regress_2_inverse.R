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

  # remove mean structure
  X <- scale(x, scale = F)

  # number of variables
  p <- ncol(x)

  # matrix storage
  mat <- matrix(0, p, p)

  for(i in 1:p){

    # fit regression model
    coef <- solve(t(X[,-i]) %*%  X[,-i]) %*% t(X[,-i]) %*% X[,i]

    # fitted values
    fitted <- X[,-i] %*% coef

    # residual variance
    res_var <- as.numeric(var(X[,i] - fitted))

    # off-diagonal elements
    off_diagonal <- (coef / res_var) * -1

    # diagonal elements
    mat[i,i] <- diagonal <- 1 / res_var

    # off diagonal elements
    mat[i,-i] <- off_diagonal
}

list(inverse_cov = mat, cov_matrix = solve(mat))


}



