#' Title
#'
#' @param X
#' @param alpha
#' @param FDR
#'
#' @return
#' @export
#'
#' @examples
GGM_fisher_z <- function(X, alpha, FDR = NULL){


  n <- nrow(X)
  p <- ncol(X)
  I = diag(1,p)
  ## compute maximum likelihood estimator for covariance matrix
  mle_cov <- crossprod(scale(X, scale = F)) / n

  ## compute maximum likelihood estimator of precision matrix
  ## (inverse covariance matrix)
  mle_inv <- solve(mle_cov)

  ## standardize and revese sign = partial correaltions
  pc  <-  as.matrix(cov2cor(mle_inv))  * - 1

  mle_parcors <- ci_par_cor(alpha = alpha, par_cors = pc, n = n, s = p - 1)
  mle_inv <- mle_parcors$sig_mat * mle_inv



  if(isTRUE(FDR)){
    mat_fdr <- matrix(0, p, p)
    pvalues <- mle_parcors$pmat[upper.tri(mle_parcors$pmat)]
    mat_fdr[upper.tri(mat_fdr)] <- ifelse(stats::p.adjust(pvalues, method = "fdr") < alpha, 1, 0)
    mat_fdr <- as.matrix(Matrix::forceSymmetric(mat_fdr))

    list(mle_parcors = mle_parcors, mle_inv = mle_inv , mat_fdr = mat_fdr)
  } else{

    list(mle_parcors = mle_parcors, mle_inv = mle_inv)

  }




}
