or_helper <- function(x, y){
  ifelse(x * y == 0, max(abs(x), abs(y)), sqrt(abs(x) * abs(y)))
}

svd_inv_helper <- function(X){
  n = nrow(X)
  X <- scale(X, scale = F)
  mle_cov <- t(X) %*% X / n
  svd_decomp <- svd(mle_cov)
  inv_cov <- svd_decomp$v  %*% (1/ svd_decomp$d * t(svd_decomp$u))
  return(inv_cov)
}

boot_inv_helper <- function(X){
  n <- nrow(X)
  boot_sample <- X[sample(1:n, size = n, replace = T),  ]
  inv_cov <- svd_inv_helper(boot_sample)
  pcor <- cov2cor(inv_cov) * -1
  diag(pcor) <- 1
  list(inv_cov = inv_cov, pcor = pcor, upper_tri = pcor[upper.tri(pcor)])

}

quantile_helper <- function(x, y){
  ifelse(x < 0 &  y > 0, 0, 1)
}

