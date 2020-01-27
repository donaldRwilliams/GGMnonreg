#' @importFrom stats na.omit quantile cov2cor pnorm qnorm sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
power_z <- function(r, n, c = 0,
                    type = "pearson",
                    compare = FALSE,
                    alpha = 0.05){

  # abs r
  r <- abs(r)

  # fisher z
  z <- atanh(r)

  if(type == "pearson"){
    # variance of z
    var_z <- 1 / (n - c - 3)

  } else if(type == "spearman"){

    var_z <- (1 + r^2/2) / (n - c - 3)

  } else{

    stop("invalid type (must be pearson or spearman)")
  }


  # differnece ?
  if(compare == TRUE){

    var_z <- var_z * 2

  }

  # z score
  z_score <- z/ sqrt(var_z)

  # quantile
  q <- stats::qnorm(1 - alpha / 2)

  # power
  1 - stats::pnorm(q - z_score)

}


or_helper <- function(x, y){
  ifelse(x * y == 0, max(abs(x), abs(y)), sqrt(abs(x) * abs(y)))
}

# fisher z to r
z2r <- function (z) {
  (exp(2 * z) - 1)/(1 + exp(2 * z))
}

fisher_z <- function(rho){
  .5 * log(( 1 + rho )/ ( 1 - rho ))
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

ci_par_cor <- function(alpha, par_cors, s, n) {

  # n: sample size
  # s: p - 1 (controlled for)
  # alpha: confidence level
  # par_cors: partial correlations
  mat <- pmat <- matrix(0,nrow =  s + 2, ncol = s + 2)
  CI_ls <- list()
  par_cor <- par_cors[upper.tri(par_cors)]
  cov <- list()
  pvalues <-list()
  for(i in 1:length(par_cor)){
    # crtiical value
    z_crit <- qnorm(1 - alpha/2)
    # standard error
    se <- sqrt(1/((n - s - 3)))
    # z transformation
    z <- log((1 + par_cor[i])/(1 - par_cor[i]))/2
    # z lower bound
    Z_L <- z - z_crit * se
    # Z upper bound

    Z_U <- z + z_crit * se
    rho_L <- (exp(2*Z_L) - 1)/(exp(2*Z_L) + 1)
    rho_U <- (exp(2*Z_U) - 1)/(exp(2*Z_U) + 1)
    CI <- c(rho_L, rho_U)
    CI_ls[[i]] <- CI
    cov[[i]] <- ifelse(CI[1] < 0 & CI[2] > 0, 0, 1)
    pvalues[[i]] <- 2 * pnorm(abs(z / se), lower.tail = F)
  }

  c_dat <- do.call(rbind.data.frame, CI_ls)
  colnames(c_dat) <- c("low", "up")
  c_dat$pcor <- unlist(par_cor)
  diag(mat) <- 1
  mat[upper.tri(mat)] <- unlist(cov)
  pmat[upper.tri(pmat)] <- unlist(pvalues)

  mat <- as.matrix(Matrix::forceSymmetric(mat))
  pmat <- as.matrix(Matrix::forceSymmetric(pmat))

  list(sig_mat = mat, par_cors = par_cors, par_sig = mat * par_cors,
       cis = c_dat, cov_prob = unlist(cov), pmat = pmat)
}

performance <- function(Estimate, True){

  True <- as.matrix(True)
  Estimate <- as.matrix(Estimate)

  # True Negative
  TN <- ifelse(True[upper.tri(True)] == 0 & Estimate[upper.tri(Estimate)] == 0, 1, 0); TN <- sum(TN)
  # False Positive
  FP <- ifelse(True[upper.tri(True)] == 0 & Estimate[upper.tri(Estimate)] != 0, 1, 0); FP <- sum(FP)
  # True Positive
  TP <- ifelse(True[upper.tri(True)] != 0 & Estimate[upper.tri(Estimate)] != 0, 1, 0); TP <- sum(TP)
  # False Negatives
  FN <- ifelse(True[upper.tri(True)] != 0 & Estimate[upper.tri(Estimate)] == 0, 1, 0); FN <- sum(FN)

  Specificity <- TN/(TN + FP)
  Sensitivity <- TP/(TP + FN)
  Precision <- TP/(TP + FP)

  Recall <- TP / (TP + FN)

  F1_score <- 2 * ((Precision * Recall) / (Precision + Recall))

  MCC <- (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  results <- c(Specificity, Sensitivity, Precision, Recall,  F1_score, MCC)
  results_name <- c("Specificity", "Sensitivity", "Precision", "Recall",  "F1_score", "MCC")
  results <- cbind.data.frame(measure = results_name, score = results)
  list(results = results)


}
