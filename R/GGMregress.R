#' Network with Multiple Regression (Ordinary Least Squares)
#'
#' @param X Data Matrix (dimensions n by p)
#' @param IC Information criteria (AIC or BIC)
#' @param method Subsect selection method
#'
#' @return pcor_and    partial correlations matrix with "and-rule"
#' @return pcor_or     partial correlations matrix with "or-rule"
#' @return adj_and     adjacency matrix with "and-rule"
#' @return adj_or      adjacency matrix with "or-rule"
#' @export
#'
#' @examples
#'
#' fit <- GGM_regression(X, IC = "BIC", method = "forward")
#'
#' partial correlation matrix with "and-rule"
#' fit$pcor_and
#'
#' partial correlation matrix with "and-rule"
#' fit$adj_and
#'
GGM_regression <- function(X, IC, method){
  # scale data

  if( IC != "AIC" && IC != "BIC" ){
    stop("IC must be AIC or BIC")
  }

  if( method != "forward" && method != "backward"  && method != "exhaustive"){
    stop("method must be foward, backward, or exhaustive")
  }

  X <- scale(X, scale = T)
  n <- nrow(X)
  p <- ncol(X)

  mat1 <- mat2 <- mat_or <- matrix(0, p, p)
  colnames(mat1) <- paste("X", 1:p, sep = "")
  colnames(X) <- 1:p

  test <- lapply(1:p, function(x) bestglm::bestglm(cbind.data.frame(X[,-x], X[,x]),
                                  method = method, IC = IC, intercept = F)$BestModel$coefficients)

  for(i in 1:p){
    mat1[i,names(test[[i]])] <- test[[i]]

  }

  mat2[lower.tri(mat2)] <- sqrt(abs(mat1[lower.tri(mat1)]) * abs(t(mat1)[lower.tri(mat1)]))
  mat2[upper.tri(mat2)] <- t(mat2)[upper.tri(mat2)]

  mat_and <- sign(mat1) * mat2


  up <- mat1[upper.tri(mat1)]
  low <- t(mat1)[upper.tri(mat1)]

  sign_or <- sign(up)

  mat_or[upper.tri(mat_or)] <- mapply(or_helper,low, up) * sign_or
  mat_or[lower.tri(mat_or)] <- t(mat_or)[lower.tri(mat_or)]


  adj_and <- ifelse(mat_and == 0, 0, 1)
  adj_or  <- ifelse(mat_or == 0, 0, 1)

  list(pcor_and = mat_and,
       pcor_or = mat_or,
       adj_or = adj_or,
       adj_and = adj_and)

}


