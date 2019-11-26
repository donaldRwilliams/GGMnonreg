#' Precision Matrix to Multiple Regression
#' @name coef.estimate
#' @description There is a direct correspondence between the covariance matrix and multiple regression. In the case of GGMs, it is possible
#' to estimate the edge set with multiple regression (i.e., neighborhood selection). In *GGMnonreg*, the precision matrix is first sampled from, and then
#' each draws is converted to the corresponding coefficients and error variances. This results in a posterior distribution. This function can be used
#' to perform multiple regression.
#'
#' @param object object of class \code{estimate} (analytic = F)
#' @param ci confidence interval used in the summary output
#' @return list of class \code{coef.estimate}:
#'
#' @examples
#'# data
#' X <- scale(GGMnonreg::ptsd)

#'# fit model
#' fit <- GGM_bootstrap(X)
#'
#'# summary for predicting the first variable
#' coef(fit, node = 1)
#' @export
coef.GGM_bootstrap <- function(object, node = 1, ci = 0.95){

  inv_to_coef <- list()
  p <- object$p
  lw_bound <- (1 - ci) / 2
  up_bound <- 1 -   lw_bound
  nodewise <-list()
  s <- 1:p
  for(i in 1:p){
    inv_to_coef[[i]] <- lapply(1:object$iter, function(x)
      (object$boot_results[[x]]$inv_cov[i,-i] / object$boot_results[[x]]$inv_cov[i,i]) * -1)
  }

  # bootstrap mean
  coef_mu <- lapply(1:p, function(x)  apply(do.call(rbind, inv_to_coef[[x]]), 2, mean) )

  # bootstrap SE
  coef_se <- lapply(1:p, function(x)  apply(do.call(rbind, inv_to_coef[[x]]), 2, sd) )

  # bootstrap confidence interval
  coef_qi <- lapply(1:p, function(x)  apply(do.call(rbind, inv_to_coef[[x]]), 2, quantile, c(lw_bound, up_bound)))

  for(i in 1:p){
    nodewise[[i]] <- data.frame(Node = s[- which(i == 1:20)],
                                Estimate = coef_mu[[i]],
                                Est.Error = coef_se[[i]],
                                CI = t(coef_qi[[i]]))
    names(nodewise)[[i]] <- paste("node", i)
  }

  returned_object <- list(nodewise = nodewise,
                          ci  = ci, node = node)

  class(returned_object) <- "coef.GGM_bootstrap"
  return(returned_object)
}

#' @name print.coef.GGM_bootstrap
#' @title  Print method for \code{coef.GGM_bootstrap} objects
#' @param object An object of class \code{coef.GGM_bootstrap}
#' @param ... currently ignored
#'
#' @seealso \code{\link{coef.GGM_bootstrap}}
#'
#' @export
print.coef.GGM_bootstrap <- function(x,..){

  node_printed <-  x$nodewise[[x$node]]
  colnames(node_printed) <- c("Node", "Estimate", "Est.Error",
                                           paste(c("lb.", "ub."), gsub("*0.","", x$ci), "%", sep = ""))

  print(round(node_printed, 3), row.names = FALSE)
}
