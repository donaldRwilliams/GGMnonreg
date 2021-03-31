#' Gaussian graphical model: automated search
#'
#' @param x  A data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes) or
#'           a correlation matrix of dimensions \emph{p} by \emph{p}.
#'
#' @param IC Character string. The desired information criterion. Options include
#'           \code{"AIC"} and \code{"BIC"} (default).
#'
#' @param type Character string. Which search method should be used? The options included
#'             \code{"regression"} and \code{"approx_L0"}. See details.
#'
#'
#' @param method Character string. The desired subset selection method
#'               Options includes \code{"forward"} (default), \code{"backward"},
#'               and \code{"exhaustive"}.
#'
#' @param n Integer. Sample size. Required if a correlation matrix is provided.
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class \code{ggm_search}.
#'
#'
#' @note \code{type = "neighborhood_selection"} employs multiple regression to estimate
#' the graph (requires the data), whereas \code{type = "approx_L0"} directly estimates
#' the precision matrix (data or a correlation matrix are acceptable). If
#' data is provided and \code{type = "approx_L0"}, by default Pearson correlations are
#' used. For another correlation coefficient, provide the desired correlation matrix.
#'
#' \code{type = "approx_L0"} is a continuous approximation to (non-regularized)
#' best subset model selection. This is accomplished by using regularization, but
#' the penalty (approximately) mimics non-regularized estimation.
#'
#'
#' @details \code{type = "neighborhood_selection"} was described in
#' \insertCite{williams2019nonregularized;textual}{GGMnonreg}
#' and \code{type = "approx_L0"} was described in \insertCite{williams2020beyond;textual}{GGMnonreg}.
#' The penalty for \code{type = "approx_L0"} is called seamless L0 \insertCite{dicker2013variable}{GGMnonreg}
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- ptsd
#'
#' # search data
#' fit <- ggm_search(Y)
#' }
#'
#' @importFrom bestglm bestglm
#' @importFrom GGMncv ggmncv
ggm_search <- function(x, IC = "BIC",
                       type = "neighborhood_selection",
                       method = "forward",
                       n = NULL) {

  if(!IC %in% c("AIC", "BIC")) {
    stop("IC must be 'AIC' or 'BIC'")
  }

  if(type == "neighborhood_selection") {
    if (isSymmetric(as.matrix(x))) {
      stop("data required for 'neighborhood_selection'")
    }

    if (!method %in% c("forward", "backward", "exhaustive")) {
      stop("method must be 'forward', 'backward', or 'exhaustive'")
    }

    # change to X
    X <- x
    X <- scale(X, scale = T)
    n <- nrow(X)
    p <- ncol(X)

    betas <- matrix(0, p, p)

    colnames(betas) <- paste("X", 1:p, sep = "")

    colnames(X) <-  paste("X", 1:p, sep = "")

    pb <- txtProgressBar(min = 0, max = p, style = 3)

    estimates <- lapply(1:p, function(x) {
      est_i <- bestglm::bestglm(
        cbind.data.frame(X[,-x], X[, x]),
        method = method,
        IC = IC,
        intercept = F
      )$BestModel$coefficients

      setTxtProgressBar(pb, x)

      est_i
    })

    for (i in 1:p) {
      betas[i, names(estimates[[i]])] <- estimates[[i]]

    }

    wadj <- sqrt(betas * t(betas)) * sign(betas)
    adj <- ifelse(wadj == 0, 0, 1)


  } else if (type == "approx_L0") {

    p <- ncol(x)

    if (!isSymmetric(as.matrix(x))) {
      n <- nrow(x)
      R <- cor(x)
    } else {
      if (is.null(n)) {
        stop("sample size required when a correlation matrix is provided")
      }

      R <- x
    }

    if (IC == "AIC") {
      IC <- "aic"
    } else {
      IC <- "bic"
    }

    # seamless L0 penalty
    fit <-
      GGMncv::ggmncv(
        R,
        n = n,
        ic = IC,
        penalty = "selo",
        unreg = TRUE
      )

    wadj <- fit$P
    adj <- fit$adj

    }

  returned_object <- list(
    wadj = wadj,
    adj = adj,
    p = p,
    n = n,
    x = x
  )

  class(returned_object) <- c("ggmnonreg",
                              "ggm_search")
  return(returned_object)

}




