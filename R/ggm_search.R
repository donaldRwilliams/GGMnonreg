#' Gaussian graphical model: automated search
#'
#' @param Y  The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes).
#'
#' @param IC Character string. The desired information criterion. Options include
#'           \code{"AIC"} and \code{"BIC"} (default).
#'
#' @param method Character string. The desired subsect selection method
#'               Options includes \code{"forward"} (default), \code{"backward"},
#'               and \code{"exhaustive"}.
#'
#' @return An object of class \code{ggm_search}.
#'
#' @note Multiple regression is used to estimate the graph.
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
ggm_search <- function(Y, IC = "BIC", method = "forward") {

  if(!IC %in% c("AIC", "BIC")){
    stop("IC must be 'AIC' or 'BIC'")
  }

  if(!method %in% c("forward", "backward", "exhaustive")){
    stop("method must be 'forward', 'backward', or 'exhaustive'")
  }

  # change to X
  X <- Y
  X <- scale(X, scale = T)
  n <- nrow(X)
  p <- ncol(X)

  betas <- matrix(0, p, p)

  colnames(betas) <- paste("X", 1:p, sep = "")

  colnames(X) <-  paste("X", 1:p, sep = "")

  pb <- txtProgressBar(min = 0, max = p, style = 3)

  estimates <- lapply(1:p, function(x) {

    est_i <- bestglm::bestglm(
      cbind.data.frame(X[, -x], X[, x]),
      method = method,
      IC = IC,
      intercept = F
    )$BestModel$coefficients

    setTxtProgressBar(pb, x)

    est_i
    })

  for(i in 1:p){
    betas[i, names(estimates[[i]])] <- estimates[[i]]

    }

  wadj <- sqrt(betas * t(betas)) * sign(betas)
  adj <- ifelse(wadj == 0, 0, 1)

  returned_object <- list(
    wadj = wadj,
    adj = adj,
    p = p,
    n = n,
    Y = Y
  )
  return(returned_object)
}
