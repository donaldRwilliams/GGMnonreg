#' Network Predictability (R2)
#'
#' @param x An object of class \code{ggm_inference}
#'
#' @param ci Numeric. The confidence interval to be computed (defaults to \code{0.95}).
#'
#' @return An object of class \code{predictability}, including a matrix of R2.
#'
#' @references
#' \insertAllCited{}
#'
#' @note Predictability is variance explained for each node in the
#'       network \insertCite{Haslbeck2018}{GGMnonreg}.
#'
#' @export
#'
#' @examples
#' # data
#' Y <- ptsd
#'
#' # estimate graph
#' fit <- ggm_inference(Y, boot = FALSE)
#'
#' # predictability
#' r2 <- predictability(fit)
#'
#' # print
#' r2
#'
#' @importFrom corpcor pcor2cor
predictability <- function(x, ci = 0.95) {

  if (!is(x, "ggm_inference")) {
    stop("must be a 'ggm_inference' object")
  }

  p <- x$p
  n <- x$n

  # scale y
  ynew <- scale(x$Y)

  # pcor to cor to inverse
  diag(x$pcors) <- 1
  cors <- corpcor::pcor2cor(x$pcors)
  inv <- solve(cors)

  crt <- qnorm((1 - ci) / 2, lower.tail = F)

  r2_ls <- list()

  for(i in 1:p){

    # var of y 1 by construction
    r2 <- var(ynew[,-i] %*%  inv[i,-i] / inv[i, i])

    num <- 4 * r2 * (1 - r2)^2 * (n - p - 1)^2

    den <- (n ^ 2 - 1) * (n + 3)

    # standard error
    se <- sqrt(num / den)

    # ci lower and upper
    upper <- r2 + se * crt
    lower <- r2 - se * crt


   r2_ls[[i]] <- data.frame(
      Estimate = r2,
      Est.Error =  se,
      Ci.lb = lower,
      Ci.ub = upper
    )
  }

  returned_object <- list(r2 = do.call(rbind.data.frame,  r2_ls))

  class(returned_object) <- c("ggmnonreg",
                              "predictability")
  return(returned_object)

}

print_r2 <- function(x, ...){
  print(x,...)
}

