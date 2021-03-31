#' Ising: automated search
#'
#' @param Y A data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes)
#'
#' @param IC Character string. The desired information criterion. Options include
#'           \code{"AIC"} and \code{"BIC"} (default).
#'
#'
#' @details Only backwards selection is currently implemented.
#'
#' @return An object of class \code{ising_search}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- ifelse( ptsd[,1:5] == 0, 0, 1)
#'
#' # search data
#' fit <- ising_search(Y)
#' }
#' @importFrom stats binomial coef glm step
ising_search <- function(Y, IC = "BIC"){

  if (!IC %in% c("AIC", "BIC")) {
    stop("IC must be 'AIC' or 'BIC'")
  }

  X <- Y
  n <- nrow(X)
  p <- ncol(X)

  betas <- matrix(0, p, p)

  colnames(betas) <- paste("X", 1:p, sep = "")

  colnames(X) <-  paste("X", 1:p, sep = "")

  pb <- txtProgressBar(min = 0, max = p, style = 3)

  if (IC == "AIC") {
    k <- 2
  } else {
    k <- log(n)
  }

  estimates <- lapply(1:p, function(x) {
    dat <- cbind.data.frame(X[,-x], y = X[, x])
    fit <- glm(y ~ ., data = dat, family = binomial())
    est_i <-
      coef(step(
        fit,
        direction = "backward",
        k = k,
        trace = FALSE
      ))[-1]

    setTxtProgressBar(pb, x)

    est_i
  })

  for (i in 1:p) {
    betas[i, names(estimates[[i]])] <- estimates[[i]]

  }

  # taken from isingFit
  adj <- betas
  adj <- (adj != 0) * 1
  EN.weights <- adj * t(adj)
  EN.weights <- EN.weights * betas
  meanweights.opt <- (EN.weights + t(EN.weights)) / 2
  diag(meanweights.opt) <- 0
  wadj <- meanweights.opt
  adj <- ifelse(wadj == 0, 0, 1)


  returned_object <- list(
    wadj = wadj,
    adj = adj,
    p = p,
    n = n,
    Y = Y
  )

  class(returned_object) <- c("ggmnonreg",
                              "ising_search")

  return(returned_object)
}


