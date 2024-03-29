#' Ising: automated search
#'
#' @description Data mining to learn the graph of binary variables with an Ising model
#'              \insertCite{lenz1920beitrvsge,ising1925beitrag}{GGMnonreg}.
#'
#' @param Y A data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes).
#'
#' @param IC Character string. The desired information criterion. Options include
#'           \code{"AIC"} and \code{"BIC"} (default).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE})?
#'
#' @references
#' \insertAllCited{}
#'
#' @details Currently only backwards selection is currently implemented.
#'
#' @note For an excellent overview of the Ising model see \insertCite{marsman2018introduction;textual}{GGMnonreg}.
#'
#' @return An object of class \code{ising_search} including:
#'
#' \itemize{
#'
#' \item{\strong{wadj}}: Weighted adjacency matrix, corresponding
#' to the partial correlation network.
#'
#' \item{\strong{adj}}: Adjacency matrix (detected effects).
#'
#' \item{\strong{pcors}}: Partial correlations.
#'
#' \item{\strong{n}}: Sample size.
#'
#' \item{\strong{p}}: Number of nodes.
#'
#' \item{\strong{Y}}: Data.
#'
#' }
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
ising_search <- function(Y, IC = "BIC", progress = TRUE) {
  if (!IC %in% c("AIC", "BIC")) {
    stop("IC must be 'AIC' or 'BIC'")
  }

  X <- Y

  n <- nrow(X)

  p <- ncol(X)

  betas <- matrix(0, p, p)

  colnames(betas) <- paste("X", 1:p, sep = "")

  colnames(X) <-  paste("X", 1:p, sep = "")

  if (progress) {
    pb <- txtProgressBar(min = 0,
                         max = p,
                         style = 3)
  }

  if (IC == "AIC") {
    k <- 2
  } else {
    k <- log(n)
  }

  estimates <- lapply(1:p, function(x) {
    dat <- cbind.data.frame(X[, -x], y = X[, x])

    fit <- glm(y ~ ., data = dat, family = binomial())

    est_i <-
      coef(step(
        fit,
        direction = "backward",
        k = k,
        trace = FALSE
      ))[-1]

    if (progress) {
      setTxtProgressBar(pb, x)
    }

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
