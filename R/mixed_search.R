#' Mixed Graphical Model: automated search
#'
#' @description Data mining to learn the graph.
#'
#' @param Y A data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes)
#'
#' @param data_type Vector of length \emph{p}. The type of data, with options of
#'                   "b" (binary), "p" (Poisson), and "g" (Gaussian).
#'
#' @param IC Character string. The desired information criterion. Options include
#'           \code{"AIC"} and \code{"BIC"} (default).
#'
#'
#' @return An object of class \code{mixed_search} including
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
#' @details Only backwards selection is currently implemented.
#'          Only an adjacency matrix is provided.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- ifelse( ptsd[,1:5] == 0, 0, 1)
#'
#' # search data (ising model)
#' fit <- mixed_search(Y, data_type = rep("b", 5))
#' }
#' @importFrom stats binomial coef glm step gaussian poisson
mixed_search <- function(Y, data_type = NULL, IC = "BIC"){

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

    if (data_type[x] == "b") {

      fit <- glm(y ~ ., data = dat, family = binomial())

    } else if (data_type[x] == "p") {

      fit <- glm(y ~ ., data = dat, family = poisson())

    } else if (data_type[x] == "g") {

      fit <- glm(y ~ ., data = dat, family = gaussian())

    } else {

      stop("'data_type' must be 'b', 'p' or 'g'")

    }

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
    wadj  = wadj,
    adj = adj,
    p = p,
    n = n,
    Y = Y
  )

  class(returned_object) <- c("ggmnonreg",
                              "mixed_search")

  return(returned_object)
}


