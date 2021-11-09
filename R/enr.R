#' Expected Network Replicability
#'
#' @description Investigate network replicability for any kind of
#'              partial correlation, assuming there is an analytic
#'              solution for the standard error (e.g., Pearson's or Spearman's).
#'
#' @param net True network of dimensions \emph{p} by \emph{p}.
#'
#' @param n Integer. The samples size, assumed equal in the replication
#'          attempts.
#'
#' @param alpha The desired significance level (defaults to \code{0.05}). Note that
#'              1 - alpha corresponds to specificity.
#'
#' @param replications Integer. The desired number of replications.
#'
#' @param type Character string. Which type of correlation coefficients
#'             to be computed. Options include \code{"pearson"} (default)
#'             and \code{"spearman"}.
#'
#' @references
#' \insertAllCited{}
#'
#' @return An list of class \code{enr} including the following:
#'
#' \itemize{
#'
#' \item{\strong{ave_power}}: Average power.
#'
#' \item{\strong{cdf}}: cumulative distribution function.
#'
#' \item{\strong{p_s}}: Power for each edge, or the probability
#' of success for a given trial.
#'
#' \item{\strong{p}}: Number of nodes.
#'
#' \item{\strong{n_nonzero}}: Number of edges.
#'
#' \item{\strong{n}}: Sample size.
#'
#' \item{\strong{replication}}: Replication attempts.
#'
#' \item{\strong{var_pwr}}: Variance of power.
#'
#' \item{\strong{type}}: Type of correlation coefficient.
#'
#' }
#'
#' @note This method was introduced in
#' \insertCite{williams2020learning;textual}{GGMnonreg}.
#'
#' The basic idea is to determine the replicability of edges in a
#' partial correlation network. This requires defining the true
#' network, which can include edges of various sizes, and then
#' solving for the proportion of edges that are expected
#' to be replicated (e.g. in two, three, or four replication attempt).
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' # (1) define partial correlation network
#'
#' # correlations from ptsd symptoms
#' cors <- cor(GGMnonreg::ptsd)
#'
#' # inverse
#' inv <- solve(cors)
#'
#' # partials
#' pcors <-  -cov2cor(inv)
#'
#' # set values to zero
#' # (this is the partial correlation network)
#' pcors <- ifelse(abs(pcors) < 0.05, 0, pcors)
#'
#'
#' # compute ENR in two replication attempts
#' fit_enr <- enr(net = pcors,
#'                n = 500,
#'                replications = 2)
#'
#'
#' # intuition for the method:
#' # The above did not require simulation, and here I use simulation
#' # for the same purpose.
#'
#' # location of edges
#' # (where the edges are located in the network)
#' index <- which(pcors[upper.tri(diag(20))] != 0)
#'
#' # convert network a into correlation matrix
#' # (this is needed to simulate data)
#' diag(pcors) <- 1
#' cors_new <- corpcor::pcor2cor(pcors)
#'
#' # replicated edges
#' # (store the number of edges that were replicated)
#' R <- NA
#'
#' # simulate how many edges replicate in two attempts
#' # (increase 100 to, say, 5,000)
#' for(i in 1:100){
#'
#'   # two replications
#'   Y1 <- MASS::mvrnorm(500, rep(0, 20), cors_new)
#'   Y2 <- MASS::mvrnorm(500, rep(0, 20), cors_new)
#'
#'   # estimate network 1
#'   fit1 <- ggm_inference(Y1, boot = FALSE)
#'
#'   # estimate network 2
#'   fit2 <- ggm_inference(Y2, boot = FALSE)
#'
#'   # number of replicated edges (detected in both networks)
#'   R[i] <- sum(
#'     rowSums(
#'       cbind(fit1$adj[upper.tri(diag(20))][index],
#'             fit2$adj[upper.tri(diag(20))][index])
#'     ) == 2)
#' }
#'
#'
#' # combine simulation and analytic
#' cbind.data.frame(
#'   data.frame(simulation = sapply(seq(0, 0.9, 0.1), function(x) {
#'     mean(R > round(length(index) * x) )
#'   })),
#'   data.frame(analytic = round(fit_enr$cdf, 3))
#' )
#'
#' # now compare simulation to the analytic solution
#' # average replicability (simulation)
#' mean(R / length(index))
#'
#' # average replicability (analytic)
#' fit_enr$ave_pwr
#' }
#'
#' @importFrom poibin ppoibin
enr <- function(net, n,
                alpha = 0.05,
                replications = 2,
                type = "pearson"){

  # variables
  p <- ncol(net)

  # control variables
  c <- p - 2

  # off diagonals
  pcs <- net[upper.tri(net)]

  # which are edges
  which_nonzero <- which(pcs != 0)

  # how many non zeroes
  n_nonzero <- length(which_nonzero)

  qs <- round(seq(0, 0.9, 0.1) * length(which_nonzero))
    # quantile(1:length(which_nonzero),
                        # probs = seq(0, 0.9, 0.1)))

  # probability of success
  p_s <- power_z(pcs[which_nonzero],
                 n = n,
                 c = c,
                 type = type,
                 alpha = alpha)^replications

  # average power
  ave_pwr <- mean(p_s)

  #var
  var_pwr <- sum((1 - p_s) * p_s)

  # CDF greater
  cdf <- 1 - poibin::ppoibin(qs, p_s)

  returned_object <- list(ave_pwr = ave_pwr,
                          cdf = cdf,
                          p_s = p_s,
                          p = p,
                          n_nonzero = n_nonzero,
                          n = n,
                          replications = replications,
                          var_pwr = var_pwr,
                          type = type)

  class(returned_object) <- c("ggmnonreg", "enr")
  returned_object
}

print_enr <- function(x,...){

  cat("Average Replicability:", round(x$ave_pwr, 2), "\n")

  cat(
    "Average Number of Edges:",
    round(round(x$ave_pwr, 2) * x$n_nonzero),
    paste0("(SD = ", round(sqrt(x$var_pwr), 2), ")"),
    "\n\n"
  )

  cat("----\n\n")

  cat("Cumulative Probability:" , "\n\n")

  dat <- data.frame(prop = seq(0, .90, .10),
                    round(x$n_nonzero *  seq(0, 0.9, 0.1)),
                    prob = round(x$cdf, 2))

  colnames(dat) <- c("prop.edges", "edges", "Pr(R > prop.edges)")

  print(dat, row.names = F, right = T)

  cat("----\n")

  cat(
    paste0(
      "Pr(R > prop.edges):\n",
      "probability of replicating more than the\n",
      "correpsonding proportion (and number) of edges"
    )
  )
  }


#' @title Plot \code{enr} Objects
#'
#' @description Plot the probability mass function for ENR.
#'
#' @param x An object of class \code{enr}.
#'
#' @param iter Integer. How many draws from the
#'             Poisson-binomial distribution (defaults to 1,000)?
#'
#' @param fill Which color to fill the density?
#'
#' @param alpha Numeric (between 0 and 1). The transparency
#'              for the density.
#'
#' @param ... Currently ignored
#'
#' @return An object of class \code{ggplot}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # correlations
#' cors <- cor(GGMnonreg::ptsd)
#'
#' # inverse
#' inv <- solve(cors)
#'
#' # partials
#' pcors <-  -cov2cor(inv)
#'
#' # set values to zero
#' pcors <- ifelse(abs(pcors) < 0.05, 0, pcors )
#'
#' est <- enr(net = pcors, n = 500, replications = 2)
#'
#' # plot
#' plot_enr(est)
#'}
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_density scale_x_continuous xlab scale_y_continuous
#'
plot_enr <- function(x, iter = 100000,
                     fill = "#009E73",
                     alpha = 0.5, ...){

  # random variable
  y <- poibin::rpoibin(iter, pp = x$p_s)

  # data frame
  dat <- data.frame(y = y)

  # plot
  ggplot(dat, aes(x = y)) +
    geom_vline(xintercept = mean(y),
               linetype = "dotted") +
    geom_density(adjust = 2,
                 fill =  fill,
                 alpha = alpha) +
    scale_x_continuous(
      limits = c(min(y), max(y)),
      breaks = seq(min(y), max(y),
                   length.out = 5),
      labels = paste0(round(
        seq(min(y), max(y),
            length.out = 5) / x$n_nonzero, 2
      ) * 100, "%")
    ) +
    xlab("Replicated Edges") +
    scale_y_continuous(expand = c(0, 0))

}


