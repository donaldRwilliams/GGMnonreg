#' Expected Network Replicability
#'
#' @param net true network
#' @param n samples size
#' @param replications number of networks
#' @param type pearson or spearman
#'
#' @return
#' @export
#'
#' @examples
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
#' enr(net = pcors, n = 500, replications = 2)
enr <- function(net, n, replications, type = "pearson"){

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

  qs <- round( quantile(1:length(which_nonzero),
                        probs = seq(0, 0.9, 0.1)))

  # probability of success
  p_s <- power_z(pcs[which_nonzero],
                 n = n,
                 c = c,
                 type = type)^replications

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

  class(returned_object) <- "enr"
  returned_object
}

#' Print Method for \code{enr} Objects
#'
#' @param x object of class \code{enr}
#' @param ... currently ignored
#'
#' @return
#' @export
print.enr <- function(x,...){
  cat("GGMnonreg: Nonregularized GGMs \n")
  cat("Method: Expected Network Replicability \n")
  cat("Nodes:", x$p, "\n")
  cat("Networks:", x$replications, "\n")
  cat("Sample Size:", x$n, "\n")
  cat("Number of Edges:",  x$n_nonzero, "\n")
  cat("----\n\n")
  cat("Average Replicability:", round(x$ave_pwr, 2), "\n")
  cat("Average Number of Edges:",
      round(round(x$ave_pwr, 2) * x$n_nonzero),
      paste0( "(SD = ", round(sqrt(x$var_pwr), 2), ")"), "\n\n")
  cat("----\n\n")
  cat("Cumulative Probability:" , "\n\n")

  dat <- data.frame(prop = seq(0, .90, .10),
                    round(quantile(1:x$n_nonzero,
                                   probs = seq(0, 0.9, 0.1))),
                    prob = round(x$cdf, 2))

  colnames(dat) <- c("Proportion", "Edges", "Probability")
  print(dat, row.names = F, right = T)

  cat("----\n")
  cat(paste0( "note: \nProbability that more edges than corresponding proportion \n",
              "and number of edges are replicated \n"))
  }


#' Summary Method for \code{enr} Objects
#'
#' @param object object of class \code{enr}
#' @param ... currently ignored
#'
#' @return
#' @export
summary.enr <- function(object,...){
  print(object)
}


#' Plot \code{enr} Objects
#'
#' @param x Object of class \code{enr}
#' @param iter number of samples
#' @param fill fill color for density
#' @param alpha transparency
#' @param ... currently ignored
#'
#' @return
#' @export
#' @examples
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
#' plot(est)
plot.enr <- function(x, iter = 100000,
                     fill = "#009E73",
                     alpha = 0.5, ...){

  # random variable
  y <- poibin::rpoibin(iter, pp = x$p_s)

  # data frame
  dat <- data.frame(y = y)

  # plot
  ggplot(dat, aes(x = y )) +
    geom_vline(xintercept = mean(y),
               linetype = "dotted") +
    geom_density(adjust = 2,
                 fill =  fill, alpha = alpha) +
    scale_x_continuous(limits = c(min(y), max(y) ),
                       breaks = seq(min(y), max(y), length.out = 5),
                       labels = paste0(round(seq(min(y), max(y),
                                                  length.out = 5) / x$n_nonzero, 2)*100, "%")) +
  xlab("Replicated Edges") +
  scale_y_continuous(expand = c(0,0))

}


