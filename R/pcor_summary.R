#' @title S3 pcor_summary method
#' @name pcor_summary
#' @description S3 pcor_summary method
#' @param object object of class \code{GGM_fisher_z} or \code{GGM_bootstrap}
#' @param ... not currently used
#'
#'
#' @return \code{select} works with the following methods:
#' \itemize{
#' \item \code{\link{GGM_fisher_z}}
#' \item \code{\link{GGM_bootstrap}}
#' }
#' @export
pcor_summary <- function(object,...){
  UseMethod("pcor_summary", object)
}

#' @name pcor_summary.GGM_fisher_z
#' @title Summary Method for \code{pcor_summary.GGM_fisher_z} objects
#' @param object object of class \code{pcor_summary.GGM_fisher_z}
#' @param ... currently ignored
#'
#' @return data frame of partial correlations with confidence intervals
#' @export
#'
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_fisher_z(X)
#'
#'# summarise the partial correlations
#'pcor_summ <- pcor_summary(fit)
#'pcor_summ
pcor_summary.GGM_fisher_z <- function(object,...){
  x <- object
  mat <- matrix(0, x$p, x$p)
  res <- x$fisher_z_results$cis
  mat[] <-unlist( lapply(1:x$p, function(z) paste(1:x$p,z, sep = "--")))
  edge_names <- mat[upper.tri(mat)]
  ci <- 1 - x$alpha
  edge_summary <- data.frame(Edge = edge_names,
                             Estimate = round(x$fisher_z_results$cis$pcor,3),
                             low = round(x$fisher_z_results$cis$low,3),
                             up = round(x$fisher_z_results$cis$up, 3))
  colnames(edge_summary) <- c("Edge", "Estimate",
                              paste(c("lb.", "ub."), gsub("*0.","", x$ci), "%", sep = ""))

  returned_object <- list(edge_summary = edge_summary, object = object)
  class(returned_object) <- "pcor_summary.GGM_fisher_z"
  return(returned_object)

}

#' @name print.pcor_summary.GGM_fisher_z
#' @title Print Method for \code{pcor_summary.GGM_fisher_z} objects
#' @param x object of class \code{pcor_summary.GGM_fisher_z}
#' @param ... currently ignored
#'
#' @export
#'
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_fisher_z(X)
#'
#'# summarise the partial correlations
#'pcor_summ <- pcor_summary(fit)
#'pcor_summ
print.pcor_summary.GGM_fisher_z <- function(x,...){
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Fisher z Transformation \n")
  cat(paste("Alpha level:", x$object$alpha, "\n"))
  cat(paste("FDR:", x$object$FDR, "\n"))
  cat("----\n\n")
  cat("Edge Summary\n")
  print(x$edge_summary, row.names = FALSE)
}

#' Plot \code{pcor_summary.GGM_fisher_z} Edges
#'
#' @param x object of class \code{pcor_summary.GGM_fisher_z }
#' @param size point size
#' @param color point color
#' @param customize \code{"GGMnonreg"} for customized plot and \code{user} for a bare plot to be further
#' customized
#' @param ... currently ignored
#' @import ggplot2
#' @return a \code{ggplot} object
#' @export
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_fisher_z(X)
#'
#'# summarise the partial correlations
#'pcor_summ <- pcor_summary(fit)
#'
#'# plot
#'plot(pcor_summ)

plot.pcor_summary.GGM_fisher_z <- function(x,
                                           size = 1,
                                           color = "black",
                                           customize = "GGMnonreg",
                                           ...){

  # visible binding
  Edge <- NA
  Estimate <- NA
  sig <- NA

  object <- x
  dat <- object$edge_summary

  dat <- dat[order(dat$Estimate),]
  dat$Edge <- factor(dat$Edge, labels = dat$Edge, levels = dat$Edge)
  dat$sig <- as.factor(ifelse(dat[,3] < 0 & dat[,4]> 0, 0, 1))

  plt <- ggplot(dat, aes(x = Edge,
                         y = Estimate)) +
    geom_errorbar(aes(ymin = dat[,3],
                      ymax = dat[,4], color = sig),
                  width = 0) +
    geom_point(size = size, color = color)
  if(customize == "GGMnonreg"){

    plt <- plt + theme_bw() + theme(legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle =90)) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      ylab("Partial Correlation") +
      scale_color_manual(values = c("grey75", "grey25")) +
      scale_x_discrete(expand = c(0.01, 0.01))

  } else if(customize == "user"){

    plt <- plt
  } else{
    stop("customize must be 'GGMnonreg' or 'user'")
  }

  plt
}

#' @name pcor_summary.GGM_bootstrap
#' @title Summary Method for a \code{pcor_summary.GGM_bootstrap} object
#' @param object object of class \code{pcor_summary.GGM_bootstrap}
#' @param ... currently ignored
#'
#' @return data frame of partial correlations with confidence intervals
#' @export
#'
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_bootstrap(X)
#'
#'# summarise the partial correlations
#'pcor_summ <- pcor_summary(fit)
#'pcor_summ

pcor_summary.GGM_bootstrap <- function(object,...){
  x <- object
  mat <- matrix(0, x$p, x$p)
  res <- x$fisher_z_results$cis
  mat[] <-unlist( lapply(1:x$p, function(z) paste(1:x$p,z, sep = "--")))
  edge_names <- mat[upper.tri(mat)]
  ci <- 1 - x$alpha


  edge_summary <- data.frame(Edge = edge_names,
                             Estimate = round(x$dat_res$mean,3),
                             low = round(x$dat_res[,2],3),
                             up = round(x$dat_res[,3], 3))
  colnames(edge_summary) <- c("Edge", "Estimate",
                              paste(c("lb.", "ub."), gsub("*0.","", ci),
                                    "%", sep = ""))

  returned_object <- list(edge_summary = edge_summary, object = object)
  class(returned_object) <- "pcor_summary.GGM_bootstrap"
  return(returned_object)

}




#' @name print.pcor_summary.GGM_bootstrap
#' @title Print Method for \code{pcor_summary.GGM_bootstrap} objects
#' @param x object of class \code{pcor_summary.GGM_bootstrap}
#' @param ... currently ignored
#'
#' @export
#'
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_bootstrap(X)
#'
#'# summarise the partial correlations
#'pcor_summ <- pcor_summary(fit)
#'pcor_summ
print.pcor_summary.GGM_bootstrap <- function(x,...){
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Non-parameteric Bootstrap \n")
  cat(paste("Alpha level:", x$object$alpha, "\n"))
  cat("----\n\n")
  cat("Edge Summary\n")
  print(x$edge_summary, row.names = FALSE)
}



#' Plot \code{pcor_summary.GGM_bootstrap} Edges
#'
#' @param x object of class \code{pcor_summary.GGM_bootstrap}
#' @param size point size
#' @param color point color
#' @param customize \code{"GGMnonreg"} for customized plot and \code{user} for a bare plot to be further
#' customized
#' @param ... currently ignored
#'
#' @return a \code{ggplot} object
#' @export
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_bootstrap(X)
#'
#'# summarise the partial correlations
#'pcor_summ <- pcor_summary(fit)
#'
#'# plot
#'plot(pcor_summ)

plot.pcor_summary.GGM_bootstrap <- function(x,
                                           size = 1,
                                           color = "black",
                                           customize = "GGMnonreg",
                                           ...){

  # visible binding
  Edge <- NA
  Estimate <- NA
  sig <- NA

  object <- x
  dat <- object$edge_summary
  dat <- dat[order(dat[,2]),]
  dat$Edge <- factor(dat$Edge, labels = dat[,1], levels = dat[,1])
  dat$sig <- as.factor(ifelse(dat[,3] < 0 & dat[,4] > 0, 0, 1))

  plt <- ggplot(dat, aes(x = Edge,
                         y = Estimate)) +
    geom_errorbar(aes(ymin = dat[,3],
                      ymax = dat[,4], color = sig),
                  width = 0) +
    geom_point(size = size, color = color)
  if(customize == "GGMnonreg"){

    plt <- plt + theme_bw() + theme(legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle =90)) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      ylab("Partial Correlation") +
      scale_color_manual(values = c("grey75", "grey25")) +
      scale_x_discrete(expand = c(0.01, 0.01))

  } else if(customize == "user"){

    plt <- plt
  } else{
    stop("customize must be 'GGMnonreg' or 'user'")
  }

  plt
}
