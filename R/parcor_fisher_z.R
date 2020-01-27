#' Network with Fisher Z Transformation
#'
#' @param X data matrix (dimensions n by p)
#' @param alpha desired type I error rate (correspondig to approximately 1 - specificity)
#' @param FDR false discovery rate control (see notes)
#' @param ... currently ignored
#'
#' @importFrom Matrix forceSymmetric
#' @importFrom stats p.adjust
#'
#' @return list of class GGM_fisher_z
#' \itemize{
#' \item \code{pcor_selected} selected partial correlation matrix
#' \item \code{adj_selected} adjacency matrix
#' \item \code{pcor_mean} saturated matrix (bootstrap mean)
#' \item \code{boot_results} matrix of bootstrap samples
#' \item \code{dat} data
#' }
#'
#' @note When \code{FDR = FALSE}, this results in a customary significance test with the desired alpha level. This can be used to control
#' specificity, that is 1- alpha, which is the false positive rate. On the other hand, when \code{FDR = FALSE},
#' the false discovery rate is controlled at the desired alpha level. This is
#' the proporiton of "discoveries" (included edges) that are actually false positives.

#'
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_fisher_z(X)
#'
#'# selected partials
#'fit$pcor_selected
#'
#'# selected adjacency
#'fit$adj_selected
#'
#'# plot results
#'plot(fit)
#' @export
GGM_fisher_z <- function(X, alpha = 0.05, FDR = FALSE,...){


  n <- nrow(X)
  p <- ncol(X)
  mle_hat <- mle(X)

  ## compute maximum likelihood estimator for covariance matrix
  mle_cov <- mle_hat$cov_mle

  ## compute maximum likelihood estimator of precision matrix
  ## (inverse covariance matrix)
  mle_inv <- mle_hat$inv_cov_mle

  ## standardize and revese sign = partial correaltions
  pc  <-  as.matrix(cov2cor(mle_inv))  * - 1

  mle_pcors <- ci_par_cor(alpha = alpha, par_cors = pc, n = n, s = p - 2)

  mle_inv <- mle_pcors$sig_mat * mle_inv



  if(isTRUE(FDR)){
    mat_fdr <- matrix(0, p, p)
    pvalues <- mle_pcors$pmat[upper.tri(mle_pcors$pmat)]
    mat_fdr[upper.tri(mat_fdr)] <- ifelse(stats::p.adjust(pvalues, method = "fdr") < alpha, 1, 0)
    mat_fdr <- as.matrix(Matrix::forceSymmetric(mat_fdr))
    pcor_selected <- pc * mat_fdr
    diag(pcor_selected) <- 0

    returned_object <- list(pcor_selected = pcor_selected,
                            adj_selected = mat_fdr, p = p,
                            alpha = alpha, dat = X, fisher_z_results = mle_pcors,
                            FDR = FDR)
  } else{

    diag(mle_pcors$par_sig) <- 0
    diag(mle_pcors$sig_mat) <- 0

    returned_object <- list(pcor_selected = mle_pcors$par_sig,
                            adj_selected = mle_pcors$sig_mat,
                            p = p, alpha = alpha, dat = X,
                            fisher_z_results = mle_pcors, FDR = FDR)

  }

  class(returned_object) <- "GGM_fisher_z"
  returned_object


}
#' @name summary.GGM_fisher_z
#' @title Summary method for a \code{GGM_fisher_z} object
#' @param object An object of class \code{GGM_fisher_z}
#' @param ... currently ignored
#' @export
#' @examples
#' X <- GGMnonreg::ptsd[, 1:5]
#' fit <- GGM_fisher_z(X)
#' summary(fit)
summary.GGM_fisher_z <- function(object, ...){
  x <- object
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Fisher z Transformation \n")
  cat(paste("Alpha level:", x$alpha, "\n"))
  cat(paste("FDR:", x$FDR, "\n"))
  cat("----\n\n")
  cat("Selected Network:\n\n")
  colnames(x$pcor_selected) <- 1:x$p
  row.names(x$pcor_selected) <- 1:x$p
  print(x$pcor_selected)
}

#' @name print.GGM_fisher_z
#' @title Print method for a \code{GGM_fisher_z} object
#' @param x An object of class \code{GGM_fisher_z}
#' @param ... currently ignored
#' @export
#' @examples
#' X <- GGMnonreg::ptsd[, 1:5]
#' fit <- GGM_fisher_z(X)
#' fit
print.GGM_fisher_z <- function(x,...){
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Fisher z Transformation \n")
  cat(paste("Alpha level:", x$alpha, "\n"))
  cat(paste("FDR:", x$FDR, "\n"))
  cat("----\n")
  cat(date())
}


#' Plot \code{GGM_fisher_z} Network
#'
#' @param x object of class \code{GGM_fisher_z}
#' @param layout network layout (\link[sna]{gplot.layout})
#' @param edge_colors color theme for positive and negative edges ("classic", "color_blind", or "vivid")
#' @param node_labels node labels
#' @param node_labels_color node labels color
#' @param node_groups node group indicator
#' @param node_outer_size node border size
#' @param node_inner_size node size
#' @param ... additional arguments (\link[GGally]{ggnet2})
#'
#' @importFrom GGally ggnet2
#' @importFrom ggplot2 ggtitle
#' @importFrom network network.vertex.names<- set.edge.value set.edge.attribute %e% %v%<-
#' @importFrom sna gplot.layout.circle
#' @return object of class \code{ggplot}
#'
#' @note See palette options here \url{http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually#use-rcolorbrewer-palettes}.
#' @export
#'
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_fisher_z(X)
#'
#'# selected partials
#'fit$pcor_selected
#'
#'# selected adjacency
#'fit$adj_selected
#' plt <- plot(fit,
#'             node_labels = letters[1:20],
#'             node_labels_color = "white",
#'             node_groups = rep(c("1", "2", "3", "4"), each = 5),
#'             edge_colors = "classic",
#'             edge.alpha = 0.5, palette = "Set2")
plot.GGM_fisher_z <-
  function(x,
           layout = "circle",
           edge_colors = "classic",
           node_labels = NULL,
           node_labels_color = "black",
           node_groups = NULL,
           node_outer_size = 12,
           node_inner_size = 11,
           ...) {
    label <- NULL
    color <- NULL
    # number of nodes
    p <- ncol(x$adj_selected)

    # network
    net <- network::network(x$adj_selected)


    # default labels
    if (is.null(node_labels)) {
      network::network.vertex.names(net) <- 1:p
      # custom labels
    } else {
      # check label length
      if (length(node_labels) != p) {
        stop("labels must be of length p (number of nodes)")
      }

      network::network.vertex.names(net) <- node_labels
    }
    # set edge weights
    network::set.edge.value(
      x = net,
      attrname =  "weights",
      value = x$pcor_selected
    )
    #
    network::set.edge.value(
      x = net,
      attrname =  "abs_weights",
      value = abs(x$pcor_selected) * 10
    )
    #
    #
    if (edge_colors == "classic") {
      network::set.edge.attribute(
        x = net,
        attrname = "edge_color",
        value = ifelse(net %e% "weights" < 0,
                       "brown3", "palegreen3")
      )
    } else if (edge_colors == "color_blind") {
      network::set.edge.attribute(
        x = net,
        attrname = "edge_color",
        value = ifelse(net %e% "weights" < 0,
                       "#009E73", "#D55E00")
      )
    } else if (edge_colors == "vivid") {
      network::set.edge.attribute(
        x = net,
        attrname = "edge_color",
        value = ifelse(net %e% "weights" < 0,
                       "darkorange1", "darkorchid4")
      )

    }

    if (is.null(node_groups)) {
      plt <- ggnet2(
        net = net,
        mode = layout,
        node.color = "white",
        edge.color = "edge_color",
        edge.size = "abs_weights",
        label = TRUE,
        ...
      ) +
        geom_point(color = "black",
                   size = node_outer_size,
                   alpha = 1) +
        geom_point(color = "white",
                   size = node_inner_size,
                   alpha = 1) +
        geom_text(aes(label = label),
                  color = node_labels_color)

      plt <- list(plt = plt)

    } else {
      if (length(node_groups) != p) {
        stop("labels must be of length p (number of nodes)")
      }

      net %v% "group" <- node_groups

      plt <- ggnet2(
        net = net,
        mode = layout,
        node.color = "group",
        edge.color = "edge_color",
        edge.size = "abs_weights",
        label = TRUE,
        ...
      ) +
        geom_point(aes(color = color),
                   size = node_outer_size,
                   alpha = 0.5) +
        geom_point(aes(color = color),
                   size = node_inner_size,
                   alpha = 1) +
        geom_text(aes(label = label),
                  color = node_labels_color)

      plt <- list(plt = plt)
    }

    return(plt)
  }
