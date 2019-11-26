#' Network with a Non-parametric Bootstrap
#' @param X data matrix (dimensions n by p)
#' @param iter number of bootstrap samples
#' @param alpha desired type I error rate (correspondig to approximately 1 - specificity)
#'
#' @return list of class GGM_bootstrap
#' \itemize{
#' \item \code{pcor_selected} selected partial correlation matrix
#' \item \code{adj_selected} adjacency matrix
#' \item \code{pcor_mean} saturated matrix (bootstrap mean)
#' \item \code{boot_results} matrix of bootstrap samples
#' \item \code{dat} data
#' }
#' @export
#'
#'
#' @examples
#'
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_bootstrap(X)
#'
#'# selected partials
#'fit$pcor_selected
#'
#'# selected adjacency
#'fit$adj_selected
#'
#'# plot results
#'plot(fit)
GGM_bootstrap.default <- function(X, iter = 1000, alpha = 0.05){


  X <- as.matrix(scale(na.omit(X)))

  p <- ncol(X)

  mat_selected <- mat_mean <- matrix(0, p, p)

  lw_bound <- alpha / 2

  up_bound <- 1 -   lw_bound

  boot_results <- lapply(1:iter, function(x) boot_inv_helper(X))

  boot_pcor <- do.call(rbind, lapply(1:iter, function(x) boot_results[[x]]$upper_tri))

  quantiles <- t(apply(boot_pcor, MARGIN = 2,
                       quantile, probs = c(lw_bound, up_bound)))

  means <-  t(apply(boot_pcor, MARGIN = 2,  mean))

  mat_selected[upper.tri(mat_selected)] <- mapply(quantile_helper,
                                                  x = quantiles[,1],
                                                  y = quantiles[,2])

  mat_selected[lower.tri(mat_selected)] <- t(mat_selected)[lower.tri(mat_selected)]

  mat_mean[upper.tri(mat_mean)] <-  means
  mat_mean[lower.tri(mat_mean)] <- t(mat_mean)[lower.tri(mat_mean)]

  returned_object <- list(pcor_selected = mat_selected * mat_mean,
                          adj_selected = mat_selected,
                          pcor_mean =  mat_mean,
                          boot_results = boot_results,
                          dat = X, iter = iter, p = p)

  class(returned_object) <- "GGM_bootstrap"
  return(returned_object)
}



#' @title S3 GGM_bootstrap method
#' @name GGM_bootstrap
#' @param ... currently not used
#'
#' @description S3 GGM_bootstrap method
#' @seealso \code{\link{GGM_bootstrap.default}}
#' @export
GGM_bootstrap <- function(...) {
  UseMethod("GGM_bootstrap")
}



#' Plot \code{GGM_bootstrap} Network
#'
#' @param x object of class \code{GGM_bootstrap}
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
#'fit <- GGM_bootstrap(X)
#'
#'# selected partials
#'fit$pcor_selected
#'
#'# selected adjacency
#'fit$adj_selected
#' plt <- plot(E,
#'             node_labels = letters[1:20],
#'             node_labels_color = "white",
#'             node_groups = rep(c("1", "2", "3", "4"), each = 5),
#'             edge_colors = "classic",
#'             edge.alpha = 0.5, palette = "Set2")
plot.GGM_bootstrap <-
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
