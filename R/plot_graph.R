#' Network Plot for \code{graph} Objects
#'
#' @description  Visualize the conditional (in)dependence structure.
#'
#' @param x An object of class \code{graph} obtained from \code{\link[GGMncv]{get_graph}}.
#'
#' @param layout Character string. Which graph layout (defaults is \code{circle}) ?
#'               See \link[sna]{gplot.layout}.
#'
#' @param neg_col Character string. Color for the positive edges
#'                (defaults to a colorblind friendly red).
#'
#' @param pos_col Character string.  Color for the negative edges
#'                (defaults to a colorblind friendly green).
#'
#' @param edge_magnify Numeric. A value that is multiplied by the edge weights. This increases (> 1) or
#'                     decreases (< 1) the line widths (defaults to 1).
#'
#' @param node_size Numeric. The size of the nodes (defaults to \code{10}).
#'
#' @param palette  A character string sepcifying the palette for the \code{groups}.
#'                (default is \code{Set3}). See \href{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}{palette options here}.
#'
#' @param node_names Character string. Names for nodes of length \emph{p}.
#'
#' @param node_groups A character string of length \emph{p} (the number of nodes in the model).
#'               This indicates groups of nodes that should be the same color
#'               (e.g., "clusters" or "communities").
#'
#' @param ... Currently ignored.
#'
#' @return An object of class \code{ggplot}
#'
#' @export
#'
#' @importFrom ggplot2 scale_color_brewer geom_text guides geom_point
#'
#' @importFrom  network network.vertex.names<- set.edge.value set.edge.attribute %e% %v%<- network
#'
#' @importFrom sna gplot.layout.circle
#'
#' @examples
plot.graph <- function(x,
                       layout = "circle",
                       neg_col = "#D55E00",
                       pos_col = "#009E73",
                       edge_magnify = 1,
                       node_size = 10,
                       palette = 2,
                       node_names = NULL,
                       node_groups = NULL,
                       ...){


  x$pcor_adj <- x$P

  p <- ncol(x$P)

  if(is.null(node_names)){
    cn <- 1:p
  } else {
    cn <- node_names
  }

  diag(x$pcor_adj) <- 0

  net <- network::network(x$pcor_adj)

  network::set.edge.value(x = net, attrname = "weights",
                          value = x$pcor_adj)

  network::set.edge.value(x = net, attrname = "abs_weights",
                          value = abs(x$pcor_adj) * edge_magnify)

  network::set.edge.attribute(x = net, attrname = "edge_color",
                              value = ifelse(net %e% "weights" < 0, neg_col,
                                             pos_col))
  e <- abs(as.numeric(x$pcor_adj))

  plt <-GGally::ggnet2(net, edge.alpha = e[e != 0]/max(e),
                       edge.size = "abs_weights",
                       edge.color = "edge_color",
                       node.size = 1,
                       mode = layout)

  if(is.null(node_groups)){
    plt <- plt + geom_point(color = "black",
                            size = node_size + 1) +
      geom_point(size = node_size,
                 color = "white") +
      guides(color = FALSE) +
      geom_text(label = cn)

  } else {

    plt <-  plt + geom_point(aes(color = node_groups, group = node_groups),
                             size = node_size + 1, alpha = 0.5) +
      geom_point(size = node_size, aes(color = node_groups)) +
      geom_text(label = cn)  +
      scale_color_brewer(palette = palette)


  }
  plt
}


#' Get Graph
#'
#' @description Extract the necessary ingredients to visualize the conditional
#'              dependence structure.
#'
#' @param x An object of class \code{ggmncv}
#'
#' @return A list including two matrices (the weighted adjacency and adjacency matrices)
#'
#' @export
#'
#' @examples
get_graph <- function(x){
  returned_object <- list(P = x$wadj, adj = x$adj)
  class(returned_object) <- "graph"
  return(returned_object)
}
