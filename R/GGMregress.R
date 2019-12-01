#' Network with Multiple Regression (Ordinary Least Squares)
#'
#' @param X data Matrix (dimensions n by p)
#' @param IC information criterion (\code{AIC} or \code{BIC})
#' @param method subsect selection method (\code{forward}, \code{backward}, or \code{exhaustive})
#' @param rule decision rule (see notes)
#' @param ... currently ignored
#' @export
#'
#' @note The \code{rule} arguement can either be \code{or} or \code{and}. For the former, only the edge only needs to
#' be selected in one regrssion models, whereas, for the later, each edge must be selected in both regression models.
#' Setting \code{rule = "and"} should result in a network with fewer connections.
#' @examples
#'# data
#'X <- GGMnonreg::ptsd
#'
#'# fit model
#'fit <- GGM_regression(X)
#'
#'# selected partials
#'fit$pcor_selected
#'
#'# selected adjacency
#'fit$adj_selected
GGM_regression.default <- function(X, IC = "BIC",
                                   method = "exhaustive",
                                   rule = "and",...){
  # scale data

  if( IC != "AIC" && IC != "BIC" ){
    stop("IC must be AIC or BIC")
  }

  if( method != "forward" && method != "backward"  && method != "exhaustive"){
    stop("method must be foward, backward, or exhaustive")
  }

  X <- scale(X, scale = T)
  n <- nrow(X)
  p <- ncol(X)

  mat1 <- mat2 <- mat_or <- matrix(0, p, p)
  colnames(mat1) <- paste("X", 1:p, sep = "")
  colnames(X) <-  paste("X", 1:p, sep = "")

  fit_models <- lapply(1:p, function(x) bestglm::bestglm(cbind.data.frame(X[,-x], X[,x]),
                                  method = method, IC = IC, intercept = F))

  estimates <- lapply(1:p, function(x) fit_models[[x]]$BestModel$coefficients)
  for(i in 1:p){
    mat1[i,names(estimates[[i]])] <- estimates[[i]]

  }

  mat2[lower.tri(mat2)] <- sqrt(abs(mat1[lower.tri(mat1)]) * abs(t(mat1)[lower.tri(mat1)]))
  mat2[upper.tri(mat2)] <- t(mat2)[upper.tri(mat2)]

  mat_and <- sign(mat1) * mat2


  up <- mat1[upper.tri(mat1)]
  low <- t(mat1)[upper.tri(mat1)]

  sign_or <- sign(up)

  mat_or[upper.tri(mat_or)] <- mapply(or_helper,low, up) * sign_or
  mat_or[lower.tri(mat_or)] <- t(mat_or)[lower.tri(mat_or)]


  adj_and <- ifelse(mat_and == 0, 0, 1)
  adj_or  <- ifelse(mat_or == 0, 0, 1)

if(rule == "or"){
  returned_object <- list(pcor_selected = mat_or,
       adj_selected = adj_or,
       p = p, IC = IC, rule = rule, method = method,
       regression_results = fit_models, dat = X)


}  else if (rule == "and"){
  returned_object <- list(pcor_selected = mat_and,
       adj_selected = adj_and,
       p = p, IC = IC, rule = rule, method = method,
       regression_results = fit_models, dat = X)

}
  class(returned_object) <- "GGM_regression"
  returned_object


}

#' @title S3 GGM_regression method
#' @name GGM_regression
#' @param ... currently not used
#'
#' @description S3 GGM_regression method
#' @seealso \code{\link{GGM_regression.default}}
#' @export
GGM_regression <- function(...) {
  UseMethod("GGM_regression")
}


#' @name summary.GGM_regression
#' @title Summary method for a \code{GGM_regression} object
#' @param object An object of class \code{GGM_regression}
#' @param ... currently ignored
#' @export
#' @examples
#' X <- GGMnonreg::ptsd[, 1:5]
#' fit <- GGM_regression(X)
#' summary(fit)
summary.GGM_regression <- function(object, ...){
  x <- object
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Mulitple Regression \n")
  cat(paste("Information Criterion:", x$IC, "\n"))
  cat(paste("Rule:", x$rule, "\n"))
  cat("----\n\n")
  cat("Selected Network:\n\n")
  colnames(x$pcor_selected) <- 1:x$p
  row.names(x$pcor_selected) <- 1:x$p
  print(x$pcor_selected)
}


#' @name print.GGM_regression
#' @title Print method for a \code{GGM_regression} object
#' @param x An object of class \code{GGM_regression}
#' @param ... currently ignored
#' @export
#' @examples
#' X <- GGMnonreg::ptsd[, 1:5]
#' fit <- GGM_regression(X)
#' fit
print.GGM_regression <- function(x,...){
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Mulitple Regression \n")
  cat(paste("Information Criterion:", x$IC, "\n"))
  cat(paste("Rule:", x$rule, "\n"))
  cat("----\n")
  cat(date())
}


#' Plot \code{GGM_regression} Network
#'
#' @param x object of class \code{GGM_regression}
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
#'fit <- GGM_regression(X)
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
plot.GGM_regression <-
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
