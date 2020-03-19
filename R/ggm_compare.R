#' Compare Gaussian Graphical Models
#'
#' @param Y_g1 data matrix or data frame of dimensions n by p (group 1)
#' @param Y_g2 data matrix or data frame of dimensions n by p (group 1)
#' @param alpha desired type I error rate (corresponding to approximately 1 - specificity)
#' @param control_precision should precision be controlled ?
#' This corresponds to controlling 1 - the false discovery rate (FDR). Default is \code{NULL}.
#' @param precision desired precision (1 - FDR). Default is 0.90 (1 - 0.9 = FDR of 0.1).
#'
#' @return list containing the adjacency matrix (0's and 1's), partial correlation differences, and an
#' 'adjacency' matrix containing the partial correlation differences.
#' @export
#'
#' @examples
#' true <- gen_pcors(20, edge_prob = 0.3)
#' Y1 <- MASS::mvrnorm(100,
#'                     mu = rep(0, 20),
#'                     Sigma = true$cors)
#'
#' Y2 <- MASS::mvrnorm(100,
#'                    mu = rep(0, 20),
#'                    Sigma = diag(20))
#'
#' fit <- ggm_compare(Y1, Y2, precision = 0.80,
#'                    control_precision = TRUE)


ggm_compare <- function(Y_g1, Y_g2, alpha = 0.05,
                        control_precision = NULL,
                        precision = 0.90){


  n1 <- nrow(Y_g1)
  n2 <- nrow(Y_g2)
  c <- ncol(Y_g1) - 2
  p <- ncol(Y_g1)

  fit1 <- GGM_fisher_z(Y_g1)
  fit2 <- GGM_fisher_z(Y_g2)

  pcors1 <- fit1$fisher_z_results$par_cors
  pcors2 <- fit2$fisher_z_results$par_cors

  var1 <- 1 / (n1 - c - 2)
  var2 <- 1 /  (n2 - c - 2)

  var_diff <- var1 + var2
  se_diff <- sqrt(var_diff)

  fz_diff <- fisher_r2z( pcors1[upper.tri(pcors1)]) -  fisher_r2z(pcors2[upper.tri(pcors2)])

  test_stat <- fz_diff / se_diff

  pvalues <- pnorm(abs(test_stat), lower.tail = FALSE) * 2

  adj <- matrix(0, p, p)

  if(is.null(control_precision)){

    adj[upper.tri(adj)] <- ifelse(pvalues < alpha, 1, 0)
    adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]

    pcors_diff <- pcors1 - pcors2
    pcors_adj <- adj * pcors_diff

    returned_object <- list(adj = adj,
                            pcors_adj = pcors_adj,
                            pcors_diff = pcors_diff)
    returned_object

  } else {


    adj[upper.tri(adj)] <- ifelse(stats::p.adjust(pvalues,
                                                  method = "fdr", ) < (1 - precision), 1, 0)
    adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]

    pcors_diff <- pcors1 - pcors2
    pcors_adj <- adj * pcors_diff

    returned_object <- list(adj = adj,
                            pcors_adj = pcors_adj,
                            pcors_diff = pcors_diff)
    returned_object

  }

  returned_object

}
