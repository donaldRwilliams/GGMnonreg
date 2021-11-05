#' @title Gaussian graphical model: statistical inference
#'
#' @description Learn the conditional dependence structure with null hypothesis
#'              significance testing. This provides a valid measure of parameter
#'              uncertainty.
#'
#' @param Y The data matrix of dimensions \emph{n} (observations) by \emph{p} (nodes).
#'
#' @param alpha The desired significance level (defaults to \code{0.05}). Note that
#'              1 - alpha corresponds to specificity.
#'
#' @param control_precision Logical. Should precision (i.e., 1 - false discovery rate)
#'                          be controlled at the level alpha (defaults to \code{FALSE}) ?
#'
#' @param boot Logical. Should a non-parametric bootstrap be employed (defaults to \code{TRUE})?
#'
#' @param B Integer. Number of bootstrap replicates (defaults to \code{1000}).
#'
#' @param cores Integer. Number of cores to be used when executing in parallel
#'              (defaults to 1).
#'
#' @param method Character string. Which type of correlation coefficients
#'               to be computed. Options include \code{"pearson"} (default),
#'               \code{"kendall"}, \code{"spearman"}, and \code{"polychoric"}.
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE})?
#'
#' @return An object of class \code{ggm_inference}
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' Y <- ptsd
#'
#' fit <- ggm_inference(Y)
#'
#' }
#'
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach `%dopar%`
#' @importFrom corpcor cor2pcor
#' @importFrom psych polychoric
#' @importFrom stats quantile cor p.adjust
#' @importFrom MASS mvrnorm
#'
ggm_inference <- function(Y,
                          alpha = 0.05,
                          control_precision = FALSE,
                          boot = TRUE,
                          B = 1000,
                          cores = 1,
                          method = "pearson" ,
                          progress = TRUE){

  ci_lower <- alpha / 2

  ci_upper <- 1 - ci_lower

  p <- ncol(Y)

  n <- nrow(Y)

  adj <- matrix(0, p, p)

  wadj <- matrix(0, p, p)

  pcors <- matrix(0, p, p)

  if(boot){
    if(control_precision){
      warning("can't control precision when 'boot = TRUE'")
    }
    cl <- parallel::makeCluster(cores)

    doParallel::registerDoParallel(cl)

    if(progress){
      pb <- utils::txtProgressBar(max = B, style = 3)
    }

    boot_samps <- foreach::foreach(i = 1:B, .combine = f(B)) %dopar%{

                                     if(method == "polychoric"){

                                       pcors <-  corpcor::cor2pcor(
                                         psych::polychoric(Y[sample(1:n, size = n, replace = TRUE),])$rho
                                       )

                                     } else {

                                       pcors <- corpcor::cor2pcor(
                                         stats::cor(
                                           Y[sample(1:n, size = n, replace = TRUE),],
                                           method = method
                                         )
                                       )
                                     }

                                     pcors <- pcors[upper.tri(diag(p))]

                                     return(pcors)

                                   }

    parallel::stopCluster(cl)

    cis <- apply(boot_samps, 2, function(x){
      stats::quantile(x, probs = c(ci_lower, ci_upper))
    })

    cis <- t(cis)

    boot_mean <- colMeans(boot_samps)

    adj[upper.tri(diag(p))] <- ifelse(cis[,1] < 0 & cis[,2] > 0, 0, 1)

    adj <- symm_mat(adj)

    pcors[upper.tri(diag(p))] <- boot_mean
    pcors <- symm_mat(pcors)

    wadj <- pcors * adj

    returned_object <- list(wadj = wadj,
                            adj = adj,
                            pcors = pcors,
                            boot_samps = boot_samps,
                            boot_mean = boot_mean,
                            cis = cis,
                            n = n, p = p,
                            Y = Y)
  } else {

    if(!method %in% c("pearson", "spearman")){
      stop("analytic solution only available for 'pearson' and 'spearman'")
    }

    pcors <- corpcor::cor2pcor(cor(Y, method = method))

    z <- abs(GGMnonreg::fisher_r_to_z(pcors)[upper.tri(diag(p))])

    if(method == "spearman"){

      test_stat <- z / sqrt(1.06 / (n - p - 2 - 3))

    } else {

      test_stat <- z / sqrt(1 / (n - p - 2 - 3))

    }

    p_values <- pnorm(test_stat, lower.tail = F) * 2

    if(control_precision){

      p_values <-  p.adjust(p_values, method = "fdr")

    }

    adj[upper.tri(diag(p))] <- ifelse(p_values < alpha, 1, 0)

    adj <- symm_mat(adj)

    wadj <- adj * pcors

    returned_object <- list(wadj = wadj,
                            adj = adj,
                            pcors = pcors,
                            n = n, p = p,
                            Y = Y)

  }

  class(returned_object) <- c("ggmnonreg",
                              "ggm_inference")
  return(returned_object)

}

