#' Gaussian Graphical Model Meta-Analysis
#'
#' @name GGMmeta
#'
#' @description Function to fit meta-analytic fixed- or random-effects models for Gaussian graphical
#' models (i.e., a meta-analysis specifically for partial correlations)
#'
#' @param R A list containing correlation matrices of dimensions \emph{p} by \emph{p}
#'          (at least two are required). Note that this is required if \code{data} or
#'          \code{fz} is not provided.
#'
#' @param data A list containing data matrices of dimensions \emph{n} by \emph{p}
#'             (at least two are required). Note that this is required if \code{R}
#'             or \code{fz} is not provided.
#'
#' @param fz A matrix of Fisher-z transformed partial correlations. Each column should be
#'           the upper-triangular elements of a partial correlation matrix for a given
#'           study. This is required if \code{R} and \code{data} are not supplied.
#'
#' @param fz_var A matrix of the variances for the Fisher-z transformed partial correlation.
#'
#' @param n Numeric vector. The sample sizes if \code{R} is supplied.
#'
#' @param metafor_opts A list of options passed to \code{\link[metafor]{rma}}. The current default
#'                     is to fit a fixed-effects meta-analysis. See \strong{Details}.
#'
#' @param method Character string. Which correlation coefficient should be computed.
#'               One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#'
#' @param alpha Numeric. The desired type I error rate.
#'
#' @param blup Logical. Should the study specifc random effects be computed (defaults to \code{FALSE})?. This
#'             provide a graph based on the Best Linear Unbiased Predictions \code{\link[metafor]{blup}}.
#'
#' @param p_adjust Character string. Correction method for the p-values \code{\link[stats]{p.adjust}}.
#'                 The default is \code{"fdr"} (false discovery rate, i.e., 1 minus precision).
#'
#' @param store Logical. Should the meta-analyses be saved (defaults to \code{FALSE})?
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE})?
#'
#' @importFrom stats p.adjust cor
#' @importFrom metafor rma blup
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#'
#' @return An object of class \code{GGMmeta}, including
#'
#' \itemize{
#' \item \code{one}
#' }
#'
#' @details
#' It is commonly thought that random-effects models are needed in the presence of
#' between-study variance. This is incorrect, as a fixed-effects model can be used
#' and has nominal error rate. The key difference is the desired inference: a
#' fixed-effects model is concerned with only the included studies whereas a
#' random-effects model generalizes to the "population" of studies. Because this function
#' is likely to be used with few studies, generalizing is perhaps a bit ambitous. Hence,
#' the current default is to fit a fixed-effects model which provides inference about about
#' the average effect for those studies included in the analysis.
#'
#'
#' @examples
#'
#' # change to REML to get study specific graphs
#' meta_fit <- GGMmeta(R = list(ptsd_cor1,
#'                           ptsd_cor2,
#'                           ptsd_cor3,
#'                           ptsd_cor4),
#'                           n = c(526, 365, 926, 956),
#'                           blup = TRUE,
#'                           metafor_opts = list(method = "REML"))
#'
#' # overall effect (fdr correction)
#' qgraph::qgraph(meta_fit$wadj_correct)
#'
#' # shrunken networks (study one)
#' qgraph::qgraph(meta_fit$wadj_study[[1]])
#' @export
GGMmeta <- function(R = list(),
                    data = list(),
                    fz = NULL,
                    fz_var = NULL,
                    n = NULL,
                    metafor_opts = list(method = "FE"),
                    method  = "pearson",
                    alpha = 0.05,
                    blup = FALSE,
                    p_adjust = "fdr",
                    store = FALSE,
                    progress = TRUE){

  if(is.null(fz) & is.null(fz_var)){

    # assing lists
    data <- data
    R <- R

  # check for correlation or data matrices
  if (length(R) != 0) {
    groups <- length(R)
    p <- ncol(R[[1]])
    k <- p - 2;
    type <- "cor"

    if(is.null(n)){
      stop("missing a vector including the sample size for each group, e.g., `n = c(100, 500, 1000)` )")
    }

  } else {

    groups <- length(data)
    p <- ncol(data[[1]])
    k <- p - 2;
    type <- "data"
    n <- sapply(1:groups, function(x) nrow(na.omit(data[[x]])))
  }

  # off-diagonals
  relations <- p * (p-1) * 0.5

  # effects and variances
  meta_info <- sapply(1:groups, function(g){

    if(type == "data"){
      Yg <- data[[g]];

      pcors <- cor_2_pcor(stats::cor(x = Yg,
                              method = method));
    } else {

      Rg <- R[[g]];

      pcors <- cor_2_pcor(Rg);

    }

    # note: kendall is approximately normal
    if(method == "kendall"){

      vars <-  (2 * (2 * (n[g] - k) + 5)) /(9 * (n[g] -  k) * (n[g] - 1 - k));

    } else {

      pcors <- fisher_r2z(pcors[upper.tri(pcors)]);
      vars <-  1 / (n[g] - k - 3 );

    }

    # return fisher z and variances
    data.frame(y = pcors, v = vars);
    },
    simplify = FALSE)

  # info into matrices
  fz <- sapply(meta_info, "[[", 1)
  fz_var <- sapply(meta_info, "[[", 2)

  }

  if(progress){
    message("Fitting Meta-Analyses")
    pb <- utils::txtProgressBar(min = 0, max = relations, style = 3)
  }
  metas <- lapply(1:relations, function(x) {

    fit <- do.call(metafor::rma,
                   do.call(c, list(list(yi = fz[x,],
                                        vi = fz_var[x,]),
                                   metafor_opts)))
    if(progress){
    utils::setTxtProgressBar(pb, x)
    }
    fit
    })

  # results into data frame
  meta_results <-
    t(sapply(1:relations, function(x) {
      fit <- metas[[x]]
      data.frame(
        effect =  fit$b,
        ci_lb = fit$ci.lb,
        ci_ub = fit$ci.ub,
        p_value = fit$pval,
        tau = sqrt(fit$tau2),
        I2 = fit$I2,
        Qtest = fit$QEp
      )
    }, simplify = FALSE))

  # combine
  meta_results <- do.call(rbind.data.frame, meta_results)

  # adjusted p-value
  meta_results$adjusted_p_value <-
    stats::p.adjust(meta_results$p_value,
                    method = p_adjust)

  # matrices for storage
  mat_pcor <-
    adj_nocorrect <-
    adj_correct <-
    matrix(0, nrow = p, ncol = p)

  mat_pcor[upper.tri(mat_pcor)] <- meta_results$effect
  mat_pcor <- symm_mat(mat_pcor)

  if (method != "kendall") {
    mat_pcor <- fisher_z2r(mat_pcor)
  }

  adj_nocorrect[upper.tri(mat_pcor)] <-
    ifelse(meta_results$p_value < alpha, 1, 0)
  adj_nocorrect <- symm_mat(adj_nocorrect)

  adj_correct[upper.tri(adj_correct)] <-
    ifelse(meta_results$adjusted_p_value < alpha, 1, 0)
  adj_correct <- symm_mat(adj_correct)


  wadj_nocorrect <- mat_pcor * adj_nocorrect
  wadj_correct <- mat_pcor * adj_correct

  pcor_study <- list()
  adj_study <-  list()
  wadj_study <- list()


  if(blup){

    blups <- sapply(1:relations, function(x) metafor::blup(metas[[x]]), simplify = FALSE)

    for (i in 1:groups) {
      mat_mu <- mat_z <- matrix(0, p, p)
      mat_mu[upper.tri(mat_mu)] <-  sapply(blups, "[[", 1)[i, ]
      mat_z[upper.tri(mat_z)] <- sapply(blups, "[[", 1)[i, ] / sapply(blups, "[[", 2)[i, ]
      mat_mu <- symm_mat(mat_mu)
      mat_mu <- fisher_z2r(mat_mu)
      mat_z <- symm_mat(mat_z)
      adj <- ifelse(abs(mat_z) > 1.959964, 1, 0)
      wadj <- mat_mu * adj
      pcor_study[[i]] <- mat_mu
      adj_study[[i]] <- adj
      wadj_study[[i]] <- wadj
    }


  }

  if(!store){
    metas <- NULL
    blups <- NULL
  }

  returned_object <- list(
    pcors = mat_pcor,
    wadj_correct = wadj_correct,
    wadj_nocorrect = wadj_nocorrect,
    adj_correct = adj_correct,
    adj_nocorrect = adj_nocorrect,
    meta_results = meta_results,
    metas = metas,
    pcor_study = pcor_study,
    adj_study = adj_study,
    wadj_study = wadj_study,
    blup = blups)

  returned_object
}





