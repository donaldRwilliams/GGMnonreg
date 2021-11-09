#' Extract Confidence Intervals from \code{ggm_inference} Objects
#'
#' @param object An object of class \code{ggm_inference}.
#'
#' @param ... Currently ignored.
#'
#' @return A matrix including bootstrap confidence intervals.
#' @export
#'
#' @examples
#' \donttest{
#' #  data
#' Y <- ptsd
#'
#' # eip
#' fit <- ggm_inference(Y, method = "spearman",
#' boot = TRUE, B = 100)
#'
#' # cis
#' confint(fit)
#'
#' }
confint.ggm_inference <- function(object, ...){

  if(is.null(object$boot_samps)) {
    stop("currently only implemented for 'boot = TRUE'")
  }

  return(object$cis)
}
