#' Title
#'
#' @param Theta
#' @param hatTheta
#'
#' @return
#' @export
#'
#' @examples
KL = function(Theta,hatTheta){
  p = ncol(Theta)

  invTheta = solve(Theta,diag(1,p))

  kl  = sum(diag(invTheta%*%hatTheta)) - log(det(invTheta%*%hatTheta)) - p

  return(kl)

}
