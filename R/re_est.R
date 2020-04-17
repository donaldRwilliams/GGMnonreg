#' Estimate Precision Matrix with a Known Structure
#'
#' @param S Covariance matrix
#' @param adj Adjancency matrix
#' @param tol tolerance
#'
#' @return
#' @export
#'
#' @examples
re_est <- function(S, adj, tol = 1e-10 ) {

  # nodes
  p <- ncol(S)


  # set diag to zero
  diag(adj) <- 0

  # set W
  W <- W_previous <- S
  n_iter <- 0
  repeat{

    for(i in 1:p){

      # beta storage
      beta <- rep(0, p - 1)

      # pad index
      pad_index <- which(adj[i,-i] == 1)

      if(length(pad_index) == 0){

        w_12 <- beta

      } else {

        # submatrix (p-1) by (p-1)
        W_11 <- W[-i,-i]

        # (p-1) by 1
        s_12 <- S[i,-i]

        # padded submatrix
        W_11_star <- W_11[pad_index, pad_index]

        # padded covariance
        s_12_star <- s_12[pad_index]

        # s_12_star = X'y / (n - 1)
        # W_11_star^-1 = (X'X)^-1 * (n - 1)
        beta[pad_index] <- solve(W_11_star) %*% s_12_star

        # updated covariance
        w_12 <- W_11 %*% beta
      }

      W[-i,i] <- W[i,-i] <- w_12

    }

    max_diff <-  max(W_previous[upper.tri(W)] - W[upper.tri(W)])

    # check to break
    if(max_diff < tol){
      break
    } else {

      W_previous <- W
      n_iter <- n_iter + 1
    }
  }

  returned_object <- list(Theta = round(solve(W), 4),
                          Sigma = round(W, 4))
  returned_object
}


