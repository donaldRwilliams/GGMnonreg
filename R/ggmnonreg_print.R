#' Print \code{ggmnonreg} Object
#'
#' @param x An object of class \code{ggmnonreg}
#'
#' @param ... Currently ignored
#'
#' @return No return value.
#' @export
#' @importFrom methods is
print.ggmnonreg <- function(x,...){

  if(is(x, "ggm_search")){
    colnames(x$wadj) <- 1:ncol(x$wadj)
    print(as.data.frame(x$wadj), ...)
  } else if(is(x, "ggm_inference")){
    colnames(x$wadj) <- 1:ncol(x$wadj)
    print(as.data.frame(x$wadj),...)
  } else if(is(x, "enr")){
    print_enr(x,...)
  } else if(is(x, "predictability")){
    print_r2(round(x$r2, 2),...)
  } else if(is(x, "eip")){
    print_eip(x, 2)
  } else if(is(x, "ising_search")){
    colnames(x$wadj) <- 1:ncol(x$wadj)
    print(as.data.frame(x$wadj), ...)
  } else if(is(x, "mixed_search")){
    colnames(x$wadj) <- 1:ncol(x$wadj)
    print(as.data.frame(x$wadj), ...)
  } else if(is(x, "ggm_compare")){
    colnames(x$wadj) <- 1:ncol(x$wadj)
    print(as.data.frame(x$wadj), ...)
  }
}
