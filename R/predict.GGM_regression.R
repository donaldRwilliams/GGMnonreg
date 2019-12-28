#' Network Predictability
#'
#' @param object  an object of class \code{GGM_regression}
#' @param new_data optional new data for out-of-sample predictability
#' @param resample bootstrap predictability (see notes)
#' @param iter number of bootstrap iteractions
#' @param ... currently ignored
#'
#' @note Currently it is only possible to bootstrap predictability, that is variance explained and mean squared
#' error. This provides bootstrapped distribution for each. Computing a point estimate will be implemented
#' in the future
#' @return list of class GGM_bootstrap
#' \itemize{
#' \item \code{R2_boot_storage} matrix of bootstrapped estimates
#' \item \code{mse_boot_storage} matrix of bootstrapped estimates
#' \item \code{iter} number of bootstrap samples
#' \item \code{p} number of nodes
#' \item \code{n} number of observations
#' \item \code{newdata} \code{TRUE} or \code{FALSE}
#' \item \code{reg_object} object of class \code{GGM_regression}
#' \item \code{resample} \code{TRUE} or \code{FALSE}
#' }

#' @export
#'
#' @examples
#'# data
#'X <- GGMnonreg::ptsd[,1:5]
#'
#'# fit model
#'fit <- GGM_regression(X)
#'
#'# predict
#'net_pred <- predict(fit, iter = 10)
#'net_pred
predict.GGM_regression <- function(object,
                                   new_data = NULL,
                                   resample = TRUE,
                                   iter = 500,
                                   ...){

  dat <- object$dat
  p <- ncol(dat)
  n <- nrow(dat)

  # bootstrap
  if(isTRUE(resample)){
    newdata <- FALSE
    pb <- txtProgressBar(min = 0, max = iter, style = 3)

    R2_boot_storage <- mse_boot_storage <- matrix(0, nrow = iter, ncol = p)

    if(is.null(new_data)){

      for(i in 1:iter){

        setTxtProgressBar(pb, i)

        boot_dat <- dat[sample(1:n, replace = T, size = n), ]

        fit <- GGM_regression(boot_dat,
                              IC = object$IC,
                              method = object$method,
                              rule = object$rule)

        boot_r2 <- lapply(1:p, function(x) summary(fit$regression_results[[x]]$BestModel)$r.squared)
        boot_mse <- lapply(1:p, function(x) mean(fit$regression_results[[x]]$BestModel$residuals^2))

        R2_boot_storage[i, ] <-  unlist(boot_r2)
        mse_boot_storage[i, ] <- unlist(boot_mse)
      }

      returned_object <- list(R2_boot_storage = R2_boot_storage,
                              mse_boot_storage = mse_boot_storage,
                              iter = iter,
                              p = p, n = n,
                              newdata = newdata,
                              reg_object = object,
                              resample = resample)


    } else {

      newdata <- TRUE
      message("currently only MSE is available for new data")
      colnames(new_data) <-  paste("X", 1:p, sep = "")
      new_data <- scale(new_data)


      for(i in 1:iter){
        setTxtProgressBar(pb, i)
        boot_dat <- dat[sample(1:n, replace = T, size = n), ]
        fit <- GGM_regression(boot_dat,
                              IC = object$IC,
                              method = object$method,
                              rule = object$rule)

        boot_mse <- lapply(1:p, function(x) {
          preds<- as.matrix(new_data[,names(fit$regression_results[[x]]$BestModel$coefficients )])  %*%
            fit$regression_results[[x]]$BestModel$coefficients;
          mean((preds - new_data[,x])^2)
        })
        mse_boot_storage[i,] <- unlist(boot_mse)

      }

      returned_object <- list(mse_boot_storage = mse_boot_storage,
                              iter = iter,
                              p = p, n = n,
                              newdata = newdata,
                              reg_object = object,
                              resample = resample)
    }
  } else {
    stop("currently only bootstrapping network predictability is possible (no point estimates)")
  }

  class(returned_object) <- "predict.GGM_regression"
  returned_object


}

#' @title Print method for a \code{predict.GGM_regression} object
#' @param x object of class \code{predict.GGM_regression}
#' @param ... currently ignored
#'
#' @export
#'
#' @examples
#' # data
#' X <- GGMnonreg::ptsd[, 1:5]
#'
#' # fit model
#' fit <- GGM_regression(X)
#'
#' # predict
#' net_pred <- predict(fit, iter = 10)
#' net_pred
print.predict.GGM_regression <- function(x,...){
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Network Predictability \n")
  cat(paste("Information Criterion:", x$reg_object$IC, "\n"))
  cat(paste("Rule:", x$reg_object$rule, "\n"))
  cat(paste("New Data:", x$newdata, "\n"))
  cat("----\n")
  cat(date())
}

#' @name summary.predict.GGM_regression
#' @title Summary method for a \code{GGM_regression} object
#' @param object object of class \code{GGM_regression}
#' @param ci confidence interval
#' @param ... currently ignored
#' @return data frame(s) including the summarized bootstrap samples

#' @export
#' @examples
#' # data
#' X <- GGMnonreg::ptsd[, 1:5]
#'
#' # fit model
#' fit <- GGM_regression(X)
#'
#' # predict
#' net_pred <- predict(fit, iter = 10)
#' net_pred_summ <- summary(net_pred)
#' net_pred_summ
summary.predict.GGM_regression <- function(object, ci = 0.95, ...){
  x <- object$reg_object
  lw_bound <- (1 - ci) / 2

  up_bound <- 1 - lw_bound

  if(isTRUE(object$resample)){
    if(isFALSE(object$newdata)){

      R2_mean <-  apply(object$R2_boot_storage, 2, mean)
      R2_sd <- apply(object$R2_boot_storage, 2, sd)
      R2_quantiles <- apply(object$R2_boot_storage, 2, quantile, c(lw_bound, up_bound))


      R2_dat <- data.frame(node = 1:object$p,
                           Estimate = R2_mean,
                           Est.Error =  R2_sd,
                           CI = t(R2_quantiles))

      colnames(R2_dat) <- c("Node", "Estimate", "Est.Error", "ci.lb", "ci.ub")



      mse_mean <-  apply(object$mse_boot_storage, 2, mean)
      mse_sd <- apply(object$mse_boot_storage, 2, sd)
      mse_quantiles <- apply(object$mse_boot_storage, 2, quantile, c(lw_bound, up_bound))


      mse_dat <- data.frame(node = 1:object$p,
                            Estimate = mse_mean,
                            Est.Error =  mse_sd,
                            CI = t(mse_quantiles))

      colnames(mse_dat) <- c("Node", "Estimate", "Est.Error", "ci.lb", "ci.ub")

      returned_object <- list(R2_dat = R2_dat,
                              mse_dat = mse_dat,
                              object = object)
    } else {

      mse_mean <-  apply(object$mse_boot_storage, 2, mean)
      mse_sd <- apply(object$mse_boot_storage, 2, sd)
      mse_quantiles <- apply(object$mse_boot_storage, 2, quantile, c(lw_bound, up_bound))


      mse_dat <- data.frame(node = 1:object$p,
                            Estimate = mse_mean,
                            Est.Error =  mse_sd,
                            CI = t(mse_quantiles))

      colnames(mse_dat) <- c("Node", "Estimate", "Est.Error", "ci.lb", "ci.ub")


      returned_object <- list(mse_dat = mse_dat,
                              object = object)
    } # end new data
  } # end resample
  class(returned_object) <- "summary.predict.GGM_regression"
  returned_object
}

#' @name print.summary.predict.GGM_regression
#' @title print method for object of class \code{summary.predict.GGM_regression}
#' @param x object of class \code{summary.predict.GGM_regression}
#' @param ... currently ignored
#' @export
#'
#' @examples
#' # data
#' X <- GGMnonreg::ptsd[, 1:5]
#'
#' # fit model
#' fit <- GGM_regression(X)
#'
#' # predict
#' net_pred <- predict(fit, iter = 10)
#' net_pred_summ <- summary(net_pred)
#' net_pred_summ
print.summary.predict.GGM_regression <- function(x,...){
  cat("GGMnonreg: Non-regularized GGMs \n")
  cat("Method: Multiple Regression \n")
  cat(paste("Information Criterion:", x$object$reg_object$IC, "\n"))
  cat(paste("Rule:", x$object$reg_object$rule, "\n"))
  cat(paste("New Data:", x$object$newdata, "\n"))
  cat("----\n\n")
  if(isFALSE(x$object$newdata)){
    cat("R2: Variance Explained\n")
    print(round(x$R2_dat, 3), row.names = F)
    cat("----\n\n")
    cat("MSE: Mean Squared Error\n")
    print(round(x$mse_dat, 3), row.names = F,...)
    cat("----\n")

  }  else {
    cat("MSE: Mean Squared Error\n")
    print(round(x$mse_dat, 3), row.names = F,...)
    cat("----\n")

  }
}


#' Plot \code{predict.GGM_regression} Object
#'
#' @param x object of class \code{predict.GGM_regression}
#' @param ci confidence interval
#' @param color point color
#' @param size point size
#' @param ... currently ignored
#'
#' @return \code{ggplot} object
#' @export
#'
#' @examples
#'# data
#' X <- GGMnonreg::ptsd[,1:5]
#'
#'# fit model
#' fit <- GGM_regression(X)
#'
#'# predict
#'net_pred <- predict(fit, iter = 10)
#'
#'#plot
#'plot(net_pred)
plot.predict.GGM_regression <- function(x,
                                        ci = 0.95,
                                        color = "black",
                                        size = 1,...){

  # summarise predict
  net_summ <- summary(x, ci = ci)

  # visible binding
  Node <- NA
  Estimate <- NA
  sig <- NA

  if(isFALSE(x$newdata)){

    R2_dat <-  net_summ$R2_dat[order(net_summ$R2_dat$Estimate),]
    mse_dat <-  net_summ$mse_dat[order(net_summ$mse_dat$Estimate),]

    R2_dat$Node  <- factor(R2_dat$Node,
                           labels = R2_dat$Node,
                           levels = R2_dat$Node)

    mse_dat$Node  <- factor(mse_dat$Node,
                            labels = mse_dat$Node,
                            levels = mse_dat$Node)

    R2_plt <-  ggplot(R2_dat, aes(x = Node,
                                  y = Estimate)) +
      geom_errorbar(aes(ymin = R2_dat[,4],
                        ymax = R2_dat[,5]),
                    width = 0) +
      geom_point(size = size, color = color) +
      coord_flip() +
      ylab("Variance Explained")

    mse_plt <-  ggplot(mse_dat, aes(x = Node,
                                    y = Estimate)) +
      geom_errorbar(aes(ymin = mse_dat[,4],
                        ymax = mse_dat[,5]),
                    width = 0) +
      geom_point(size = size, color = color) +
      coord_flip() +
      ylab("Mean Squared Error")

    plts <- list(R2_plt = R2_plt, mse_plt = mse_plt)

  } else {

    mse_dat <-  net_summ$mse_dat[order(net_summ$mse_dat$Estimate),]
    mse_dat$Node  <- factor(mse_dat$Node,
                            labels = mse_dat$Node,
                            levels = mse_dat$Node)

    mse_plt <-  ggplot(mse_dat, aes(x = Node,
                                    y = Estimate)) +
      geom_errorbar(aes(ymin = mse_dat[,4],
                        ymax = mse_dat[,5]),
                    width = 0) +
      geom_point(size = size, color = color) +
      coord_flip() +
      ylab("Mean Squared Error")

    plts <- list( mse_plt = mse_plt)
  }
  return(plts)
}


