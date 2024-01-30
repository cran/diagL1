#' Bootstrapping Linear L1 Models
#'
#' This function can be used to construct standard errors,
#' confidence intervals and tests of hypotheses regarding linear L1
#' regression models. The bootstrap method used compute xypairs bootstrap
#' for linear L1 regression.
#'
#' @param x the regression design matrix.
#' @param y the regression response vector.
#' @param R the number of bootstrap replications.
#' @returns A list consisting of four elements: a vector of lambda mle estimate for each bootstrap sample and another vector with empirical quantiles. A matrix B of dimension R by p is returned with the R resampled estimates of the vector of L1 linear regression parameters. A matrix U of sampled indices.
#' @importFrom quantreg boot.rq
#' @importFrom stats quantile
#' @export
#'
#' @seealso
#' \code{\link[quantreg]{boot.rq}} Bootstrapping Quantile Regression from package \code{quantreg}.
#'
#'
#' @examples
#' \donttest{
#' data(stackloss)
#' bt1 = boot.regL1(stack.x, stack.loss, 1000)
#' plot(bt1$lambda_mle)
#' }
#'
#'

boot.regL1 = function(x, y, R = 1000){

  N = length(y) # sample size
  x_add_intercept = cbind(1, x) # add intercept
  # compute xypairs bootstrap for L1 linear regression
  bootstrap = quantreg::boot.rq(x_add_intercept, y, tau = 0.5, R = R, bsmethod = "xy")

  lambda_mle = rep(0, R)

  # calcula resíduos pela fórmula y -X*Beta_chapeu
  for(k in 1:R){
    y_hat = x_add_intercept[bootstrap$U[,k],] %*% bootstrap$B[k,]
    residuals = y[bootstrap$U[,k]] - y_hat

    lambda_mle[k] = sum(abs(as.numeric(residuals)))/N
  }

  lambda_quantiles = stats::quantile(lambda_mle, probs = seq(0,1,0.01))

  object = list(lambda_mle = lambda_mle, lambda_quantiles = lambda_quantiles, B = bootstrap$B, U = bootstrap$U )
  return(object)
}

