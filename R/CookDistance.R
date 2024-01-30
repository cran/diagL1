#main author: Kévin Allan Sales Rodrigues

#### uma medida de influencia

#' Calculate Cook Distance
#'
#' @param model Object returned from regL1 representing the fit of the L1 model.
#'
#' @returns
#' \item{Cook Distance}{          A vector with Cook Distance for each observation.}.
#'
#' @references Sun, R.-B. and Wei, B.-C. (2004). On influence assessment for lad regression.
#' \emph{Statistics & Probability Letters}, \strong{67}, 97-110. \doi{10.1016/j.spl.2003.08.018}.
#'
#' @export
#' @examples
#' set.seed(123)
#' x = matrix(rnorm(100), ncol = 2)
#' y = x[, 1] + x[, 2] + rlaplace(50, 0, 5)
#'
#' # Fits a linear regression L1 model
#' mod1 = regL1(y ~ x)
#' CookDistance(mod1)

CookDistance = function(model){

  # complete model
  y = as.numeric(model$y)
  x = model.frame(model)[,-1]
  colnames(x) <- NULL
  x = as.matrix(x)
  beta = as.numeric(coef(model))
  lambda_mle = model$MAE
  N = model$N

  CD = NULL

  for(i in 1:N){ # calculate statistics without ith obs
    model_i = regL1(y ~ x, subset = -i)
    beta_i = as.numeric(coef(model_i))

    CD_i = (lambda_mle^(-2))* t(beta-beta_i) %*%(t(cbind(1,x))%*%cbind(1,x)) %*%(beta-beta_i)
    CD = c(CD, CD_i)
  }

  plot(CD, main="Cook Distance", xlab="Indices",ylab="Cook Distance")
  return(CD)
}

