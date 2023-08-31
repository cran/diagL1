#main author: Kévin Allan Sales Rodrigues

#### medida de influencia

#' Calculate Likelihood Displacement
#'
#' @param model Object returned from regL1 representing the fit of the L1 model.
#'
#' @returns
#' \item{Likelihood Displacement}{          A vector with Likelihood Displacement for each observation.}.
#'
#' @references Elian, S. N., André, C. D. S. and Narula, S. C. (2000). Influence Measure for the
#' \ifelse{html}{\out{L<sub>1</sub>}}{\eqn{L1}} regression.
#' \emph{Communications in Statistics - Theory and Methods}, \strong{29}(4), 837-849. \doi{10.1080/03610920008832518}.
#'
#' @importFrom stats model.frame
#' @export
#' @examples
#' set.seed(123)
#' x = matrix(rnorm(100), ncol = 2)
#' y = x[, 1] + x[, 2] + rlaplace(50, 0, 5)
#'
#' # Fits a linear regression L1 model
#' mod1 = regL1(y ~ x)
#' likelihoodD(mod1)

likelihoodD = function(model){

  # complete model
  y = as.numeric(model$y)
  x = model.frame(model)[,-1]
  colnames(x) <- NULL
  beta = as.numeric(coef(model))
  lambda_mle = model$MAE
  lambda_ros = model$lambda_ros
  N = model$N

  LD = NULL

  for(i in 1:N){ # calculate statistics without ith obs
    model_i = regL1(y ~ x, subset = -i)
    beta_i = as.numeric(coef(model_i))
    lambda_mle_i = model_i$MAE
    lambda_ros_i = model_i$lambda_ros

    LD_i = 2*(N*log(lambda_mle_i/lambda_mle) + abs(y[i] -t(as.matrix(cbind(1,x)[i,]))%*%as.matrix(beta_i))/lambda_mle_i -1)
    LD = c(LD, LD_i)
  }

  plot(LD, main="Likelihood Displacement", xlab="Indices",ylab="Likelihood Displacement")
  return(LD)
}

