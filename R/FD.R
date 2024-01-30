#main author: Kévin Allan Sales Rodrigues

#### uma medida de influencia

#' F Distance
#'
#' @param model Object returned from regL1 representing the fit of the L1 model.
#' @param norm Type of norm, there are two types, 1 (norm L1) and 2 (norm L2). The default is norm L2.
#'
#' @returns
#' \item{F Distance}{          A vector with F Distance for each observation.}.
#'
#' @references Sun, R.-B. and Wei, B.-C. (2004). On influence assessment for lad regression.
#' \emph{Statistics & Probability Letters}, \strong{67}, 97-110. \doi{10.1016/j.spl.2003.08.018}.
#'
#'
#' @export
#' @examples
#' set.seed(123)
#' x = matrix(rnorm(100), ncol = 2)
#' y = x[, 1] + x[, 2] + rlaplace(50, 0, 5)
#'
#' # Fits a linear regression L1 model
#' mod1 = regL1(y ~ x)
#' FD(mod1)

FD = function(model, norm = 2){

  # complete model
  y = as.numeric(model$y)
  x = model.frame(model)[,-1]
  colnames(x) <- NULL
  x = as.matrix(x)
  beta = as.numeric(coef(model))
  lambda_mle = model$MAE
  N = model$N

  F_distance = NULL

  if(norm != 1 && norm != 2){
    norm = 2
    warning("Invalid norm, we will continue using the L2 norm.")
  }

  if(norm == 2){
    for(i in 1:N){ # calculate statistics without ith obs
      model_i = regL1(y ~ x, subset = -i)
      beta_i = as.numeric(coef(model_i))

      F_distance_i = lambda_mle^(-2)*sqrt(sum((cbind(1,x)%*%beta -cbind(1,x)%*%beta_i)^(2)))
      F_distance = c(F_distance, F_distance_i)
    }
  }else{
    for(i in 1:N){ # calculate statistics without ith obs
      model_i = regL1(y ~ x, subset = -i)
      beta_i = as.numeric(coef(model_i))

      F_distance_i = lambda_mle^(-1)*sum(abs(cbind(1,x)%*%beta -cbind(1,x)%*%beta_i))
      F_distance = c(F_distance, F_distance_i)
    }
  }

  plot(F_distance, main="F Distance", xlab="Indices",
       ylab= paste('F Distance (norm L', norm, ')', sep = ""))
  return(F_distance)
}

