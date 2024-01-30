#main author: Kévin Allan Sales Rodrigues

#### medida de influencia condicional

#' Calculate Conditional Likelihood Displacement
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
#' likelihoodCD(mod1)

likelihoodCD = function(model){

  # complete model
  y = as.numeric(model$y)
  x = model.frame(model)[,-1]
  colnames(x) <- NULL
  x = as.matrix(x)
  #beta = as.numeric(coef(model))
  SAE = model$SAE # same that sum(abs(y -as.matrix(cbind(1,x))%*%as.matrix(beta)))
  N = model$N

  LD_conditional = NULL

  for(i in 1:N){ # calculate statistics without ith obs
    model_i = regL1(y ~ x, subset = -i)
    beta_i = as.numeric(coef(model_i))

    LD_conditional_i =  2*N*log(sum(abs(y -as.matrix(cbind(1,x))%*%as.matrix(beta_i)))/
                                SAE)
    LD_conditional = c(LD_conditional, LD_conditional_i)
  }

  plot(LD_conditional, main="Conditional Likelihood Displacement", xlab="Indices",ylab="Conditional Likelihood Displacement")
  return(LD_conditional)
}


