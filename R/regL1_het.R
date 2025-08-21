#' Fitting Heteroscedastic Linear L1 Models
#'
#' This function fits an groupwise heteroscedastic L1 regression model using the \code{\link[quantreg]{rq}} function from the 'quantreg' package.
#'
#'
#' @param x the regression design matrix.
#' @param y the regression response vector.
#' @param groups vector with the group index associated with the observation, observations with the same index belong to the same group.
#' @param tolerance threshold that determines when the iterative algorithm should stop.
#' @param max_iteration maximum number of iterations.
#' @param na.action a function to filter missing data. This is applied to the model.frame after any subset argument has been used. The default (with na.fail) is to create an error if any missing values are found. A possible alternative is na.omit, which deletes observations that contain one or more missing values.
#' @param method the algorithmic method used to compute the fit. There are several options: "br", "fn", "pfn", "sfn", "fnc", "conquer", "pfnb", "qfnb", "ppro" and "lasso". See \code{\link[quantreg]{rq}} for more details.
#' @param model if TRUE then the model frame is returned. This is essential if one wants to call summary subsequently.
#' @param ... additional arguments for the fitting routines (see \code{\link[quantreg]{rq.fit.br}} and \code{\link[quantreg]{rq.fit.fnb}}, etc. and the functions they call).
#'
#' @returns A fitted heteroscedastic L1 linear regression model object.
#'
#' @importFrom stats na.omit
#' @import quantreg
#' @importFrom MatrixModels model.Matrix
#' @import Matrix
#'
#' @export
#'
#' @details L1 regression is an important particular case of quantile regression, so this function inherits from the "rq" class of the \code{quantreg} package.
#'
#' @examples
#' set.seed(123)
#' x1 = matrix(rnorm(20), ncol = 2)
#' y1 = x1[, 1] + x1[, 2] + rlaplace(10, 0, 5)
#' x2 = matrix(rnorm(20), ncol = 2)
#' y2 = x2[, 1] + x2[, 2] + rlaplace(10, 0, 10)
#' x3 = matrix(rnorm(20), ncol = 2)
#' y3 = x3[, 1] + x3[, 2] + rlaplace(10, 0, 15)
#' x4 = matrix(rnorm(20), ncol = 2)
#' y4 = x4[, 1] + x4[, 2] + rlaplace(10, 0, 20)
#' x5 = matrix(rnorm(20), ncol = 2)
#' y5 = x5[, 1] + x5[, 2] + rlaplace(10, 0, 30)
#'
#' y = c(y1, y2, y3, y4, y5)
#' x = rbind(x1, x2, x3, x4, x5)
#' group_index = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10))
#' # Fits a heteroscedastic linear regression L1 model
#' mod1 = regL1_het(x, y, group_index)
#'
#@param formula a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right.

regL1_het = function(x, y, groups, na.action = stats::na.omit, method="br", model = TRUE, tolerance = 1e-3, max_iteration = 2000, ...)
{

  # Transforms NAs in zeros, useful for empty lambda indices
  replace_na_with_zero <- function(vector) {
    vector[is.na(vector)] <- 0
    return(vector)
  }

  mod = regL1(y ~ as.matrix(x), na.action = stats::na.omit, method="br", model = TRUE)

  #print("erro no primeiro modelo")

  #groups é o vetor com o índice do grupo associado à observação, observações com o mesmo índice pertencem ao mesmo grupo

  #número de grupos
  gg = length(unique(groups))

  #vetor de betas (inclusive intercepto)
  beta_old = as.numeric(mod$coefficients)
  residuos = cbind(as.numeric(mod$residuals),groups)

  #lambda_old = rep(mean(abs(residuos[,1])), gg) #lambda homoscedástico
  lambda_old = NULL
  for( i in unique(groups)){ # se der problema, trocar "unique(groups)" por "1:gg"
    lambda_old[i] = mean(abs(residuos[residuos[,2]==i,1]))
  }
  lambda_new = lambda_old
  beta_new = beta_old
  weighted_lambda = 1/lambda_new[residuos[,2]]

  #contador de iteração
  iteration = 1

  #Já calculei beta e lambda de 1º passo
  #Preciso atualizar beta

  #weighted L1 regression
  mod2 = regL1(y ~ x, weights=weighted_lambda,ci=FALSE)

  beta_new = as.numeric(mod2$coefficients)
  residuos = cbind(as.numeric(mod2$residuals),groups)

  for( i in unique(groups)){
    lambda_new[i] = mean(abs(residuos[residuos[,2]==i,1]))
  }
  #pesos do novo modelo
  weighted_lambda = 1/lambda_new[residuos[,2]]

  #contador de iteração
  iteration = 2

  #vetor dos valores intermediários de beta e lambda
  beta_list = rbind(beta_old,beta_new)
  lambda_list = rbind(lambda_old,lambda_new)

  #Removing NAs
  lambda_old = replace_na_with_zero(lambda_old)
  lambda_new = replace_na_with_zero(lambda_new)


  while( sum(abs(beta_new -beta_old)) +sum(abs(lambda_new -lambda_old)) > tolerance && iteration < max_iteration ){
    lambda_old = lambda_new
    beta_old = beta_new

    #Já calculei beta e lambda de 2º passo
    #Preciso atualizar beta

    #weighted L1 regression
    mod2 = regL1(y ~ x, weights=weighted_lambda,ci=FALSE)

    beta_new = as.numeric(mod2$coefficients)
    residuos = cbind(as.numeric(mod2$residuals),groups)

    for( i in unique(groups)){
      lambda_new[i] = mean(abs(residuos[residuos[,2]==i,1]))
    }
    #pesos do novo modelo
    weighted_lambda = 1/lambda_new[residuos[,2]]

    beta_list = rbind(beta_list,beta_new)
    lambda_list = rbind(lambda_list,lambda_new)

    iteration = iteration +1
  }
  beta_final = beta_new
  lambda_final = lambda_new

  FIT = list(coefficients = beta_final, x = mod2$x, y = mod2$y,
             residuals = residuos, dual = mod2$dual,
             fitted.values = mod2$fitted.values, formula = mod2$formula,
             terms = mod2$terms, xlevels = mod2$xlevels, call = mod2$call,
             tau = mod2$tau, rho = mod2$rho, method = mod2$method,
             model = mod2$model, call_rq = mod2$call_rq, SAE = mod2$SAE,
             N = mod2$N, lambda = lambda_final, groups = groups,
             iteration = iteration, beta_list = beta_list,
             lambda_list = lambda_list)
  class(FIT) = "regL1" # troquei a classe regL1
  return(FIT)
}

