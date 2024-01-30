#' Bootstrapping Linear L1 Models
#'
#' This function can be used to construct standard errors,
#' confidence intervals and tests of hypotheses regarding linear L1
#' regression models. The bootstrap method used compute xypairs bootstrap
#' for linear L1 regression.
#'
#' @param x the regression design matrix.
#' @param y the regression response vector.
#' @param groups vector with the group index associated with the observation, observations with the same index belong to the same group.
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
#'
#' bt1 = boot_het.regL1(x, y, group_index, 1000)
#' bt1$lambda_mle_het
#' }
#'
#'

boot_het.regL1 = function(x, y, groups, R = 1000){

  # Transforms NAs in zeros, useful for empty lambda indices
  replace_na_with_zero <- function(vector) {
    vector[is.na(vector)] <- 0
    return(vector)
  }

  N = length(y) # sample size
  x_add_intercept = cbind(1, x) # add intercept
  # compute xypairs bootstrap for L1 linear regression
  RR = R +200
  bootstrap = quantreg::boot.rq(x_add_intercept, y, tau = 0.5, R = RR, bsmethod = "xy")
  # Os 50 a mais no RR servem de reserva ara o caso de descarte de amostragens bootstrap

  ## Parte heterocedástica
  lambda_mle_het = NULL
  #número de grupos
  gg = length(unique(groups))

  discarded_samples = NULL
  for(k in 1:RR){

    #Loop de teste para evitar "Singular design matrix"
    # aplicando o subset à matriz x e aos vetores y e groups
    if(is.null(dim(x))){
      x_subseted = x[bootstrap$U[,k]]
    }else{
      x_subseted = x[bootstrap$U[,k],]
    }
    y_subseted = y[bootstrap$U[,k]]
    groups_subseted = groups[bootstrap$U[,k]]

    mod_het_try = try(
      suppressWarnings({
        mod_het = regL1_het(x_subseted, y_subseted, groups_subseted)
      })
      , silent = TRUE)
    if(inherits(mod_het_try, "try-error")){
      discarded_samples = c(discarded_samples, k)
    }else{
      mod_het = mod_het_try
      #Prossegue com os cálculos
      y_hat = x_add_intercept[bootstrap$U[,k],] %*% mod_het$coefficients
      residuals = y[bootstrap$U[,k]] - y_hat

      residuals_by_group = cbind(residuals,groups[bootstrap$U[,k]])

      lambda_mle_het_temp = NULL
      for( i in unique(groups)){
        lambda_mle_het_temp[i] = mean(abs(residuals_by_group[residuals_by_group[,2]==i,1]))
      }

      lambda_mle_het = rbind(lambda_mle_het, lambda_mle_het_temp)
    }

  }

  # using the R bootstrap sample requested
  lambda_mle_het = lambda_mle_het[1:R,]
  rownames(lambda_mle_het) = 1:dim(lambda_mle_het)[1]
  #Removing NAs
  lambda_mle_het = replace_na_with_zero(lambda_mle_het)

  # descartando amostras bootstrap ("Singular design matrix")
  # caso o vetor de índices descartados não seja nulo
  if(!is.null(discarded_samples)){
    bootstrap$B = bootstrap$B[-discarded_samples,]
    bootstrap$U = bootstrap$U[,-discarded_samples]
  }


  # print(discarded_samples)
  # print(dim(bootstrap$B))
  # print(dim(bootstrap$U))

  ## Parte homocedástica
  lambda_mle = rep(0, R)

  # calcula resíduos pela fórmula y -X*Beta_chapeu
  for(k in 1:R){
    y_hat = x_add_intercept[bootstrap$U[,k],] %*% bootstrap$B[k,]
    residuals = y[bootstrap$U[,k]] - y_hat

    lambda_mle[k] = sum(abs(as.numeric(residuals)))/N
  }


  #estatísticas alternativas
  lambda_mle_Ha = apply(lambda_mle_het, 1, max)
  max_lambda_ratio = lambda_mle_Ha/lambda_mle
  max_log_lambda_ratio = log(lambda_mle_Ha/lambda_mle)

  # LRT = Likelihood Ratio Test. É o teste padrão.
  test_LRT = 0
  for(k in 1:R){
    test_LRT[k] = 0
    group_temp = groups[bootstrap$U[,k]]
    for(i in (unique(group_temp))){
      test_LRT[k] = test_LRT[k] +length(group_temp[group_temp==i])*log(lambda_mle_het[k,i]/lambda_mle[k])
    }
  }
  test_LRT = -2*test_LRT

  lambda_mle_quantiles = stats::quantile(lambda_mle, probs = seq(0,1,0.01))
  test_LRT_quantiles = stats::quantile(test_LRT, probs = seq(0,1,0.01))
  max_lambda_ratio_quantiles = stats::quantile(max_lambda_ratio, probs = seq(0,1,0.01))
  max_log_lambda_ratio_quantiles = stats::quantile(max_log_lambda_ratio, probs = seq(0,1,0.01))

  object = list(test_LRT = test_LRT, test_LRT_quantiles = test_LRT_quantiles,
                max_lambda_ratio = max_lambda_ratio,
                max_lambda_ratio_quantiles = max_lambda_ratio_quantiles,
                max_log_lambda_ratio = max_log_lambda_ratio,
                max_log_lambda_ratio_quantiles = max_log_lambda_ratio_quantiles,
                lambda_mle = lambda_mle, lambda_mle_quantiles = lambda_mle_quantiles,
                lambda_mle_het = lambda_mle_het,
                B = bootstrap$B, U = bootstrap$U )
  return(object)
}

