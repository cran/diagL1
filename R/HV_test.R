#main author: Kévin Allan Sales Rodrigues

#' Homogeneity of Variance Tests for Linear L1 Models
#'
#' @param y A vector with response variables.
#' @param x A matrix with a single explanatory variable.
#' @param groups Vector containing the group index to which the observation belongs.
#' @param alpha Significance level of the test, must be between 0 and 1.
#' @param tolerance threshold that determines when the iterative algorithm should stop (for \code{regL1_het}).
#' @param max_iteration maximum number of iterations (for \code{regL1_het}).
#'
#' @returns
#' A list with results from 3 homogeneity of variance tests
#' \item{alpha}{alpha argument.}
#' \item{asymptotic_critical_value}{asymptotic alpha-based test critical value.}
#' \item{HV_LRT}{LF1 statistic value using MLE (Maximum Likelihood Estimator).}
#' \item{HV_max_lambda_ratio}{p-value of LF1 statistic using MLE.}
#' \item{HV_max_log_lambda_ratio}{LF1 statistic value using ROS (Residuals Order Statistics).}
#' \item{number_of_groups}{number of groups.}
#' \item{N}{overall sample size.}
# #' \item{beta_heteroscedastic}{heteroscedastic estimate of beta.}
#' \item{lambda_heteroscedastic}{heteroscedastic estimation of the scaling parameters.}
#' \item{lambda_homocedastic}{homoscedastic estimation of the scale parameter.}
#'
#' @importFrom stats qchisq
#' @export
#' @details The 3 statistics to test homogeneity of variance are discussed in Rodrigues, Elian and Pereira (2025), for more details see this reference. In practice, use the HV_LRT statistic results. If possible, use the HV_LRT statistic with a critical value obtained via bootstrap simulation.
#'
#' @references Rodrigues, K. A. S., Elian, S. N., & Pereira, G. H. A. (2025). Homoscedasticity tests for L1 regression and their performance evaluation through simulations. \emph{Statistics}, Advance online publication. https://doi.org/10.1080/02331888.2025.2536097
#'
#' Rodrigues, K. A. S., & Elian, S. N. (2025). Influence measures for L1 regression: an analysis with the R package diagL1. \emph{Journal of Applied Statistics}, Advance online publication. https://doi.org/10.1080/02664763.2025.2510691
#'
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x1 = matrix(rnorm(20), ncol = 1)
#' y1 = x1 + rlaplace(20, 0, 1)
#' x2 = matrix(rnorm(20), ncol = 1)
#' y2 = x2 + rlaplace(20, 0, 1.5)
#' x3 = matrix(rnorm(20), ncol = 1)
#' y3 = x3 + rlaplace(20, 0, 2)
#' x4 = matrix(rnorm(20), ncol = 1)
#' y4 = x4 + rlaplace(20, 0, 2.5)
#' x5 = matrix(rnorm(20), ncol = 1)
#' y5 = x5 + rlaplace(20, 0, 3)
#'
#' y = c(y1, y2, y3, y4, y5)
#' x = rbind(x1, x2, x3, x4, x5)
#' group_index = c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))
#'
#' # Application of the homogeneity of variance test
#' test_result = HV_test(y, x, group_index)
#' test_result
#' }
#'

HV_test = function(y, x, groups, alpha = 0.05, tolerance = 1e-3, max_iteration = 2000){
  mod_ha = regL1_het(x, y, groups = groups, tolerance = tolerance, max_iteration = max_iteration)
  residuos = cbind(abs(as.numeric(mod_ha$residuals)),groups)

  mod_h0 = regL1(y ~ as.matrix(x), na.action = na.omit, method="br", model = TRUE)

  lambda0 = mod_h0$MAE

  # LRT = Likelihood Ratio Test. (It's the default test.)
  test_LRT = 0
  gg = length(unique(groups)) #number of groups
  for (i in 1:gg){
    test_LRT = test_LRT +length(groups[groups==i])*log(mod_ha$lambda[i]/lambda0)
  }
  test_LRT = -2*test_LRT

  #alternative statistics
  max_lambda_ratio = max(mod_ha$lambda/lambda0)
  max_log_lambda_ratio = log(max(mod_ha$lambda/lambda0))

  #heteroscedastic estimates
  beta_final   = mod_ha$beta
  lambda_final = mod_ha$lambda

  HV = NULL
  HV$alpha = alpha
  HV$asymptotic_critical_value = stats::qchisq(1- alpha, gg -1, ncp=0)
  HV$HV_LRT = test_LRT
  HV$HV_max_lambda_ratio = max_lambda_ratio
  HV$HV_max_log_lambda_ratio = max_log_lambda_ratio
  HV$number_of_groups = gg
  HV$N = length(y)
  #HV$beta_heteroscedastic = beta_final
  HV$lambda_heteroscedastic = lambda_final
  HV$lambda_homocedastic = lambda0

  return(HV)
}


