#main author: Kévin Allan Sales Rodrigues

#' Lack of Fit Tests for Linear L1 Models
#'
#' @param y A vector with response variables.
#' @param x A matrix with a single explanatory variable.
#' @param groups Vector containing the group index to which the observation belongs.
#' @param alpha Significance level of the test, must be between 0 and 1.
#'
#' @returns
#' A list with results from 3 lack of fit tests
#' \item{alpha}{alpha argument.}
#' \item{critical_value}{alpha-based test critical value.}
#' \item{LF1_MLE}{LF1 statistic value using MLE (Maximum Likelihood Estimator).}
#' \item{LF1_MLE}{p-value of LF1 statistic using MLE.}
#' \item{LF1_ROS}{LF1 statistic value using ROS (Residuals Order Statistics).}
#' \item{LF2}{LF2 statistic value.}
#' \item{modelo_H0}{model fitted under H0.}
#' \item{modelo_Ha}{model fitted under Ha.}
#' \item{MLE}{estimation of the scale parameter of the estimator model via MLE.}
#' \item{ROS}{estimation of the scale parameter of the estimator model via ROS.}
#' \item{SAE_H0}{SAE (Sum of Absolute Errors) of the adjusted model under H0.}
#' \item{SAE_Ha}{SAE (Sum of Absolute Errors) of the adjusted model under Ha.}
#' \item{matrix_mean_x}{average of the explanatory variable per group of observations.}
#' \item{number_of_groups}{number of groups.}
#'
#' @importFrom stats qchisq
#' @export
#' @details The 3 statistics to test lack of fit are discussed in Rodrigues (2024), for more details see this reference. In practice, use the LF1_MLE statistic results. These tests were developed with just one explanatory variable in mind, which is why we include an error if there is more than one explanatory variable.
#'
#' @references Rodrigues, K. A. S. (2024). \strong{Analysis of the adjustment of the L1 regression model}.
#' Phd dissertation, University of São Paulo, BR.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x1 = matrix(rnorm(20), ncol = 1)
#' y1 = x1 + rlaplace(20, 0, 5)
#' x2 = matrix(rnorm(20), ncol = 1)
#' y2 = x2 + rlaplace(20, 1, 5)
#' x3 = matrix(rnorm(20), ncol = 1)
#' y3 = x3 + rlaplace(20, 2, 5)
#' x4 = matrix(rnorm(20), ncol = 1)
#' y4 = x4 + rlaplace(20, 3, 5)
#' x5 = matrix(rnorm(20), ncol = 1)
#' y5 = x5 + rlaplace(20, 4, 5)
#'
#' y = c(y1, y2, y3, y4, y5)
#' x = rbind(x1, x2, x3, x4, x5)
#' group_index = c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))
#'
#' # Application of the lack of fit test
#' test_result = LF_test(y, x, group_index)
#' test_result
#' }
#'

LF_test = function(y, x, groups, alpha = 0.05){

  if (!is.null(dim(x)) && dim(x)[2] > 1){
    stop("Currently this test only works for a single explanatory variable.")
  }

  x = as.matrix(x)
  explanatory_variables_number = 1 #dim(x)[2]
  N = length(y) #sample size
  group_index = unique(groups)

  # matriz com médias das variáveis explicativas por cluster
  matrix_mean_x =  matrix(, nrow = 0, ncol = explanatory_variables_number + 1)

  for(i in group_index){

    matrix_mean_x = rbind(matrix_mean_x,
                          c(i, as.numeric(mean(x[groups == i])))
                          #c(i, as.numeric(colMeans(x[groups == i,]))) # caso com mais de 1 variável explicativa
    )

  }

  #colocando média das variáveis explicativas dos clusters nas observações
  x_new = x

  for(i in 1:N){

    # trocando os valores das variáveis explicativas da i-ésima observação pela média das variáveis explicativas do cluster em que a observação se encontra
    x_new[i] = matrix_mean_x[groups[i], 2:(explanatory_variables_number+1)]
    #x_new[i, ] = matrix_mean_x[groups[i], 2:(explanatory_variables_number+1)] # caso com mais de 1 variável explicativa
  }
  # o número de colunas é o número de variáveis + 1 porque uma das colunas serve para registrar o índice do grupo

  fit_linear_H0 = regL1(y ~ as.matrix(x_new), na.action = na.omit, method="br", model = TRUE)
  fit_linear_HA = regL1(y ~ as.factor(groups), na.action = na.omit, method="br", model = TRUE)

  #number of groups
  k = length(group_index)

  LF1_linear_MLE = (2/lambda_mle(fit_linear_H0))*(SAE(fit_linear_H0) - SAE(fit_linear_HA))
  LF1_linear_ROS = (2/lambda_ros(fit_linear_H0))*(SAE(fit_linear_H0) - SAE(fit_linear_HA))
  LF2_linear = 2*(SAE(fit_linear_H0) - SAE(fit_linear_HA))/(SAE(fit_linear_HA)/(N-k))

  LF = NULL
  LF$alpha = alpha
  LF$critical_value = stats::qchisq(1- alpha, k -explanatory_variables_number -1, ncp=0)
  LF$LF1_MLE = LF1_linear_MLE
  LF$LF1_MLE_p_value = 1 -stats::pchisq(LF1_linear_MLE, k -explanatory_variables_number -1, ncp=0)
  LF$LF1_ROS = LF1_linear_ROS
  LF$LF2 = LF2_linear
  LF$modelo_H0 = fit_linear_H0
  LF$modelo_Ha = fit_linear_HA
  LF$MLE = lambda_mle(fit_linear_H0)
  LF$ROS = lambda_ros(fit_linear_H0)
  LF$SAE_H0 = SAE(fit_linear_H0)
  LF$SAE_Ha = SAE(fit_linear_HA)
  LF$matrix_mean_x = matrix_mean_x
  LF$number_of_groups = k

  return(LF)
}


