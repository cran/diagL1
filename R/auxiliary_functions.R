#main author: Kévin Allan Sales Rodrigues

# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#Funções auxiliares para quantreg

#' Function to calculate robust lambda estimator
#'
#' Function to compute a robust lambda estimator, which is based on the residuals order statistics.
#'
#' @param model A linear model L1 fitted with the \code{rq} function from the \code{quantreg} package.
#' @returns Robust lambda estimator, which is based on the residuals order statistics.
#' @export
#'
#' @examples
#' require(quantreg)
#' data(stackloss)
#' model_L1 = regL1(stack.loss ~ stack.x)
#' lambda_ros(model_L1)
#'
lambda_ros = function(model){
  return(model$lambda_ros)
}

#' Function to estimate lambda via MLE
#'
#'
#' @param model A linear model L1 fitted with the \code{rq} function from the \code{quantreg} package.
#' @returns lambda MLE (Maximum Likelihood Estimator) estimator.
#' @export
#'
#' @examples
#' require(quantreg)
#' data(stackloss)
#' model_L1 = regL1(stack.loss ~ stack.x)
#' lambda_mle(model_L1)
#'
lambda_mle = function(model){
  return(model$MAE)
}


#' Function to calculate SAE
#'
#' Function to calculate SAE (Sum of Absolute Errors) of the adjusted model.
#'
#' @param model A linear model L1 fitted with the \code{rq} function from the \code{quantreg} package.
#' @returns SAE (Sum of Absolute Errors).
#' @export
#'
#' @examples
#' require(quantreg)
#' data(stackloss)
#' model_L1 = regL1(stack.loss ~ stack.x)
#' SAE(model_L1)
#'
SAE = function(model){
 return(model$SAE)
}
