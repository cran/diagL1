
#' Probability density function (PDF) of the Laplace distribution
#'
#' @param x Values for which to calculate the probability density.
#' @param location Location parameter of the Laplace distribution (default = 0).
#' @param scale Scale parameter of the Laplace distribution (default = 1).
#' @returns Vector of calculated probability densities.
#' @examples
#' dlaplace(0, 0, 1)
#' @export
dlaplace = function(x, location = 0, scale = 1) {
  density = 0.5 * scale * exp(-abs(x - location) / scale)
  return(density)
}

#' Cumulative distribution function (CDF) of the Laplace distribution
#'
#' @param q Values for which to calculate the cumulative distribution.
#' @param location Location parameter of the Laplace distribution (default = 0).
#' @param scale Scale parameter of the Laplace distribution (default = 1).
#' @returns Vector of calculated cumulative distributions.
#' @examples
#' plaplace(0, 0, 1)
#' @export
plaplace = function(q, location = 0, scale = 1) {
  cdf = 0.5 + 0.5 * sign(q - location) * (1 - exp(-abs(q - location) / scale))
  return(cdf)
}

#' Quantile function (inverse of the CDF) of the Laplace distribution
#'
#' @param p Values for which to calculate quantiles.
#' @param location Location parameter of the Laplace distribution (default = 0).
#' @param scale Scale parameter of the Laplace distribution (default = 1).
#' @returns Vector of calculated quantiles.
#' @examples
#' qlaplace(0.5, 0, 1)
#' @export
qlaplace = function(p, location = 0, scale = 1) {
  quantile = location - scale * sign(p - 0.5) * log(1 - 2 * abs(p - 0.5))
  return(quantile)
}

#' Function to generate random values from the Laplace distribution
#'
#' @param n Number of values to generate.
#' @param location Location parameter of the Laplace distribution (default = 0).
#' @param scale Scale parameter of the Laplace distribution (default = 1).
#' @returns Vector of generated random values.
#'
#' @importFrom stats runif
#'
#' @examples
#' rlaplace(10, 0, 1)
#' @export
rlaplace = function(n, location = 0, scale = 1) {
  u = stats::runif(n)
  values = location - scale * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5))
  return(values)
}
