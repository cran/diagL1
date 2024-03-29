% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forwardSearch_regL1.R
\name{forwardSearch_regL1}
\alias{forwardSearch_regL1}
\title{Forward Search in Linear L1 Models}
\usage{
forwardSearch_regL1(
  formula,
  data,
  nsamp = "best",
  intercept = TRUE,
  trace = TRUE,
  subset,
  weights,
  na.action,
  method = "br",
  model = TRUE,
  contrasts = NULL,
  ...
)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right.}

\item{data}{a data.frame in which to interpret the variables named in the formula, or in the subset and the weights argument. If this is missing, then the variables in the formula should be on the search list. This may also be a single number to handle some special cases – see below for details.}

\item{nsamp}{the initial subset for the forward search in linear regression is found by fitting the regression model with the R function \code{\link{lmsreg}}. This argument allows to control how many subsets are used in the Least Median of Squares regression. The choices are: the number of samples or "best" (the default) or "exact" or "sample". For details see \code{\link{lmsreg}}.}

\item{intercept}{logical for the inclusion of the intercept (if no formula is provided).}

\item{trace}{logical, if TRUE a message is printed for every ten iterations completed during the forward search.}

\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}

\item{weights}{vector of observation weights; if supplied, the algorithm fits to minimize the sum of the weights multiplied into the absolute residuals. The length of weights must be the same as the number of observations. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous.}

\item{na.action}{a function to filter missing data. This is applied to the model.frame after any subset argument has been used. The default (with na.fail) is to create an error if any missing values are found. A possible alternative is na.omit, which deletes observations that contain one or more missing values.}

\item{method}{the algorithmic method used to compute the fit. There are several options: "br", "fn", "pfn", "sfn", "fnc", "conquer", "pfnb", "qfnb", "ppro" and "lasso". See \code{\link{rq}} for more details.}

\item{model}{if TRUE then the model frame is returned. This is essential if one wants to call summary subsequently.}

\item{contrasts}{a list giving contrasts for some or all of the factors default = NULL appearing in the model formula. The elements of the list should have the same name as the variable and should be either a contrast matrix (specifically, any full-rank matrix with as many rows as there are levels in the factor), or else a function to compute such a matrix given the number of levels.}

\item{...}{additional arguments for the fitting routines (see \code{\link{rq.fit.br}} and \code{\link{rq.fit.fnb}}, etc. and the functions they call).}
}
\value{
A fitted forward search in linear L1 regression model object.
}
\description{
This function applies the forward search approach to robust analysis in linear L1 models. This function is based on function \code{fwdlm} of package 'forward'.
}
\examples{
\donttest{
# applies the forward search approach to robust analysis in a linear L1 model
mod = forwardSearch_regL1(Concentration ~ Age, data = bile)
}
}
\references{
Atkinson, A.C. and Riani, M. (2000). \emph{Robust Diagnostic Regression Analysis}. New York: Springer.
}
