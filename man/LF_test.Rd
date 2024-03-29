% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LF_test.R
\name{LF_test}
\alias{LF_test}
\title{Lack of Fit Tests for Linear L1 Models}
\usage{
LF_test(y, x, groups, alpha = 0.05)
}
\arguments{
\item{y}{A vector with response variables.}

\item{x}{A matrix with a single explanatory variable.}

\item{groups}{Vector containing the group index to which the observation belongs.}

\item{alpha}{Significance level of the test, must be between 0 and 1.}
}
\value{
A list with results from 3 lack of fit tests
\item{alpha}{alpha argument.}
\item{critical_value}{alpha-based test critical value.}
\item{LF1_MLE}{LF1 statistic value using MLE (Maximum Likelihood Estimator).}
\item{LF1_MLE}{p-value of LF1 statistic using MLE.}
\item{LF1_ROS}{LF1 statistic value using ROS (Residuals Order Statistics).}
\item{LF2}{LF2 statistic value.}
\item{modelo_H0}{model fitted under H0.}
\item{modelo_Ha}{model fitted under Ha.}
\item{MLE}{estimation of the scale parameter of the estimator model via MLE.}
\item{ROS}{estimation of the scale parameter of the estimator model via ROS.}
\item{SAE_H0}{SAE (Sum of Absolute Errors) of the adjusted model under H0.}
\item{SAE_Ha}{SAE (Sum of Absolute Errors) of the adjusted model under Ha.}
\item{matrix_mean_x}{average of the explanatory variable per group of observations.}
\item{number_of_groups}{number of groups.}
}
\description{
Lack of Fit Tests for Linear L1 Models
}
\details{
The 3 statistics to test lack of fit are discussed in Rodrigues (2024), for more details see this reference. In practice, use the LF1_MLE statistic results. These tests were developed with just one explanatory variable in mind, which is why we include an error if there is more than one explanatory variable.
}
\examples{
\donttest{
set.seed(123)
x1 = matrix(rnorm(20), ncol = 1)
y1 = x1 + rlaplace(20, 0, 5)
x2 = matrix(rnorm(20), ncol = 1)
y2 = x2 + rlaplace(20, 1, 5)
x3 = matrix(rnorm(20), ncol = 1)
y3 = x3 + rlaplace(20, 2, 5)
x4 = matrix(rnorm(20), ncol = 1)
y4 = x4 + rlaplace(20, 3, 5)
x5 = matrix(rnorm(20), ncol = 1)
y5 = x5 + rlaplace(20, 4, 5)

y = c(y1, y2, y3, y4, y5)
x = rbind(x1, x2, x3, x4, x5)
group_index = c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))

# Application of the lack of fit test
test_result = LF_test(y, x, group_index)
test_result
}

}
\references{
Rodrigues, K. A. S. (2024). \strong{Analysis of the adjustment of the L1 regression model}.
Phd dissertation, University of São Paulo, BR.
}
