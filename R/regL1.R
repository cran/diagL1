#main author: Kévin Allan Sales Rodrigues

#' Fitting Linear L1 Models
#'
#' This function fits an L1 regression model using the \code{\link{rq}} function from the 'quantreg' package. L1 regression allows dealing with outliers and non-normal distributions in the data.
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right.
#' @param data a data.frame in which to interpret the variables named in the formula, or in the subset and the weights argument. If this is missing, then the variables in the formula should be on the search list. This may also be a single number to handle some special cases – see below for details.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights vector of observation weights; if supplied, the algorithm fits to minimize the sum of the weights multiplied into the absolute residuals. The length of weights must be the same as the number of observations. The weights must be nonnegative and it is strongly recommended that they be strictly positive, since zero weights are ambiguous.
#' @param na.action a function to filter missing data. This is applied to the model.frame after any subset argument has been used. The default (with na.fail) is to create an error if any missing values are found. A possible alternative is na.omit, which deletes observations that contain one or more missing values.
#' @param method the algorithmic method used to compute the fit. There are several options: "br", "fn", "pfn", "sfn", "fnc", "conquer", "pfnb", "qfnb", "ppro" and "lasso". See \code{\link{rq}} for more details.
#' @param model if TRUE then the model frame is returned. This is essential if one wants to call summary subsequently.
#' @param contrasts a list giving contrasts for some or all of the factors default = NULL appearing in the model formula. The elements of the list should have the same name as the variable and should be either a contrast matrix (specifically, any full-rank matrix with as many rows as there are levels in the factor), or else a function to compute such a matrix given the number of levels.
#' @param ... additional arguments for the fitting routines (see \code{\link{rq.fit.br}} and \code{\link{rq.fit.fnb}}, etc. and the functions they call).
#'
#' @returns A fitted L1 linear regression model object.
#'
#' @import stats
#' @import quantreg
#' @export
#'
#' @details L1 regression is an important particular case of quantile regression, so this function inherits from the "rq" class of the \code{quantreg} package.
#'
#' @examples
#' set.seed(123)
#' x = matrix(rnorm(100), ncol = 2)
#' y = x[, 1] + x[, 2] + rlaplace(50, 0, 5)
#'
#' # Fits a linear regression L1 model
#' mod1 = regL1(y ~ x)
#'
# @importFrom stats as.formula
# @importFrom stats na.omit
# Define a constructor function for the "regL1" class
regL1 = function(formula, data, subset = NULL, weights, na.action = na.omit,
                 method = "br", model = TRUE, contrasts = NULL, ...) {
  tau = 0.5 # L1 reg
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  if (method == "model.frame")
    return(mf)
  mt <- attr(mf, "terms")
  weights <- as.vector(model.weights(mf))
  tau <- sort(unique(tau))
  eps <- .Machine$double.eps^(2/3)
  if (any(tau == 0))
    tau[tau == 0] <- eps
  if (any(tau == 1))
    tau[tau == 1] <- 1 - eps

  # Fit the model using rq() with the provided formula
  model_formula = stats::as.formula(match.call()$formula)
  # Ajuste do modelo de regressão usando rq() com tau = 0.5
  model = quantreg::rq(formula = model_formula, data = data, subset, weights, tau=.5,
                       na.action = na.omit, method="br", model = TRUE,
                        contrasts = NULL, ...)

  model$call_rq = model$call
  model$call = call
  model$SAE = sum(abs(model$residuals))
  model$N = length(model$y)
  model$MAE = model$SAE / model$N # lambda_mle

  # lambda_ros calculations | begin
  res = cbind(sort(as.numeric(model$residuals)))
  res1 = matrix(0, nrow(res) -ncol(model$x),1)
  n = nrow(model$x) - ncol(model$x)
  m = (n+1)/2 -sqrt(n)
  j=1
  for(i in 1:nrow(res)) if(res[i] !=0) ((res1[j] = res[i]) & (j = j+1))
  pos1 = n-m +1
  pos2 = m
  e1 = round(pos1)
  e2 = round(pos2)
  if (e1 > n) (e1 = n)
  if(e2 == 0) (e2 = 1)
  lambda_robust = sqrt(n)*(res[e1] -res1[e2])/4
  # lambda_ros calculations | end

  model$lambda_ros = lambda_robust

  #model = structure(model, class = c('regL1', 'rq')) # alterando a classe
  class(model) = "regL1"

  # Create an object of class "regL1"
  #object = list(model = model)
  #class(object) = "regL1"

  return(model)
}


#' Print Linear L1 Regression Summary Object
#'
#' Print summary of linear L1 regression object.
#'
#' @param x This is an object of class "\code{summary.regL1}" produced by a call to \code{summary.regL1()}.
#' @param digits Significant digits reported in the printed table.
#' @param ... Optional arguments.
#' @returns No return value, called for side effects.
#'
#' @import greekLetters
#' @export
print.summary.regL1 = function(x, digits = max(5, .Options$digits - 2), ...)
  {
    cat("\nCall: ")
    dput(x$call)
    coef <- x$coef
    ## df <- x$df
    ## rdf <- x$rdf
    #tau <- x$tau
    #cat("\ntau: ")
    #print(format(round(tau,digits = digits)), quote = FALSE, ...)

    # important statistics
    SAE = x$SAE
    N = x$N
    MAE = x$MAE
    lambda_ROS = x$lambda_ros
    cat("\n", greekLetters::greeks('lambda'),"(MLE) | MAE: ", format(round(MAE,digits = digits)), "    SAE: ", format(round(SAE,digits = digits)), "     N: ", N, sep = "")
    cat("\n", greekLetters::greeks('lambda'),"(ROS)      : ", format(round(lambda_ROS,digits = digits)), sep ="")

    cat("\n\nCoefficients:\n")
    print(format(round(coef, digits = digits)), quote = FALSE, ...)
    invisible(x)
  }


#' Summary methods for L1 Regression
#'
#' Returns a summary list for a L1 regression fit. A null value will be returned if printing is invoked
#'
#' @param object Object returned from regL1 representing the fit of the L1 model.
#' @param se specifies the method used to compute standard standard errors. There are currently seven available methods: "rank", "iid", "nid", "ker", "boot", "BLB", "conquer" and "extreme".
#' @param covariance logical flag to indicate whether the full covariance matrix of the estimated parameters should be returned.
#' @param hs Use Hall Sheather bandwidth for sparsity estimation If false revert to Bofinger bandwidth.
#' @param U Resampling indices or gradient evaluations used for bootstrap, see \code{\link{summary.rq}} and \code{\link{boot.rq}} for more details.
#' @param gamma parameter controlling the effective sample size of the'bag of little bootstrap samples that will be b = n^gamma where n is the sample size of the original model.
#' @param ... Optional arguments. See \code{\link{summary.rq}} for more details.
#'
#' @returns No return value, called for side effects.
#' @import quantreg
#' @importFrom utils modifyList
#' @importFrom conquer conquer
#' @export
#'
#' @seealso
#' \code{\link{regL1}} for fitting linear L1 models.
#' \code{\link{summary.rq}} summary methods for Quantile Regression.
#'
# Define a summary method for the "regL1" class
summary.regL1 = function(object, se = NULL, covariance = FALSE, hs = TRUE, U = NULL, gamma = 0.7, ...){
  class(object) = "rq" # voltando a classe do objeto
  if (object$method == "lasso")
    stop("no inference for lasso'd rq fitting: try rqss (if brave, or credulous)")
  if (object$method == "conquer")
    se = "conquer"
  mt <- terms(object)
  m <- model.frame(object)
  y <- model.response(m)
  dots <- list(...)
  method <- object$method
  if (object$method == "sfn") {
    x <- object$model$x
    vnames <- names(object$coef)
    ctrl <- object$control
  }
  else {
    x <- model.matrix(mt, m, contrasts = object$contrasts)
    vnames <- dimnames(x)[[2]]
  }
  wt <- as.vector(model.weights(object$model))
  tau <- object$tau
  eps <- .Machine$double.eps^(1/2)
  coef <- coefficients(object)
  if (is.matrix(coef))
    coef <- coef[, 1]
  resid <- object$residuals
  n <- length(y)
  p <- length(coef)
  rdf <- n - p
  if (!is.null(wt)) {
    resid <- resid * wt
    x <- x * wt
    y <- y * wt
  }
  if (is.null(se)) {
    if (n < 1001 & covariance == FALSE)
      se <- "rank"
    else se <- "nid"
  }
  if (se == "rank") {
    f <- rq.fit.br(x, y, tau = tau, ci = TRUE, ...)
  }
  if (se == "iid") {
    xxinv <- diag(p)
    xxinv <- backsolve(qr(x)$qr[1:p, 1:p, drop = FALSE],
                       xxinv)
    xxinv <- xxinv %*% t(xxinv)
    pz <- sum(abs(resid) < eps)
    h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = hs)))
    ir <- (pz + 1):(h + pz + 1)
    ord.resid <- sort(resid[order(abs(resid))][ir])
    xt <- ir/(n - p)
    sparsity <- rq(ord.resid ~ xt)$coef[2]
    cov <- sparsity^2 * xxinv * tau * (1 - tau)
    scale <- 1/sparsity
    serr <- sqrt(diag(cov))
  }
  else if (se == "nid") {
    h <- bandwidth.rq(tau, n, hs = hs)
    while ((tau - h < 0) || (tau + h > 1)) h <- h/2
    bhi <- rq.fit(x, y, tau = tau + h, method = method)$coef
    blo <- rq.fit(x, y, tau = tau - h, method = method)$coef
    dyhat <- x %*% (bhi - blo)
    if (any(dyhat <= 0))
      warning(paste(sum(dyhat <= 0), "non-positive fis"))
    f <- pmax(0, (2 * h)/(dyhat - eps))
    fxxinv <- diag(p)
    if (method == "sfn") {
      D <- t(x) %*% (f * x)
      D <- chol(0.5 * (D + t(D)), nsubmax = ctrl$nsubmax,
                nnzlmax = ctrl$nnzlmax, tmpmax = ctrl$tmpmax)
      fxxinv <- backsolve(D, fxxinv)
    }
    else {
      fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p,
                                             drop = FALSE], fxxinv)
      fxxinv <- fxxinv %*% t(fxxinv)
    }
    xx <- t(x) %*% x
    cov <- tau * (1 - tau) * fxxinv %*% xx %*% fxxinv
    scale <- mean(f)
    serr <- sqrt(diag(cov))
  }
  else if (se == "ker") {
    h <- bandwidth.rq(tau, n, hs = hs)
    while ((tau - h < 0) || (tau + h > 1)) h <- h/2
    uhat <- c(y - x %*% coef)
    h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)),
                                                 (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
    f <- dnorm(uhat/h)/h
    fxxinv <- diag(p)
    fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p, drop = FALSE],
                        fxxinv)
    fxxinv <- fxxinv %*% t(fxxinv)
    cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*%
      fxxinv
    scale <- mean(f)
    serr <- sqrt(diag(cov))
  }
  else if (se == "boot") {
    if ("cluster" %in% names(dots)) {
      bargs <- utils::modifyList(list(x = x, y = y, tau = tau),
                          dots)
      if (length(object$na.action)) {
        cluster <- dots$cluster[-object$na.action]
        bargs <- utils::modifyList(bargs, list(cluster = cluster))
      }
      if (class(bargs$x)[1] == "matrix.csr")
        bargs <- utils::modifyList(bargs, list(control = ctrl))
      B <- do.call(boot.rq, bargs)
    }
    else B <- boot.rq(x, y, tau, coef = coef, ...)
    cov <- cov(B$B)
    serr <- sqrt(diag(cov))
  }
  else if (se == "BLB") {
    n <- length(y)
    b <- ceiling(n^gamma)
    S <- n%/%b
    V <- matrix(sample(1:n, b * S), b, S)
    Z <- matrix(0, NCOL(x), S)
    for (i in 1:S) {
      v <- V[, i]
      B <- boot.rq(x[v, ], y[v], tau, bsmethod = "BLB",
                   blbn = n, ...)
      Z[, i] <- sqrt(diag(cov(B$B)))
    }
    cov <- cov(B$B)
    serr <- apply(Z, 1, mean)
  }
  else if (se == "extreme") {
    tau0 <- tau
    if (tau > 0.5) {
      y <- -y
      tau <- 1 - tau
    }
    if (length(dots$mofn))
      mofn = dots$mofn
    else mofn = floor(n/5)
    if (length(dots$mofn))
      kex = dots$kex
    else kex = 20
    if (length(dots$alpha))
      alpha = dots$alpha
    else alpha = 0.1
    if (length(dots$R))
      R = dots$R
    else R = 200
    m <- (tau * n + kex)/(tau * n)
    taub <- min(tau * n/mofn, tau + (0.5 - tau)/3)
    xbar <- apply(x, 2, mean)
    b0 <- rq.fit(x, y, tau, method = method)$coef
    bm <- rq.fit(x, y, tau = m * tau, method = method)$coef
    An <- (m - 1) * tau * sqrt(n/(tau * (1 - tau)))/c(crossprod(xbar,
                                                                bm - b0))
    bt <- rq.fit(x, y, tau = taub, method = method)$coef
    s <- matrix(sample(1:n, mofn * R, replace = T), mofn,
                R)
    mbe <- (taub * mofn + kex)/(taub * mofn)
    bmbeb <- rq.fit(x, y, tau = mbe * taub, method = method)$coef
    B0 <- boot.rq.pxy(x, y, s, taub, bt, method = method)
    Bm <- boot.rq.pxy(x, y, s, tau = mbe * taub, bmbeb, method = method)
    B <- (mbe - 1) * taub * sqrt(mofn/(taub * (1 - taub))) *
      (B0 - b0)/c((Bm - B0) %*% xbar)
    if (tau0 <= 0.5) {
      bbc <- b0 - apply(B, 2, quantile, 0.5, na.rm = TRUE)/An
      ciL <- b0 - apply(B, 2, quantile, 1 - alpha/2, na.rm = TRUE)/An
      ciU <- b0 - apply(B, 2, quantile, alpha/2, na.rm = TRUE)/An
    }
    else {
      bbc <- -(b0 - apply(B, 2, quantile, 0.5, na.rm = TRUE)/An)
      ciL <- -(b0 - apply(B, 2, quantile, alpha/2, na.rm = TRUE)/An)
      ciU <- -(b0 - apply(B, 2, quantile, 1 - alpha/2,
                          na.rm = TRUE)/An)
    }
    B <- R - sum(is.na(B[, 1]))
    coef <- cbind(b0, bbc, ciL, ciU)
    if (tau0 > 0.5) {
      coef <- -coef
      tau <- tau0
    }
    dimnames(coef) = list(dimnames(x)[[2]], c("coef", "BCcoef",
                                              "ciL", "ciU"))
  }
  else if (se == "conquer") {
    if (length(dots$R))
      R = dots$R
    else R = 200
    Z <- conquer::conquer(x[, -1], y, tau, ci = TRUE, B = R)
    coef <- cbind(Z$coef, Z$perCI)
    cnames <- c("coefficients", "lower bd", "upper bd")
    dimnames(coef) <- list(vnames, cnames)
    resid <- y - x %*% Z$coef
  }
  if (se == "rank") {
    coef <- f$coef
  }
  else if (!(se %in% c("conquer", "extreme"))) {
    coef <- array(coef, c(p, 4))
    dimnames(coef) <- list(vnames, c("Value", "Std. Error",
                                     "t value", "Pr(>|t|)"))
    coef[, 2] <- serr
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- if (rdf > 0)
      2 * (1 - pt(abs(coef[, 3]), rdf))
    else NA
  }

  # important statistics
  SAE = object$SAE
  N = object$N
  MAE = object$MAE
  lambda_ROS = object$lambda_ros
  object <- object[c("call", "terms")]
  if (covariance == TRUE) {
    if (se != "rank")
      object$cov <- cov
    if (se == "iid")
      object$scale <- scale
    if (se %in% c("nid", "ker")) {
      object$Hinv <- fxxinv
      object$J <- crossprod(x)
      object$scale <- scale
    }
    else if (se == "boot") {
      object$B <- B$B
      object$U <- B$U
    }
  }
  object$coefficients <- coef
  object$residuals <- resid
  object$rdf <- rdf
  object$tau <- tau
  object$SAE = SAE
  object$N = N
  object$MAE = MAE
  object$lambda_ros = lambda_ROS
  class(object) = "summary.regL1" # troquei a classe regL1
  object
}


#' Print an regL1 object
#'
#' Print an object generated by regL1
#'
#' @param x Object returned from regL1 representing the fit of the L1 model.
#' @param ... Optional arguments.
#'
#' @returns No return value, called for side effects.
#'
#' @import stats
#' @export
#'
#' @seealso
#' \code{\link{regL1}} for fitting linear L1 models.
# Define a summary method for the "regL1" class
print.regL1 = function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  coef <- stats::coef(x)
  cat("\nCoefficients:\n")
  print(coef, ...)
  rank <- x$rank
  nobs <- length(residuals(x))
  if (is.matrix(coef))
    p <- dim(coef)[1]
  else p <- length(coef)
  rdf <- nobs - p
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
  if (!is.null(attr(x, "na.message")))
    cat(attr(x, "na.message"), "\n")
  invisible(x)
}


#' Change object class from "regL1" to "rq"
#'
#' Changing the object's class from "regL1" to "rq" allows you to use functions from the 'quantreg' package normally.
#'
#' @param object Object from "regL1" class.
#' @returns Object with class "rq".
#'
#' @export
#'
#' @seealso
#' \code{\link{regL1}} for fitting linear L1 models.
#' \code{\link{rq}} for fitting linear L1 models.
class_to_rq = function(object) {
  class(object) = "rq"
  return(object)
}

#' Change object class from "rq" to "regL1"
#'
#' Changing the object's class from "rq" to "regL1" allows you to use functions from the 'diagL1' package normally.
#'
#' @param object Object from "rq" class.
#' @returns Object with class "regL1".
#'
#' @export
#'
#' @seealso
#' \code{\link{regL1}} for fitting linear L1 models.
#' \code{\link{rq}} for fitting linear L1 models.
class_to_regL1 = function(object) {
  class(object) = "regL1"
  return(object)
}
