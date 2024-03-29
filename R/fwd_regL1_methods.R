
#' Summarizing Fit of Forward Search in Linear L1 Regression
#'
#' Returns a summary list for a forward search in linear L1 regression fit.
#'
#'
#' @param object Object returned from forwardSearch_regL1.
#' @param steps the number of forward steps to show.
#' @param ... Optional arguments.
#'
#' @returns No return value, called for side effects.
#' @export
#'
#' @references Atkinson, A.C. and Riani, M. (2000). \emph{Robust Diagnostic Regression Analysis}. New York: Springer.
#'
#'
#' @method summary forwardSearch_regL1
#' @rdname summary.forwardSearch_regL1
#'
#' @seealso
#' \code{\link{forwardSearch_regL1}} for apply forward search in linear L1 regression model.
#'
#' @examples
#' \donttest{
#' # applies the forward search approach to robust analysis in a linear L1 model
#' mod = forwardSearch_regL1(Concentration ~ Age, data = bile)
#' summary(mod)
#' }
#'

summary.forwardSearch_regL1 = function(object, steps = "auto", ...)
{
  n = dim(object$Residuals)[1]
  p = dim(object$Coefficients)[2]

  if (steps == "auto")
    steps <- max(5, as.integer(0.05 * n))
  else
    if (steps == "all")
      steps <- n-p #-1

  if (steps + p > n)
    stop("Too many steps")
  steps <- steps - 1

  if (any(!is.na(object$Unit[n-p-steps, 2])))
    Unit <- object$Unit[(n-p-steps+1):(n-p+1), ]
  else
    Unit <- t(object$Unit[(n-p-steps+1):(n-p+1), 1])

  int <- intersect(1:n, as.vector(object$Unit[(n - p - steps):(n - p),  ]))
  Residuals <- object$Residuals[int, (n - p - steps + 1):(n - p + 1)]

  MinDelRes <- object$MinDelRes[(n - p - steps - 1):(n - p - 1)]
  names(MinDelRes) <- paste("m=", as.character((n - steps - 1):(n - 1)), sep = "")

  Coefficients <- object$Coefficients[(n - p - steps + 1):(n - p + 1),  ]


  m_start = n-dim(object$tStatistics)[1]
  print(dim(object$tStatistics))
  tStatistics <- object$tStatistics[(n - m_start - steps):(n - m_start), ]


  MAE <- object$MAE[(n - p - steps):(n - p)]
  names(MAE) <- paste("m=", as.character((n - steps):(n)), sep = "")

  #R2 <- object$R2[(n - p - steps):(n - p)]
  #names(R2) <- paste("m=", as.character((n - steps):(n)), sep = "")

  last_summary <- list(call = object$call, Unit = Unit,
               Residuals = Residuals,
               MinDelRes = MinDelRes, Coefficients = Coefficients,
               tStatistics = tStatistics, MAE = MAE)
               #, R2 = R2)

  class(last_summary) <- "summary.forwardSearch_regL1"
  last_summary
}

#' Print an forwardSearch_regL1 object
#'
#' Print an object generated by forwardSearch_regL1
#'
#'
#' @param x Object returned from forwardSearch_regL1.
#' @param ... Optional arguments.
#'
#' @returns No return value, called for side effects.
#'
#' @export
#'
#' @seealso
#' \code{\link{forwardSearch_regL1}} for apply forward search in linear L1 regression model.
#'
#'
#' @method print forwardSearch_regL1
#' @rdname print.forwardSearch_regL1
#'
#' @examples
#' \donttest{
#' # applies the forward search approach to robust analysis in a linear L1 model
#' mod = forwardSearch_regL1(Concentration ~ Age, data = bile)
#' mod # or print(mod)
#' }
#'

print.forwardSearch_regL1 = function(x, ...)
{
  cat("Call:\n")
  print(x$call)

  n = dim(x$Residuals)[1]
  p = dim(x$Coefficients)[2]

  cat("\nLast 5 units included in the forward search:\n")
  units <- x$Unit[(n-p-3):(n-p+1), 1, drop = FALSE]
  dimnames(units)[[2]] <- ""
  print(t(units))

  cat("\nEstimated coefficients for 50% of the data and for the full model:\n")

  print(x$Coefficients[c(ceiling(n/2) - p + 1, n - p + 1), ])
  cat("\n")

  invisible(x)
}

#' Print Forward Search in Linear L1 Model Summary Object
#'
#' Print summary of forward search in linear L1 model object.
#'
#'
#' @param x This is an object of class "\code{summary}" produced by a call to \code{summary.regL1()}.
#' @param digits Significant digits reported in the printed table.
#' @param ... Optional arguments.
#' @returns No return value, called for side effects.
#'
#' @export
#'
#' @seealso
#' \code{\link{forwardSearch_regL1}} for apply forward search in linear L1 regression model.
#'
#'
#' @import greekLetters
#' @method print summary.forwardSearch_regL1
#' @rdname print.summary.forwardSearch_regL1
#'
#' @examples
#' \donttest{
#' # applies the forward search approach to robust analysis in a linear L1 model
#' mod = forwardSearch_regL1(Concentration ~ Age, data = bile)
#' summary(mod)
#' }

print.summary.forwardSearch_regL1 = function(x, digits = 4, ...)
{
  cat("Call:\n")
  print(x$call)

  steps <- dim(x$Unit)[2]
  cat(paste("\nUnits included in the last", steps, "steps:\n"))
  if(dim(x$Unit)[1] == 1)
    dimnames(x$Unit)[1] <- "Unit"
  print(x$Unit)

  cat("\nResiduals:\n")
  print(x$Residuals, digits = digits)

  cat("\nMinimum Deletion Residuals:\n")
  print(x$MinDelRes, digits = digits)

  cat("\nEstimated Coefficients:\n")
  print(x$Coefficients, digits = digits)

  cat("\nt Value:\n")
  print(x$tStatistics, digits = digits)

  cat("\nMAE - Mean Absolute Error:\n")
  print(x$MAE, digits = digits)

  #cat("\nR^2:\n")
  #print(x$R2, digits = digits)

  invisible(x)
}


#' Forward Search in Linear L1 Models
#'
#' This function plots the results of a forward search in linear L1 models.
#'
#' @param x a "forwardSearch_regL1" object.
#' @param type.plot select which plots to draw, by default all. Each graph is addressed by an integer:
#' 1) scaled residuals
#' 2) minimum deletion residuals
#' 3) coefficients
#' 4) statistics
#' 5) MAE (Mean Absolute Error) values
# #' 6) R2 values
#'
#' @param squared logical, if TRUE plots squared residuals.
#' @param scaled logical, if TRUE plots scaled coefficient estimates.
#' @param ylim a two component vector for the min and max of the y axis.
#' @param xlim a two component vector for the min and max of the x axis.
#' @param th.Res numerical, a threshold for labelling the residuals.
#' @param th.Lev numerical, a threshold for labelling the leverages.
#' @param sig.Tst numerical, a value (on the scale of the t statistics) used to draw the confidence interval on the plot of the t statistics.
#' @param labels.in.plot logical, if TRUE units are labelled in the plots when required.
#'
#' @param ... additional arguments.
#'
#' @returns No return value, just plots the results of a forward search in linear L1 models.
#'
#' @importFrom graphics par lines text abline
#'
#' @export
#'
#' @seealso
#' \code{\link{forwardSearch_regL1}} for apply forward search in linear L1 regression model.
#'
#'
#' @examples
#' \donttest{
#'
#' # applies the forward search approach to robust analysis in a linear L1 model
#' mod = forwardSearch_regL1(Concentration ~ Age, data = bile)
#' plot(mod, 1)
#' }

"plot.forwardSearch_regL1" = function(x, type.plot = 1:5, squared = FALSE,
                        scaled = FALSE, ylim = NULL, xlim = NULL,
                        th.Res = 2, th.Lev = 0.25, sig.Tst = 2.58,
                        labels.in.plot = TRUE, ...)
  {
    if(!is.numeric(type.plot) || max(type.plot) > 5 || min(type.plot) < 1)
      stop("type.plot must be an integer vector containing only 1 to 5")

    oldpar <- graphics::par()
    #  on.exit(par(oldpar))
    if (length(type.plot)>1)
    { graphics::par(ask = TRUE)
      on.exit(par(ask=oldpar$ask)) }

    n <- dim(x$Residuals)[1]
    p <- dim(x$Coefficients)[2]

    ## Plot Residuals ##

    if (is.element(1, type.plot))
    {
      x1 <- x$Residuals

      if (squared)
      { x1 <- x1^2
      ylab <- "Squared Scaled Residuals" }
      else
      { ylab <- "Scaled Residuals" }

      if (is.null(ylim))
        y.lim <- c(ifelse(squared, 0, min(x1, na.rm = TRUE)),
                   max(x1, na.rm = TRUE))
      else
        y.lim <- ylim

      if (is.null(xlim))
        x.lim <- c(p - 1, n + 1)
      else
        x.lim <- xlim

      plot(p:n, x1[1, ],
           type = "n",
           ylim = y.lim,
           xlim = x.lim,
           xlab = "Subset Size",
           ylab = ylab)

      for (i in 1:nrow(x1))
        graphics::lines(p:n, x1[i, ], lty = i)

      ### Write numbers for selected units which
      ### have residuals bigger than a certain threshold

      ma <- apply(abs(x1), 1, max, na.rm = TRUE)
      nam <- which(ma > th.Res)

      if (lnam <- length(nam))
      {
        graphics::text(rep(p, lnam), x1[nam, 1], paste(nam, " ", sep = ""), adj = 1)
        graphics::text(rep(n, lnam), x1[nam, ncol(x1)], paste(" ", nam, sep = ""), adj = 0)
      }
    }


    ## Plot minimum deletion residuals ##

    if (is.element(2, type.plot))
    {

      x1 <- x$MinDelRes

      if (nrow(x1) > 3)
      { a <- p + 4
      x1 <- x1[4:nrow(x1),] }
      else
        a <- p + 1

      if (is.null(ylim))
        y.lim <- c(0, max(x1, na.rm = TRUE))
      else
        y.lim <- ylim

      if (is.null(xlim))
        x.lim <- c(a-1, n)
      else
        x.lim <- xlim

      plot(a:(n-1), as.numeric(x1),
           type = "l",
           xlab = "Subset Size",
           ylab = "Minimum Deletion Residuals",
           ylim = y.lim,
           xlim = x.lim
           #col ='blue'
           )
    }

    ## Plot Coefficients    ##

    if (is.element(3, type.plot))
    {
      x1 <- na.omit(x$Coefficients)
      n.dropped <- n - p + 1 - dim(x1)[1]

      if (scaled)
      { x1 <- apply(x1, 2, function(x) sign(mean(x))*x/mean(x))
      ylab <- "Scaled Coefficient Estimates" }
      else
      { ylab <- "Coefficient Estimates" }

      if (is.null(ylim))
        y.lim <- range(x1, na.rm = TRUE)
      else
        y.lim <- ylim

      if (is.null(xlim))
        x.lim <- c(p, ifelse(labels.in.plot, p+1.1*(n-p), n))
      else
        x.lim <- xlim

      plot(0, 0,
           type = "n",
           xlab = "Subset Size",
           ylab = ylab,
           ylim = y.lim,
           xlim = x.lim)

      for (i in 1:ncol(x1))
        graphics::lines((p + n.dropped):n, x1[, i], lty = i, lwd = 1)

      if (labels.in.plot)
        # graphics::text(rep(n,p), x1[dim(x1)[1],], labels = paste(" C", 0:(p-1), sep = ""), cex = 1.1, adj = 0)
        # Mathematical annotation
        for (i in 1:p)
        { graphics::text(rep(n,p), x1[dim(x1)[1],i],
               substitute(hat(beta)[b], list(b=i-1)), pos=4) }
    }

    ## Plot t statistics ##

    if (is.element(4, type.plot))
    {
      x1 <- x$tStatistics

      if (is.null(ylim))
      { y.lim <- as.integer(.667 * dim(x1)[1]):(dim(x1)[1])
      y.lim <- range(x1[y.lim, , drop = FALSE], na.rm = TRUE)
      #y.lim <- max(abs(range(y.lim)))
      y.lim <- max(abs(range(x1)))
      y.lim <- c(-1*y.lim, y.lim)
      }
      else
        y.lim <- ylim

      if (is.null(xlim))
        x.lim <- c(p, ifelse(labels.in.plot, p+1.1*(n-p), n))
      else
        x.lim <- xlim

      #         x1[x1 <= y.lim[1] | x1 >= y.lim [2]] <- NA

      plot(0, 0,
           type = "n",
           xlab = "Subset Size",
           ylab = "t Value",
           ylim = y.lim,
           xlim = x.lim)

      #Significance lines
      if (is.numeric(sig.Tst))
      { graphics::abline(h = sig.Tst, col = "lightgrey")
        graphics::abline(h = -sig.Tst, col = "lightgrey") }

      m_start = n-dim(x1)[1]
      for (i in 1:ncol(x1))
        graphics::lines((m_start +1):n, x1[ ,i], lty = i, lwd = 1)


      if (labels.in.plot)
        #      graphics::text(rep(n,p), x1[n-p,], labels = paste(" C", 0:(p-1), sep = ""), cex = 1.1, adj = 0)
        # Mathematical annotation
        for (i in 1:p)
        { graphics::text(rep(n,p), x1[dim(x1)[1],i],
               substitute(hat(beta)[b], list(b=i-1)), pos=4) }

    }


    ## Plot MAE - Mean Absolute Error   ##

    if (is.element(5, type.plot))
    {
      x1 <- x$MAE

      if (is.null(ylim))
        y.lim <- c(0, max(x1[2:nrow(x1)], na.rm = TRUE))
      else
        y.lim <- ylim

      if (is.null(xlim))
        x.lim <- c(p - 1, n + 1)
      else
        x.lim <- xlim

      plot(p:n, x1[,1],
           type = "l",
           xlab = "Subset Size",
           ylab = "MAE - Mean Absolute Error",
           ylim = y.lim,
           xlim = x.lim)
    }

    ## Plot R^2     ##

    # if (is.element(6, type.plot))
    # {
    #   x1 <- x$R2
    #
    #   if (is.null(ylim))
    #     y.lim <- c(0, max(x1[2:nrow(x1)], na.rm = TRUE))
    #   else
    #     y.lim <- ylim
    #
    #   if (is.null(xlim))
    #     x.lim <- c(p - 1, n + 1)
    #   else
    #     x.lim <- xlim
    #
    #   plot(p:n, x1[,1],
    #        type="l",
    #        xlab = "Subset Size",
    #        ylab = expression(R^2),
    #        ylim = y.lim,
    #        xlim = x.lim)
    #
    # }

    invisible(x)
  }
