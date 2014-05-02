# new version of move.1 (Simple LOC as opposed to SLR)
#
# Coding history:
#    2006Dec05 DLLorenz Initial version
#    2009Mar19 DLLorenz Dated version
#    2010Dec08 DLLorenz Conversion to R
#    2011Aug05 DLLorenz Added print and plot (diagPlot) files
#    2011Oct25 DLLorenz Update for package
#    2012Aug28 DLLorenz Change from diagPlot to plot
#    2012Nov21 DLLorenz Added Q-Q plot to diagnostics
#    2012Dec21 DLLorenz Change gaussian to normal
#    2013Apr09 DLLorenz Added setGD to plot
#

move.1 <- function(formula, data, subset, na.action, distribution="normal") {
  ## Arguments:
  ##  formula (a formula with 1 left and 1 right entry) the model formula
  ##  data (a data.frame or envirnment) where to find the data in the formula
  ##  subset (any subset expression) subset the data?
  ##  na.action (a function) how to hande missing values
  ##  distribution (character scalar) use bivariate "normal" or
  ##    "lognormal" or "commonlog" distribution
  ##
  ## Process formula as with regular linear regression!
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$distribution <- m$modeluse.log10 <- NULL
  m$drop.unused.levels <- TRUE
  m[[1L]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  if(ncol(m) != 2L)
    stop("Must specify a single response and a single explanatory variable in the formula")
  Y <- m[, 1L]
  yname <- names(m)[1L]
  X <- m[, 2L]
  xname <- names(m)[2L]
  fit <- list(coefficients=double(2L))
  ## Extract missing value info
  fit$na.action <- attr(m, "na.action")
  ## Adjust for transforms
  distribution <- match.arg(distribution, c("normal", "lognormal", "commonlog"))
  if(distribution == "lognormal") {
    Y <- log(Y)
    X <- log(X)
    xname <- paste("log(", xname, ")", sep='')
    yname <- paste("log(", yname, ")", sep='')
  }
  else if(distribution == "commonlog") {
    Y <- log10(Y)
    X <- log10(X)
    xname <- paste("log10(", xname, ")", sep='')
    yname <- paste("log10(", yname, ")", sep='')
  }
  ybar <- mean(Y)
  yvar <- var(Y)
  xbar <- mean(X)
  xvar <- var(X)
  retcor <- cor.test(X, Y)
  fit$R <- retcor$estimate
  fit$p.value <- retcor$p.value
  fit$coefficients[2L] <- sqrt(yvar/xvar)
  if(fit$R < 0)
    fit$coefficients[2L] <-  - fit$coefficients[2L]
  fit$coefficients[1L] <- ybar - fit$coefficients[2L] * xbar
  names(fit$coefficients) <- c("(Intercept)", xname)
  ## Complete steps for the model
  fit$call <- call
  fit$fitted.values <- cbind(1, X) %*% as.matrix(fit$coefficients)
  yresiduals <- Y - fit$fitted.values
  ## Now reverse the sense of the fit to get X-residuals
  junk <- double(2L)
  junk[2L] <- 1/fit$coefficients[2L]
  junk[1L] <- xbar - junk[2L] * ybar
  xfitted <- cbind(1, Y) %*% as.matrix(junk)
  xresiduals <- X - xfitted
  fit$residuals <- cbind(yresiduals, xresiduals)
  colnames(fit$residuals) <- c(yname, xname)
  ## Pack the object with other necessary information
  fit$x <- X
  fit$y <- Y
  fit$xstats <- c(mean=xbar, sd=sqrt(xvar))
  fit$ystats <- c(mean=ybar, sd=sqrt(yvar))
  fit$var.names <- c(yname, xname)
  fit$model <- m
  ## set class to move.1
  oldClass(fit) <- "move.1"
  return(fit)
}

print.move.1 <- function(x, digits=4, ...) {
  ## Arguments:
  ##  x (a move.1 object) the object to print
  ##  digits (integer scalar) the number of significant digits to print
  ##  ... (dots) any other arguments to print
  ##
  cat("Call:\n")
  dput(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, digits=digits, ...)
  cat("\nStatistics of the variables:\nResponse (", x$var.names[1], "):\n",
      sep="")
  print(x$ystats, digits=digits, ...)
  cat("Predictor (", x$var.names[2], "):\n", sep="")
  print(x$xstats, digits=digits, ...)
  cat("Correlation coefficient: ", round(x$R, digits), 
      "\n                p-value: ", round(x$p.value, digits), "\n", sep="")
  if(!is.null(x$na.action)) {
    n.na <- length(x$na.action)
    if(n.na > 1)
      cat("  (", n.na, " observations deleted due to missing values)\n", sep="")
    else
      cat("  (", n.na, " observation deleted due to missing values)\n", sep="")
  }
  if(!is.null(x$cx) && !is.null(x$cy))
    cat(sum(x$cx | x$cy), " observations were left-consored\n", sep="")
  invisible(x)
}

plot.move.1 <- function(x, which="All", set.up=TRUE, span=0.8, ...) {
  ## Arguments:
  ##  x (a move.1 object) the object to plot
  ##  which (character scalar or numeric vector) which plots?
  ##  span (numeric scalar) the span for the x on y, y on x plots
  ##  set.up and ... (dots) not used, required for method functions
  ##
  ## Identify which plots to do:
  ## 1 Q-Q plot
  ## 2 x on y, y on x
  ## 
  ## Set up graphics page
  if(set.up) 
    setGD("MOVE.1")
  ## Set up to do all plots
  doPlot <- c(rep(TRUE, 2L) , FALSE)
  if(is.numeric(which)) {
    if(min(which) > 0) # select which to plot
      doPlot[seq(2L)[-which]] <- FALSE
    else # which not to plot
      doPlot[-which] <- FALSE
  }
  xname <- x$var.names[2L]
  yname <- x$var.names[1L]
  if(doPlot[3L]) {
    ## This is done only by special request
    AA <- scalePlot(x$x, x$y, scale=x$coef[2L], Plot=list(what="points"),
                    xtitle=xname, ytitle=yname)
    refLine(coefficients=x$coef, current=AA)
    cov <- x$R*(x$xstats[2L]*x$ystats[2L])
    cov <- matrix(c(x$xstats[2L]^2, cov, cov, x$ystats[2L]^2), 2)
    el <- cov2Ellipse(cov, c(x$xstats[1L], x$ystats[1L]))
    addXY(el$x, el$y, Plot=list(what="lines", color="blue"), current=AA)
    el <- cov2Ellipse(cov, c(x$xstats[1L], x$ystats[1L]), scale=2)
    addXY(el$x, el$y, Plot=list(what="lines", color="blue"), current=AA) 
  }
  if(doPlot[2L]) {
    ## Set up for 2 graphs
    AA.lo <- setLayout(width=6, height=c(3,3))
    AA.gr <- setGraph(1, AA.lo)
    AA <- xyPlot(x$x, x$y, Plot=list(what="points"), xtitle=xname, ytitle=yname)
    refLine(coefficients=lsfit(x$x, x$y)$coef, Plot=list(color="green"), current=AA)
    addSmooth(x$x, x$y, span=span, Plot=list(color="cyan"), current=AA)
    AA.gr <- setGraph(2, AA.lo)
    AA <- xyPlot(x$y, x$x, Plot=list(what="points"), xtitle=yname, ytitle=xname)
    refLine(coefficients=lsfit(x$y, x$x)$coef, Plot=list(color="green"), current=AA)
    addSmooth(x$y, x$x, span=span, Plot=list(color="cyan"), current=AA)
  }
  if(doPlot[1L]) {
    qqPlot(x$x, x$y, xtitle=xname, ytitle=yname, Line1.1=list(what="none"))
  }
  invisible(x)
}

predict.move.1 <- function(object, newdata, type = c("response", "link")
                           , ...) {
  ## Arguments:
  ##  object (a move.1 object) the model from which predictions are made
  ##  newdata (a data.frame) the new data, if missing use old data
  ##  type (character scalar) what to predict
  ##  ... (dots) not used, required for method function
  ##
  if(missing(newdata))
    newdata <- object$model
  xlab <- attr(attr(object$model, "terms"), "term.labels")
  xindex <- newdata[, xlab] # xindex is simple vector
  ckdist <- FALSE
  if(!is.null(object$call$distribution)) {
    dist <- match.arg(object$call$distribution, c("normal", "lognormal",
                                                  "commonlog"))
    if(dist == "lognormal") {
      ckdist <- TRUE
      xindex <- log(xindex)
    }
    else if(dist == "commonlog") {
      ckdist <- TRUE
      xindex <- log10(xindex)
    }
  }
  out <- cbind(1, xindex) %*% object$coef
  type <- match.arg(type)
  if(type == "response" && ckdist) {
    if(dist == "commonlog")
      out <- 10^out
    else # Must be natural log
      out <- exp(out)
  }
  as.vector(out) # Strip matrix info
}

## The default extraction methods: coef, residuals, fitted work

