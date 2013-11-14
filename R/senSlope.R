# The Sen slope estimator.
#
# Coding History:
#    2000Oct27 JRSlack  Original coding as kensen
#    2011May02 DLLorenz Conversion to R--as sen slope only
#    2011Oct25 DLLorenz Update for package
#    2013Apr30 DLLorenz Bug fixes
#

senSlope <- function(formula, data, subset, na.action, intercept='Ac',
                     CI=.95) {
  ## Arguments:
  ##  formula (a formula with 1 left and 1 right entry) the model formula
  ##  data (a data.frame or envirnment) where to find the data in the formula
  ##  subset (any subset expression) subset the data?
  ##  na.action (a function) how to hande missing values
  ##  intercept (character scalar) how to compute the intercept term.
  ##
  ## Two options are available for the computation of the intercept:
  ##   'Ac' the default, is the Conover method, which passes through the
  ##   median of x and the median of y
  ##   'A1m' forces the median of the residuals to be 0.
  ##
  ## Define the variance of Kendall's S function, required for conf. int.
  vark <- function(y, x) {
    ties.y <- rle(sort(y))$lengths
    ties.x <- rle(sort(x))$lengths
    n <- length(y)
    t1 <- n * (n - 1) * (2 * n + 5)
    ty2 <- sum(ties.y * (ties.y - 1) * (2 * ties.y + 5))
    tx2 <- sum(ties.x * (ties.x - 1) * (2 * ties.x + 5))
    v <- (t1 - ty2 -tx2)/18
    v1 <- sum(ties.y * (ties.y - 1)) * sum(ties.x * (ties.x - 1)) /
      (2 * n * (n - 1))
    v2 <- sum(ties.y * (ties.y - 1) * (ties.y - 2)) *
      sum(ties.x * (ties.x - 1) * (ties.x - 2)) /
        (9 * n * (n - 1) * (n - 2))
    v <- v + v1 + v2
    return (v)
  }
  ## Process formula as with regular linear regression!
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$intercept <- m$CI <- NULL
  m$drop.unused.levels <- TRUE
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  if(ncol(m) != 2L)
    stop("Must specify a single response and a single explanatory variable in the formula")
  ## Extract missing value info
  na.action <- attr(m, "na.action")
  y <- m[, 1L]
  yname <- names(m)[1L]
  x <- m[,2]
  xname <- names(m)[2L]
  n <- length(y)
  if (n < 3L)
    stop("model data should effectively be longer than 2.")  
  ## Calculate the statistics
  ## A fast, efficient way to compute all possible slopes:
  slopes <- unlist(lapply(seq(along = y), function(i, y, t)
                                    ((y[i] - y[1L:i]) / (t[i] - t[1L:i])), y, x))
  slopes <- sort(slopes) # removes missings (due to ties in x)
  sen.slope <- median(slopes)
  ## Compute the variance of S, accounting only for ties in x
  varS <- vark(seq(along=y), x) # Forces no ties in y
  Calpha <- -qnorm((1-CI)/2) * sqrt(varS)
  M1 <- as.integer((length(slopes) - Calpha)/2)
  M2 <- as.integer((length(slopes) + Calpha)/2) + 1
  sen.CI <- slopes[c(M1, M2)]
  names(sen.CI) <- paste(c("Lower", "Upper"),
                         substring(sprintf("%.2f", CI), 2), sep="") # drop 0
  ## Median of the data values.
  median.y <- median(y)         
  ## Median of the time values.
  median.x <- median(x)
  
  ## A line representing the trend of the data then is given by
  ##
  ##    y(t) = sen.slope*(t-med.time)+med.data
  ##
  ##    This line has the slope of the trend and passes through
  ##       the point (t=med.time, y=med.data)
  ## Compute the coefficients for the line: intercept and slope
  coef <- c(sen.slope*(-median.x)+median.y, sen.slope)
  fits <- coef[1L] + coef[2L]*x
  resids <- y - fits
  if(intercept == 'A1m') {
    coef[1L] <- coef[1L] + median(resids)
    resids <- y - coef[1L] - coef[2L]*x
  }
  ## Return the statistics.
  retval <- list(call=call, coefficients=coef, slope.CI=sen.CI, residuals=resids, 
                 fitted.values=fits, na.action=na.action, x=x, y=y, var.names=c(yname, xname),
                 model=m)
  oldClass(retval) <- "senSlope"
  return(retval)
}
