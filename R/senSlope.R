#' Compute the Sen Slope
#' 
#' Computes the Sen slope with confidence interval and an intercept for paired
#' data.
#' 
#' The argument \code{intercept} may be either "Ac" or "A1m." If it is "Ac,"
#' then the intercept is computed from the median of \code{y} and \code{x},
#' also known as the Conover method. If it is "A1m," then the intercept is
#' chosen so that the median of the residuals is zero.
#' 
#' @param formula a model formula with exactly one explanatory variable.
#' @param data the data.
#' @param subset any descriptiopn to subset \code{data}.
#' @param na.action the function to handle missing values.
#' @param intercept a character string indicating the method to compute the
#' intercept. See \bold{Details}.
#' @param CI the desired confidence interval for the slope.
#' @return An object of class "senSlope" with these components: \item{call}{the
#' matched call.} \item{coefficients}{the intercept ans Sen slope.}
#' \item{slope.CI}{the lower and upper confidence limits of the Sen slope.}
#' \item{residuals}{the residuals of the regression.} \item{fitted.values}{the
#' fitted values.} \item{na.action}{information about any missing values.}
#' \item{x}{the explanatory variable.} \item{y}{the response variable.}
#' \item{var.names}{the response and explanatory variable names.}
#' \item{model}{the model frame.}
#' @seealso \code{\link{kensen.test}}, \code{\link{serial.test}}
#' @references Dietz, E.J., 1989, Teaching regression in a nonparametric
#' statistics course: The American Statstician, v. 43, p. 35--40\cr
#' 
#' Helsel, D.R., and Hirsch, R.M., 2002, Statistical methods in water
#' resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' 
#' Sen, P.K., 1968, Estimates of the regression coefficient based on Kendall's
#' tau: Journal of the American Statistical Association, v. 63 p. 1379--1389\cr
#' @keywords regression nonparametric robust
#' @examples
#' 
#' \dontrun{
#' library(smwrData)
#' data(SaddlePeaks)
#' senSlope(Flow ~ Year, data=SaddlePeaks)
#' }
#' 
#' @export senSlope
senSlope <- function(formula, data, subset, na.action, intercept='Ac',
                     CI=.95) {
	# Coding History:
	#    2000Oct27 JRSlack  Original coding as kensen
	#    2011May02 DLLorenz Conversion to R--as sen slope only
	#    2011Oct25 DLLorenz Update for package
	#    2013Apr30 DLLorenz Bug fixes
	#    2014Dec29 DLLorenz Conversion to roxygen header
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
