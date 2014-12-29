#' Test for a Trend
#' 
#' Tests for a temporal trend using Kendall's tau and computes the Sen slope
#' estimate of the trend.
#' 
#' 
#' @param y the data collected over time. Missing values (NAs) are allowed and
#' removed before computations.
#' @param t the time corresponding to each observations in \code{y}. Missing
#' values are allowed only where \code{y} is missing. These should be expressed
#' as Julian or decimal time and must be strictly increasing.
#' @param n.min the minimum number of observations for adjusting the p-value
#' for serial correlation.  Used when \code{t} are uniformly spaced to adjust
#' the p-value for serial correlation. Any value larger than the number of
#' observations in \code{y} or \code{Inf} will suppress the adjustment.
#' @return An object of class "htest" containing the following components:
#' \item{method}{ a description of the method.  } \item{statistic}{ the value
#' of Kendall's tau. } \item{p.value}{ the p-value. } \item{estimate}{ a named
#' vector containing the Sen estimate of the slope in units per year, the
#' median value of the data, the median value of time, the number of
#' observations, and if the serial correction is applied, the effective number
#' of observations (n*). } \item{data.name}{ a string containing the actual
#' name of the input series. } \item{coef}{ a vector of an estimate of the
#' intercept and slope. } \item{alternative}{ a character string describing
#' alternative to the test ("two.sided"). } \item{null.value}{ the value for
#' the hypothesized slope (0).  }
#' @note A straight line of the form \cr \preformatted{ ytrend = sen.slope * (
#' t - median.time ) + median.data} may be used as a trend line for graphically
#' portraying or detrending the data. It goes through the point\cr
#' \preformatted{ (t,y) = ( median.time , median.data )} with slope sen.slope.
#' \cr
#' 
#' Many tied values can cause misleading results.\cr
#' 
#' The p-values for uniformly spaced data (\code{t} values unit value like
#' years) are adjusted for lag-1 autoregressive serial correlation according to
#' the method described by Yue and Wang (2004) that adjusts for trend. In
#' keeping with the logic of \code{seaken}, the p-value adjustment is never
#' performed for fewer than 10 observations. The user can suppress the
#' adjustment by setting the value of \code{n.min} to \code{Inf}.
#' @seealso \code{\link{dectime}}, \code{\link{seaken}}
#' @references Conover, W.J., 1980, Practical nonparametric statistics (2d
#' ed.): New York, Wiley, 512 p.\cr
#' 
#' Helsel, D.R., and Hirsch, R.M., 2002, Statistical methods in water
#' resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' 
#' Hirsch, R.M., Alexander, R.B. , and Smith, R.A., 1991, Selection of methods
#' for the detection and estimation of trends in water quality: Water Resources
#' Research, v. 27 p. 803--813.\cr
#' 
#' Kendall, M.G., 1938, A new measure of rank correlation: Biometrika v. 30, p.
#' 81--89.\cr
#' 
#' Kendall, M.G., 1976, Rank correlation methods (4th ed.): London, Griffin,
#' 202 p.\cr
#' 
#' Sen, P.K., 1968, Estimates of regression coefficient based on Kendall's tau:
#' Journal of the American Statisical Association, v. 63, p. 1379--1389.\cr
#' 
#' Yue, S. and Wang. C., 2004, The Mann-Kendall test modified by effective
#' sample size to detect trend in serially correlated hydrological series:
#' Water Resources Management v. 18, p. 201-218.
#' @keywords htest
#' @examples
#' 
#' \dontrun{
#' library(smwrData)
#' data(SaddlePeaks)
#' with(SaddlePeaks, kensen.test(Flow, Year))
#' }
#' 
#' @export kensen.test
kensen.test <- function(y, t, n.min=10) {
	# Coding History:
	#    2000Oct27 JRSlack  Original coding.
	#    2001Feb01 JRSlack  Additional input checking.
	#    2003Apr04 JRSlack  Change variable x to y plus editorial cleanup.
	#    2005Mar24 DLLorenz Fixed removal of NAs
	#    2005May04 DLLorenz made class htest
	#    2005Jul14 DLLorenz Bug fix
	#    2006Apr10 DLLorenz Modifed to print slope and medians
	#    2007Jul12 DLLorenz Removes 1:1 match between values and time
	#    2007Aug27 DLLorenz Added S and VarS to return value
	#    2007Sep05 DLLorenz Added protection against all ties in data
	#    2008Sep17 DLLorenz Bugfix to computation on n
	#    2011May02 DLLorenz Conversion to R (major changes needed)
	#    2011Oct26 DLLorenz Renamed to indicate that a hypo test is done
	#    2013Mar26 DLLorenz Added serial correction 
	#    2014Dec22 DLLorenz Roxygen headers
	#
  ## Set names
  data.name <- paste(deparse(substitute(y)), "and", deparse(substitute(t)))
  ## Remove missing values and error checking.
  ## Make sure every data value has a time value.
  sel <- !is.na(y)
  y <- y[sel]
  t <- t[sel]
  if (any(is.na(t)))
    stop("Each data value must have a time value.")
  ## Make sure time is strictly increasing
  dift <- diff(t)
  if (any(dift <= 0))
    stop("Time vector must be strictly increasing.")
  ## Make sure we have enough data.
  n <- length(y)
  if (n < 3)
    stop("y and t should effectively be longer than 2.")  
  ## Calculate the statistics
  if(diff(range(y)) == 0) { # Protect against error message from cor.test
    retval <- list(statistic=0, parameter=NULL, p.value=1,
      estimate=c(tau=0), null.value=c(tau=0), 
      alternative="two.sided", method="", data.name="")
    class(retval) <- "htest"
    n.min <- Inf # Suppress adjustment
  }
  else
   retval <- cor.test(y, t, method='kendall', continuity=TRUE, exact=FALSE)
  retval$method <- "Kendall's tau with the Sen slope estimator"
  ## a fast, efficient way to compute slopes:
  slopes <- unlist(lapply(seq(along = y), function(i, y, t)
                                    ((y[i] - y[1:i]) / (t[i] - t[1:i])), y, t))
  ## The Sen slope estimator.
  sen.slope <- median(slopes,na.rm=TRUE)
  ## Median of the data values.
  median.data <- median(y)         
  ## Median of the time values.
  median.time <- median(t)
  ## A line representing the trend of the data then is given by
  ##
  ##    y(t) = sen.slope*(t-med.time)+med.data
  ##
  ##    This line has the slope of the trend and passes through
  ##       the point (t=med.time, y=med.data)
  ## Compute the coefficients for the line
  coef <- c(sen.slope*(-median.time)+median.data, sen.slope) # intercept and slope
  est <- c(sen.slope, median.data, median.time, n)
  names(est) <- c("slope", "median.data", "median.time", "n")
  ## Check if ajdustment can/should be made
  if(n.min < n && n > 9 && max(abs(dift - 1)) < 0.02) { # A small wiggle room
    Resids <- (y - coef[2L] * t)
    cor1 <- acf(Resids, lag.max=1L, plot=FALSE)$acf[-1L, 1L, 1L] # drop 0 lag
    ## Use eqn (7) Matalas & Langbein
    adjS <- 1 + 2 *(cor1^(n+1)-n*cor1^2+(n-1)*cor1)/(n*(cor1-1)^2)
    nstar <- n/adjS
    retval$p.value <- (1 - pnorm(abs(retval$statistic/sqrt(adjS))))*2
    retval$method <- "Kendall's tau (adjusted p-value) with the Sen slope estimator"
    est <- c(est, "n*"=nstar)
  }
  ## Return the statistics.
  retval$statistic <- retval$estimate # Replace T/z with tau
  retval$estimate <- est 
  retval$coef <- coef
  retval$data.name <- data.name
  retval$y <- y
  retval$t <- t
  class(retval) <- c("htest", "kensen")
  return(retval)
}

