#' Seasonal Peak Timing
#' 
#' Compute the timing of the seasonal peak value. The timing of the seasonal
#' peak is needed for the seasonalWave model (Vecchia and others, 2008).
#' 
#' The timing of the peak of the data is computed by identifying the largest
#' value produced by smoothing that data with \code{supsmu}. The remaining data
#' in the attributes are used by using the seasonalPeak method of confirm.
#' 
#' @param x a vector of decimal time representing dates and times. Missing
#' values are permitted and are removed before analysis.
#' @param y a vector of the data for which the peak is needed. Missing values
#' are permitted and are removed before analysis.
#' @return An object of class seasonalPeak. The unconfirmed object is a single
#' value that represents the estimate of the timing of the peak and five
#' additional attributes.\cr
#' 
#' Data: a list of the \code{x} and \code{y} values where x is the fractional
#' part of the original decimal time data. Missing values have been removed.\cr
#' Smooth: a list of the \code{x} and \code{y} smoothed values. \cr Points: a
#' list of 361 evenly spaced xout and yout values. \cr Extra: pointers to all
#' the peaks in \code{Points}. \cr Confirmed: logical indicating that the
#' object has not been confirmed.
#' @note The generic functions print and confirm have methods for object of
#' class seasonalPeak.
#' @seealso \code{\link{seasonalWave}}, \code{\link{confirm.seasonalPeak}},
#' \code{\link{supsmu}}, \code{\link{print.seasonalPeak}}
#' @references Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling
#' variability and trends in pesticide concentrations in streams: Journal of
#' the American Water Resources Association, v. 44, no. 5, p. 1308-1324
#' @keywords manip
#' @examples
#' 
#' library(smwrData)
#' data(QW05078470)
#' with(QW05078470, seasonalPeak(dectime(DATES), P00665))
#' ## Should be:
#' # Default value: 0.499 
#' # Alternate values: 0.497 
#' 
#' @export seasonalPeak
seasonalPeak <- function(x, y) {
	# Coding history:
	#    2007????? SVecchia Original Coding in seasonal wave regression
	#    2007Aug22 DLLorenz Separate function and added to USGS library
	#    2007Aug29 DLLorenz Added hlife attributes to return value
	#    2007Sep14 DLLorenz Bug fix to confirm and print
	#    2007Oct10 DLLorenz Bug fix to confirm, GUI=F
	#    2011May25 DLLorenz Begin Conversion to R and rename
	#    2012Aug11 DLLorenz Integer fixes
	#    2013Apr03 DLLorenz Final tweaks for release
	#    2014Dec29 DLLorenz Conversion to roxygen header
  ##
  ## Remove missing values
  Nas <- is.na(x) | is.na(y)
  x <- x[!Nas]
  y <- y[!Nas]
  ## Find the peak
  x <- x - floor(x)
  tmpsmu <- supsmu(x, y, periodic = TRUE)
  xtmp <- tmpsmu$x
  ytmp <- tmpsmu$y
  ntmp <- length(ytmp)
  cmax <- xtmp[order(ytmp)[ntmp]]
  ymax <- max(ytmp)
  ## find other peaks
  xout <- seq(0,1,1/360)
  yout <- approx(xtmp, ytmp, xout=xout, rule=3)$y
  ## ydiff = 0 => peak (decreasing values) or trough (increasing values)
  ydiff <- diff(c(yout[359L:361L], yout, yout[seq(3L)]), lag=6)
  ydecrease <- eventNum(ydiff <= 0, reset=TRUE)
  ## find the first value of the decrease
  xindx <- tapply(seq(361L), ydecrease, min)[-1L] # drop 0
  ## do not double count the first point
  if(xindx[1] == 1 && ydecrease[361L] > 0)
    xindx <- xindx[-1L]
  
  retval <- cmax
  attr(retval, "Data") <- list(x=x, y=y)
  attr(retval, "Smooth") <- tmpsmu
  attr(retval, "Points") <- list(xout=xout, yout=yout)
  attr(retval, "Extra") <- xindx
  attr(retval, "Confirmed") <- FALSE
  oldClass(retval) <- "seasonalPeak"
  return(retval)
}
