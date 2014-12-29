#' Basis for Piecewise Linear Trends
#' 
#' Generate a basis matrix for piecewise linear modeling of trends.
#' 
#' 
#' @param x a vector of dates/times, assumed to be in dectime format. Missing
#' values are permitted and result in corresponding missing values in the
#' output.
#' @param breaks a vector of breakpoints in the linear trends.
#' @param boundary.breaks a logical vector of length 2, indicating whether the
#' breaks include the range of the data.  If the first is TRUE, then the first
#' break denotes the beginning of a trend.  If the last is TRUE, then the last
#' break denotes the end of a trend. Otherwise, the trends begin at the
#' \code{floor} of the first value of \code{x} and end at the \code{ceiling} of
#' the last value of \code{x}.
#' @param steps a vector indicating any step trends.
#' @return A matrix of dimension length of x by number of linear trend and step
#' trends.  The breaks are included as an attribute.
#' @note Each trend is 0 prior to the break, and then increases at the rate of
#' 1 per unit and maintain that maximum value after the break. The regression
#' coefficient then reprents the trend as a rate. A step trend is 0 before the
#' step and 1 after it.
#' @seealso \code{\link{floor}}, \code{\link{ceiling}}, \code{\link{curvi}}
#' @keywords model
#' @examples
#' 
#' # model two piecewise linear trends from 2000 to 2004, with a break at 2001 and 2003
#' trends(2000 + seq(0,20)/5, breaks=c(2001, 2003))
#' 
#' @export trends
trends <- function(x, breaks, boundary.breaks=c(FALSE, FALSE), steps) {
	# Coding history:
	#    2004Nov23 DLLorenz Original
	#    2011Aug31 DLLorenz Conversion to R
	#    2011Oct25 DLLorenz Update for package
	#    2013Aug15 DLLorenz Added class
	#    2014Dec29 DLLorenz Conversion to roxygen header
  ##
  boundary.breaks <- rep(boundary.breaks, length.out=2) # just in case
  if(!boundary.breaks[1])
    breaks <- c(floor(min(x)), breaks)
  if(!boundary.breaks[2])
    breaks <- c(breaks, ceiling(max(x)))
  breaks <- sort(breaks)
  nCol <- length(breaks) - 1
  if(!missing(steps))
    nStep <- length(steps)
  else
    nStep <- 0
  retval <- matrix(0, nrow=length(x), ncol=nCol+nStep)
  for(i in 1:nCol) {
    y <- x[x > breaks[i]]
    retval[x > breaks[i], i] <- ifelse(y <= breaks[i+1], y - breaks[i], breaks[i+1] - breaks[i])
  }
  breakNames <- paste(breaks[1:nCol], breaks[-1], sep="-")
  if(nStep > 0) {
    for(i in 1:nStep)
      retval[x > steps[i], nCol + i] <- 1
    stepNames <- paste("step:", steps, sep="")
  }
  else
    stepnames <- NULL
  dimnames(retval) <- list(NULL, c(breakNames, stepnames))
  attr(retval, "breaks") <- breaks
  class(retval) <- "trends"
  return(retval)
}
