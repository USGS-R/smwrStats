#' Confirm Seasonal Peak
#' 
#' This function processes the output of function seasonalPeak and requires the
#' user to accept or change the timing of the peak, if there is a secondary
#' peak and the necessary characteristics of the peaks.
#' 
#' 
#' @param x a seasonalPeak object.
#' @param all a logical value indicating whether to accept the peak without
#' interactive user input or to force the user to process the peak. The default
#' value is \code{FALSE}, forcing the user to process the peak. Can also be set
#' to either 1 or 2, indicating the number of peaks.
#' @param plot.only a logical value indicating that only a plot is desired. If
#' \code{TRUE}, then \code{x} is returned invisibly and unchanged.
#' @param \dots not used, required for compatibility with other methods.
#' @return An object of class seasonalPeak. The confirmed object is a single
#' value that represents the estimate of the timing of the peak and four or
#' five aadditional ttributes.\cr
#' 
#' NumPeaks: the number of seasonal peaks; either 1 or 2.\cr Models: candidate
#' loading models. The number indicates the number of months of constituent
#' loading.\cr hlife: candidate half-life values. The muner indicates the
#' half-life in terms of months. \cr Confirmed: logical indicating that the
#' object has not been confirmed.\cr
#' 
#' If \code{NumPeaks} is 2, then an additional attirbute \code{Second.peak}
#' that is a data frame of candidate parameters for the second peak is
#' included. See \code{\link{seasonalWave}} for details.
#' @seealso \code{\link{seasonalWave}}, \code{\link{seasonalPeak}}
#' @references Vecchia, A.V., Martin, J.D., and Gilliom, R.J., 2008, Modeling
#' variability and trends in pesticide concentrations in streams: Journal of
#' the American Water Resources Association, v. 44, no. 5, p. 1308-1324.
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' library(smwrData)
#' data(QW05078470)
#' # Simply click on the identified peak, and enter 1 for a single peak.
#' confirm(with(QW05078470, seasonalPeak(dectime(DATES), P00665)))
#' }
#' 
#' @export
#' @method confirm seasonalPeak
confirm.seasonalPeak <- function(x, all=FALSE, plot.only=FALSE, ...) {
	# Coding history:
	#    2007Sep14 DLLorenz Start of confirm
	#    2007Oct10 DLLorenz Bug fix to confirm, GUI=F
	#    2011May25 DLLorenz Begin Conversion to R and rename
	#    2012Aug11 DLLorenz Integer fixes
	#    2013Jun17 DLLorenz  Final tweaks for release
	#    2014Dec22 DLLorenz	 Roxygen header
	#
  ## This function forces the user to select a peak identified in the
  ## seasonalPeak() function or add a new peak.
  ##
  if(attr(x, "Confirmed")) {
    warning("x already confirmed")
    return(x)
  }
  ## Plot the data and smooth:
  Sel <- as.integer(all) # 0 is interactive
  if(Sel == 0L) {
    if(!plot.only)
      setGD("Confirm")
    Data <- attr(x, "Data")
    Smooth <- attr(x, "Smooth")
    plot(Data, xlim=c(0,1), xlab='x', ylab='y')
    lines(Smooth, col = 3)
    ## plot a guideline for the second peak
    yrange <- range(Smooth$y)
    yguide <- yrange[1]  + diff(yrange)/3
    abline(h=yguide, col=3, lty=2)
    ## plot the peak and alternates
    abline(v=as.double(x), col=2)
    Points <- attr(x, "Points")
    xindx <- attr(x, "Extra")
    points(Points$xout[xindx], Points$yout[xindx], col=2, pch=3)
    ## Select peak
    title(main='Select Main Peak')
    if(plot.only)
      return(invisible(x))
    retval <- locator(1)$x
    ## Try to determine range of models
    Sel <- menu(c("YES", "NO"), title='Single peak?')
  } else
    retval <- as.vector(x)
  attr(retval, "NumPeaks") <- Sel
  if(Sel == 1L)
    attr(retval, "Models") <- seq(9L)
  else {
    attr(retval, "Models") <- c(2L, 3L)
    sp <- data.frame(la = c(3,  3,  4,  4,  5,  5,  6,  6,  7,  7),
                     lo = c(1,  2,  1,  2,  1,  2,  1,  2,  1,  2),
                     w  = c(1,.75,  1,.75,  1,.75,  1,.75,  1,.75))
    attr(retval, "Second.peak") <- sp
  }
  attr(retval, "hlife") <- seq(4)
    attr(retval, "Confirmed") <- TRUE
  oldClass(retval) <- "seasonalPeak"
  return(retval)
}

