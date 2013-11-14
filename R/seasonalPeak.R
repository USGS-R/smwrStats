# function to find a peak in noisy data
# Now that this is in stats, we can use a more sophisticated method for
# finding peaks and second peaks. Do we need to create a S4 method for
# lcens y?
#
# Coding history:
#    2007????? SVecchia Original Coding in seasonal wave regression
#    2007Aug22 DLLorenz Separate function and added to USGS library
#    2007Aug29 DLLorenz Added hlife attributes to return value
#    2007Sep14 DLLorenz Bug fix to confirm and print
#    2007Oct10 DLLorenz Bug fix to confirm, GUI=F
#    2011May25 DLLorenz Begin Conversion to R and rename
#    2012Aug11 DLLorenz Integer fixes
#    2013Apr03 DLLorenz Final tweaks for release
#

seasonalPeak <- function(x, y) {
  ## Arguments:
  ##  x (numeric vector) is time, expressed in decimal format
  ##  y (numeric vector) is noisy data, for example chemical
  ## concentrations in water
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
