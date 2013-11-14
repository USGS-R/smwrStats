# function to confirm a peak in noisy data
#
# Coding history:
#    2007Sep14 DLLorenz Start of confirm
#    2007Oct10 DLLorenz Bug fix to confirm, GUI=F
#    2011May25 DLLorenz Begin Conversion to R and rename
#    2012Aug11 DLLorenz Integer fixes
#    2013Jun17 DLLorenz  Final tweaks for release
#

confirm.seasonalPeak <- function(x, all=FALSE, plot.only=FALSE, ...) {
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
      X11(title="Confirm")
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

