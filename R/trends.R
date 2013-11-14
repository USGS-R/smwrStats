# function to compute terms for piecewise modeling of trends
#
# Coding history:
#    2004Nov23 DLLorenz Original
#    2011Aug31 DLLorenz Conversion to R
#    2011Oct25 DLLorenz Update for package
#    2013Aug15 DLLorenz Added class
#

trends <- function(x, breaks, boundary.breaks=c(FALSE, FALSE), steps) {
  ## Arguments:
  ##  x (numeric vector) the data on which to make the trends
  ##  breaks (numeric vector) the breakpoints for the piecewise segments
  ##  boundary.breaks (logical of length 2) extend breaks to boundary
  ##  steps (numeric vector) the chnage points for step change
  ## 
  ## Creates separate trends of x at breaks, x is assumed to be dectime format
  ## Each trend is 0 prior to the break, increase at 1 per unit and maintain
  ## that maximum value after the break
  ## If boundary.breaks=T, then breaks completely define the trends, otherwise
  ## the breaks are interior breaks (the default)
  ## steps indicates a step trend 0 before, 1 after
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
