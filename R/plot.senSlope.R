plot.senSlope <- function(x, which="All", set.up=TRUE, span=0.8, ...) {
  ## Coding history:
  ##    2013Apr15 DLLorenz Original Coding
  ##
  ## Arguments:
  ##  x (a senSlope object) the object to plot
  ##  which (character scalar or numeric vector) which plots?
  ##  span (numeric scalar) the span for the x on y, y on x plots
  ##  set.up and ... (dots) not used, required for method functions
  ##
  ## Identify which plots to do:
  ## y on x with regression line
  ##
  ## Set up graphics page
  if(set.up) 
    setGD("SenSlope")
  ## Set up to do all plots
  doPlot <- TRUE    
  if(is.numeric(which)) {    
    if(which[1L] == -1) # why is beyond me!
      doPlot <- FALSE    
    if(length(which) > 1 || which != 1)
    warning("Only one diagnostic plot for senSlope")
  }
  xname <- x$var.names[2L]
  yname <- x$var.names[1L]
  if(doPlot[1L]) {
    AA <- xyPlot(x$x, x$y, Plot=list(what="points"), xtitle=xname, ytitle=yname)
    refLine(coefficients=x$coefficients, Plot=list(color="green"), current=AA)
    addSmooth(x$x, x$y, span=span, Plot=list(color="cyan"), current=AA)
  }
  invisible(x)
}

