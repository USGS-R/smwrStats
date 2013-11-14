# compute the area under the ROC curve
#
# Coding history:
#    2008Apr03 DLLorenz Original Coding in plotLogistic()
#    2009Feb11 DLLorenz Moved to this file
#    2011Aug22 DLLorenz Conversion to R
#    2011Oct25 DLLorenz Update for package
#    2012Aug28 DLLorenz Change from diagPlot to plot
#    2013Apr09 DLLorenz Added setGD to plot
#    

roc <- function(object) {
  ## Argument:
  ##  object (a glm model object) the logistic regression model for ROC
  ##
  fits <- fitted(object)
  ## Force to integer, just in case there is some off-1 values
  y <- as.integer(round(object$y,3))
  ## Sort y and fits by fits
  fits.order <- order(fits)
  y.sort <- y[fits.order]
  fits.sort <- fits[fits.order]
  sens <- 1-cumsum(y.sort)/sum(y.sort)
  spec <- cumsum(1-y.sort)/sum(1-y.sort)
  roc.mat <- do.call("rbind", by(cbind(fits=fits.sort, spec=spec, sens=sens),
                               spec, function(x) x[x$sens == max(x$sens),]))
  ## Compute c
  c.tbl <- table(outer(fits[y==1], fits[y==0], function(x,y) sign(x-y)))
  if(length(c.tbl) == 3) # ties in predicted values
    c.val <- (c.tbl[3] + .5*c.tbl[2])/sum(c.tbl)
  else # no ties
    c.val <- c.tbl[2]/sum(c.tbl)
  retval <- list(c.val=c.val, table=roc.mat)
  oldClass(retval) <- "roc"
  return(retval)
}

print.roc <- function(x, digits=3, ...) {
  ## Arguments:
  ##  x (roc object) the object to print
  ##  digits (integer scalar) how many digits to use when prining
  ##   ... (dots) not used, required for method function
  ##
  cat("\nArea under the ROC curve: ", round(x$c.val, digits), "\n\n", sep='')
  invisible()
  return(x)
}

plot.roc <- function(x, which="All", set.up=TRUE, ...) {
  ## Arguments:
  ##  x (roc object) the object to plot
  ##   ... (dots) unused, required for method function
  ## 
  ## Set up graphics page
  if(set.up) 
    setGD("ROC")
  ## Set up to do all plots
  doPlot <- TRUE    
  if(is.numeric(which)) {    
    if(which[1L] == -1) # why is beyond me!
      doPlot <- FALSE    
    if(length(which) > 1 || which != 1)
    warning("Only one diagnostic plot for senSlope")
  }
  if(doPlot[1L]) {
    sens <- x$table$sens
    spec<- 1-x$table$spec
    fit <- x$table$fit
    ## Need to start at origin
    xyPlot(c(1,spec,0), c(1,sens,0),
           ytitle='True positive rate (Sensitivity)',
           xtitle='False positive rate (1-Specificity)',
           margin=c(NA,NA, 2.2, NA))
    refLine(coefficients=c(0,1))
    ## Select labels along the length of the curve and plot the predicted value
    curve.dist <- cumsum(diff(c(0, 1 - spec))^2 + diff(c(0, 1 - sens))^2)
    for(i in seq(.05, .95, by=0.05) * max(curve.dist)) {
      xdist <- abs(curve.dist - i)
      min.xdist <-  min(xdist)
      pick <- which(xdist == min.xdist)[1]
      xpick <- spec[pick]
      ypick <- sens[pick]
      if(fit[pick] > 0.008 && fit[pick] < 0.993) { # label it
        if(xpick > .04 && ypick < .96) {
          addAnnotation(xpick-0.01, ypick-0.01,
                        format(round(fit[pick], 2)), angle=-45, justification='right')
        }
        else { # Put on the left side of the line
          addAnnotation(xpick+0.01, ypick-0.03,
                        format(round(fit[pick], 2)), angle=-45, justification='left')
        }
      }
    }
    addTitle(Main=paste('ROC Analysis\nArea under curve =',
               round(x$c.val,3), sep=' '))
  }
  invisible(x)
}
