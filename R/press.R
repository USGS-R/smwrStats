# Compute PRESS statistic
#
# Coding history:
#    2005Jul14 DLLorenz Initial dated verion
#    2011Aug09 DLLorenz Conversion to R
#    2012Nov06 DLLorenz Bug fix to account for influence returning 0 for NAs
#    2012Nov06          This version.
#

press <- function(model) {
  ## Argument:
  ##  model (a regression model) the object
  ##   currently object must be lm, or the output from lsfit
  ##
  if(inherits(model, "lm")) {
    h <- influence(model)$hat
    r <- resid(model)
    retval <- sum((r/(1 - h))^2, na.rm=TRUE)
  }
  else if(!is.null(model$residuals) && !is.null(model$qr)) {
    ## Assumed lsfit
    h <- ls.diag(model)$hat
    ## lsfit returns the same length and corresponding values even if
    ## there are missing values, so we can compute the correct value!
    retval <- sum((model$residuals/(1 - h))^2, na.rm=TRUE)
  }
  else
    stop("Input must be a 'lm' or 'lsfit' object")
  return(retval)
}
