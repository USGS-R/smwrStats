# compute rmse for various objects
# includes the generic and method functions
#
# Coding history:
#   2009Oct02 DLLorenz Original Coding
#   2012May24 DLLorenz Conversion to R
#   2012Sep29          This version.
#

rmse <- function(x, ...)
  UseMethod("rmse")

rmse.default <- function(x, y, ...) {
  ## Arguments:
  ##  x, y, paired sample and duplicate
  ##
 sqrt(mean((x-y)^2)/2) # from Analytical Chemistry
}

rmse.lm <- function(x, ...) {
  ## Arguments:
  ##  x, an object that inherit class lm, such as aov.
  ##
  rdf <- x$df.resid
  if(!is.null(x$weights))
    w <- x$weights
  else
    w <- 1.0
  sqrt(sum(w * x$residuals^2)/rdf)
}

rpd <- function(x, y) {
  ## Arguments:
  ##  x, y, paired sample and duplicate
  ##
 (x - y)/(x + y) * 50 # from standard methods
}
