# Compute curvi-linear trends
#
# Coding history:
#   2012Apr20 DLLorenz Initial dated version
#   2012Sep12 DLLorenz Added style option
#   2013Mar26 DLLorenz Added to package and some tweaks
#

curvi <- function(x, ..., style=c("mw", "se")) {
  fcn <- function(x, m, w) {
    retval <- 0.5 + ((x - m)/w)/(1 + abs(x - m)/w)
    retval <- pmin(1, pmax(0, retval))
    return(retval)
  }
  style <- match.arg(style)
  dots <- list(...)
  if(style == "mw")
    nms <- sapply(dots, function(x) paste("m", x[1], ",w", x[2], sep=''))
  else if(style == "se") {
    nms <- sapply(dots, function(x) paste(x[1], "-", x[2], sep=''))
    ## but the conversion function needs m and w
    dots <- lapply(dots, function(se) return(c(mean(se), diff(se)/2)))
  }
  retval <- lapply(dots, function(tt, dt) fcn(dt, tt[1], tt[2]), dt=x)
  retval <- do.call("cbind", retval)
  colnames(retval) <- nms
  return(retval)
}
