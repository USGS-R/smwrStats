#' Curvi-Linear Trends
#' 
#' Generates a matrix for curvi-linear modeling of trends.
#' 
#' Curvi-linear trends are described in Appendix 1 (equations 1--9 through
#' 1--11) in Vecchia (2005). The original form of the trend was described in
#' terms of the midpoint of the trend and the halfwidth. Specifying the start
#' and end times of the trend are added as an easy-to-use option.
#' 
#' @param x a vector of dates/times, assumed to be in dectime format. Missing
#' values are permitted and result in corresponding missing values in the
#' output.
#' @param \dots paired vectors specifying the midpoint of the trend and the
#' halfwidth, or the start and end times.
#' @param style either "mw" for midpoint time and halfwidth, or "se" for start
#' and end times.
#' @return A matrix with one column for each of the trends sprecified in
#' \dots{}.
#' @note Curvi-linear trends provide a pleasing visual trend with a gradual
#' transition between trends in contrast to linear trends with sharp changes at
#' the endpoints.\cr Each trend is 0 prior to the start, and then increases to
#' a maximum of 1 and maintain that maximum value after the end. The overall
#' change is described by the regresison coefficient, rather than as a rate
#' described by \code{\link{trends}}, for example.
#' @seealso \code{\link{trends}},
#' @references Vecchia, A.V., 2005, Water-quality trend analysis and sampling
#' desgin for streams in the Red River fo the North basin, Minnesota, Norht
#' Dakota, and South Dakota, 1970--2001: U.S. Geolgical Survey Scientific
#' Investigations Report 2005--5224, 54 p. Available on line at
#' \url{http://pubs.usgs.gov/sir/2005/5224/}.
#' @keywords model
#' @examples
#' 
#' # Model with a single curvi-linear trend from 2001 through 2003
#' # First using midpoint and half width (default) and then start and end.
#' curvi(2000 + seq(0,20)/5, c(2002, 1))
#' curvi(2000 + seq(0,20)/5, c(2001, 2003), style="se")
#' 
#' @export curvi
curvi <- function(x, ..., style=c("mw", "se")) {
	# Coding history:
	#   2012Apr20 DLLorenz Initial dated version
	#   2012Sep12 DLLorenz Added style option
	#   2013Mar26 DLLorenz Added to package and some tweaks
	#   2014Dec22 DLLorenz Roxygen headers
	#
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
