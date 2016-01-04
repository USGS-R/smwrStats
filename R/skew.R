#' Skewness
#' 
#' Computes the skewness statistic.
#' 
#' @param x any numeric vector.
#' @param na.rm logical; if TRUE, then remove missing value before computation.
#' @param method the method to use for computing the skew, must be "fisher" or
#' "moments."
#' @return a single value representing the skewness of the data in \code{x}.
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' @keywords univariate
#' @examples
#' 
#' skew(c(1.0, 1.2, 1.5, 1.9, 2.5))
#' 
#' @export skew
skew <- function(x, na.rm=TRUE, method="fisher") {
	# Coding history:
	#    Unknown   DLLorenz Original Coding 
	#    2011Aug24 DLLorenz Conversion to R
	#    2013Feb28 DLLorenz Tweak to handling method
	#    2014Dec29 DLLorenz Conversion to roxygen header
  ##
  method <- match.arg(method, c("fisher", "moments"))
  if(na.rm)
    x <- x[!is.na(x)]
  n <- length(x)
  mn <- mean(x)
  dif.x <- x - mn
  m2 <- sum(dif.x^2)/n
  m3 <- sum(dif.x^3)/n
  b1 <- (m3/(m2^(3/2)))
  if(method == "moments")
    return(b1)
  if(n < 3)
    g1 <- NA
  else
    g1 <- (sqrt(n * (n - 1)) * b1)/(n - 2)
  return(g1)
}
