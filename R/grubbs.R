#' The Grubbs Distribution
#' 
#' The one-sided critical value and the attained signicance level for the
#' Grubbs test for an outlier in a sample from Normal distribution.
#' 
#' 
#' @aliases Grubbs pgrubbs qgrubbs
#' @param G the maximum or minimum scaled difference from the mean.
#' @param N the number of values in the sample.
#' @param alpha the significance level.
#' @return The function \code{pgrubbs} returns the attained p-value, not the
#' cumulative distribution, for the maximum scaled distance from the mean. the
#' funciton \code{qgrubbs} returns the one-sided critical value for the maximum
#' scaled difference from the mean.
#' @seealso \code{\link{grubbs.test}}
#' @references Grubbs, F., 1969, Procedures for Detecting Outlying Observations
#' in Samples, Technometrics, v. 11, no. 1, pp. 1-21.
#' @keywords distribution
#' @examples
#' 
#' # The difference is due to rounding errors
#' pgrubbs(c(.9, .95, .99), 32)
#' qgrubbs(c(5.905348, 5.483234, 5.159097), 32)
#' 
#' @rdname Grubbs
#' @export
qgrubbs <- function(alpha, N) {
	# Coding history:
	#    2008Apr24 DLLorenz Original
	#    2012May11 DLLorenz Conversion to R
	#    2012May21 DLLorenz Roxygen headers
  ##
  ## return the one-sided critical values for alpha
  ## to get two-sided divied alpha by 2
  if(N < 3)
    return(as.double(NA))
  G <- qt(alpha/N, N-2)
  return((N-1)/sqrt(N)*sqrt(G^2/(N-2+G^2)))
}

#' @rdname Grubbs
#' @export
pgrubbs <- function(G, N) {
  ## Arguments:
  ##  G (numeric vector) the maximum or minimum scaled difference from the mean
  ##  N the number of values in the sample
  ##
  ## return the one-sided attained p-value for the observed statistic
  ## G = (Ymax - Ymean)/s or
  ## G = (Ymean - Ymin)/s
  TT <- sqrt((N-2)/((N-1)^2/(N*G^2) - 1))
  pt(-TT, N-2)*N
}
