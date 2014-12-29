#' Quantiles with Confidence Limits
#' 
#' Compute sample quantiles and confidence limts for specified probabilities.
#' 
#' 
#' @param x numeric vector to compute the sample quantiles.
#' @param probs numeric vector of desired probabilities with values between 0
#' and 1.
#' @param CI the minimum desired confidence interval for each level specifed in
#' probs.
#' @param bound a character string indicating the desired bounds, "two.sided"
#' means the two-sided interval, "upper" means the upper bound of the interval,
#' and "lower" means the lower bound of the interval. Only a single character
#' is needed. The lower confidence limit is \code{-Inf} when \code{bound} is
#' "upper" and the upper confidence limit is \code{Inf} when \code{bound} is
#' "lower."
#' @param na.rm logical; if TRUE, then missing values are removed before
#' computation.
#' @return A matrix of sample quantiles, the lower confidence limit, the upper
#' confidence limit, and the probability represented by the confidence interval
#' corresponding to the probs levels in the sorted x data.
#' @seealso \code{\link{quantile}}
#' @references Helsel, D.R. and Hirsch, R.M., 2002, Statistical methods in
#' water resources: U.S. Geological Survey Techniques of Water-Resources
#' Investigations, book 4, chap. A3, 522 p.\cr
#' @keywords univar
#' @examples
#' 
#' ## Generate a random sample
#' set.seed(222)
#' XX.rn <- rexp(32)
#' qtiles.CI(XX.rn, probs=c(.25, .5, .75))
#' 
#' @export qtiles.CI
qtiles.CI <- function(x, probs=0.5, CI=0.90, bound=c("two.sided", "upper", "lower"),
											na.rm=TRUE) {
	# Coding history
	#    2007Mar21 DLLorenz Original coding
	#    2007Mar26 DLLorenz modified for USGS lib 3.2
	#    2008Feb28 DLLorenz Bug Fix for na.rm=F, and added type=2.
	#    2011Oct25 DLLorenz Update for package
	#    2012Aug13 DLLorenz Bug fix for short vectors
	#    2014Jun26 DLLorenz Added bound.
	#    2014Dec29 DLLorenz Convert to roxygen headers
  ##
  if(na.rm)
    x <- x[!is.na(x)]
  if(any(is.na(x))) { # return all NA if any NA
    retval <- cbind(estimate=NA*probs, lcl=NA, ucl=NA, ci=NA)
    rownames(retval) <- format(probs)
    return(retval)
  }
  bound <- match.arg(bound)
  if(bound == "two.sided") {
  	alpha <- (1 - CI) / 2
  } else
  	alpha <- (1 - CI)
  N <- length(x)
  x <- sort(x)
  qtiles <- quantile(x, probs=probs, type=2)
  lci <- qbinom(alpha, N, probs)
  ## Temporarily replace 0 with N + 1 to extract correctly
  lci <- ifelse(lci == 0, N + 1, lci)
  x.lci <- x[lci]
  ## Fix lci
  lci <- ifelse(lci == N + 1, 0, lci)
  uci <- qbinom(1 - alpha, N, probs) + 1	
  x.uci <- x[uci]
  if(bound == "two.sided") {
  	ci.actual <- pbinom(uci-1, N, probs) - pbinom(pmax(0, lci - 1), N, probs)
  } else if(bound == "upper") {
  	ci.actual <- pbinom(uci-1, N, probs)
  	x.lci <- -Inf
  } else { # must be lower
  	ci.actual <- 1 - pbinom(pmax(0, lci - 1), N, probs)
  	x.uci <- Inf
  }
  ## Fix 0s for very small probabilities, and NAs
  ci.actual <- ifelse(ci.actual < 0.1e-6, 1, ci.actual)
  ci.actual <- ifelse(is.na(x.lci + x.uci), NA, ci.actual)
  retval <- cbind(estimate=qtiles, lcl=x.lci, ucl=x.uci, ci=ci.actual)
  rownames(retval) <- names(qtiles)
  return(retval)
}
