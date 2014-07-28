# Compute quantiles and corresponding nonparametric confidence limits
# Note this should become a generic function when the USGSqw package is released
#
# Coding history
#    2007Mar21 DLLorenz Original coding
#    2007Mar26 DLLorenz modified for USGS lib 3.2
#    2008Feb28 DLLorenz Bug Fix for na.rm=F, and added type=2.
#    2011Oct25 DLLorenz Update for package
#    2012Aug13 DLLorenz Bug fix for short vectors
#    2014Jun26 DLLorenz Added bound.
#

qtiles.CI <- function(x, probs=0.5, CI=0.90, bound=c("two.sided", "upper", "lower"),
											na.rm=TRUE) {
  ## Arguments:
  ##  x (numeric vector) the data for which to compute the stats
  ##  probs (numeric vector) The probability levels desired
  ##  CI (numeric scalar) the approximate confidence interval
  ##  na.rm (logical scalar) remove missing values
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
