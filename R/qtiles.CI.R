# Compute quantiles and corresponding nonparametric confidence limits
# Note this should become a generic function when the USGSqw package is released
#
# Coding history
#    2007Mar21 DLLorenz Original coding
#    2007Mar26 DLLorenz modified for USGS lib 3.2
#    2008Feb28 DLLorenz Bug Fix for na.rm=F, and added type=2.
#    2011Oct25 DLLorenz Update for package
#    2012Aug13 DLLorenz Bug fix for short vectors
#    2012Aug13          This version.
#

qtiles.CI <- function(x, probs=0.5, CI=0.90, na.rm=TRUE) {
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
  N <- length(x)
  x <- sort(x)
  qtiles <- quantile(x, probs=probs, type=2)
  lci <- qbinom((1 - CI) / 2, N, probs)
  ## Temporarily replace 0 with N + 1 to extract correctly
  lci <- ifelse(lci == 0, N + 1, lci)
  x.lci <- x[lci]
  ## Fix lci
  lci <- ifelse(lci == N + 1, 0, lci)
  uci <- qbinom(0.5 + CI/2, N, probs) + 1	
  x.uci <- x[uci]
  ci.actual <- pbinom(uci-1, N, probs) - pbinom(pmax(0, lci - 1), N, probs)
  ## Fix 0s for very small probabilities, and NAs
  ci.actual <- ifelse(ci.actual < 0.1e-6, 1, ci.actual)
  ci.actual <- ifelse(is.na(x.lci + x.uci), NA, ci.actual)
  retval <- cbind(estimate=qtiles, lcl=x.lci, ucl=x.uci, ci=ci.actual)
  rownames(retval) <- names(qtiles)
  return(retval)
}
