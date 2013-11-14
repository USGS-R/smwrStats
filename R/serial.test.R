# a nonparametric test for serial correlation
#
# Coding history:
#    2007Apr20 DLLorenz Original Coding
#    2007May08 DLLorenz Added to USGS library
#    2011Oct25 DLLorenz Update for package
#    2013Mar25 DLLorenz fixed case in metthod output
#

serial.test <- function(x, method = "wilcoxon") {
  ## Arguments:
  ##  x (numeric vector) the data to test
  ##  method (character scalar) the name of the test
  ##
  data.name <- deparse(substitute(x))
  if(any(is.na(x)))
    stop("Missing values not permitted")
  ## x is assumed to be an annual series with no missing data
  method <- match.arg(method, c("wilcoxon", "runs"))
  ## Step 1 detrend and normalize the data
  N <- length(x)
  t <- seq(N)
  slopes <- unlist(lapply(seq(along = x), function(i, y, t)
                          ((y[i] - y[1:i])/(t[i] - t[1:i])), x, t))
  slope <- median(slopes, na.rm = TRUE)
  x <- x - slope * t
  x <- x - median(x)
  ## Step 2 compute Z1 and the test
  Z <- x * shiftData(x)
  ## Remove the first (NA) value
  Z <- Z[-1]
  if(method == "wilcoxon") {
    method <- "Wilcoxon" # Fix case for output
    R <- rank(abs(Z))
    S <- sum((Z >= 0) * R)
    m <- sum(Z >= 0)
    ## Number >= 0
    n <- N - m - 1
    ## Continuity adjustment
    if(max(m, n) < 50 && min(m, n) > 0)
      p.value <- pwilcox(S - 1, m, n)
    else {
      ## Use large-sample approximation
      varS <- sum(R^2)/4
      ES <- sum(R)*m/(m + n)
      p.value <- pnorm((S - ES)/sqrt(varS))
    }
  } # end of wilcoxon
  else {
    method <- "Runs" # Fix case for output
    ## Test the number of times sequential values of x have the same sign
    S <- sum(Z >= 0)
    p.value <- pbinom(N - S, N, 0.5)
  }
  ## Correct p-value for two-sided test and build the components
  p.value <- 1 - abs(p.value - 0.5) * 2
  names(S) <- "S"
  zero <- 0
  names(zero) <- "lag-1 serial dependence"
  ## Return the htest
  retval <- list(data.name = data.name,
                 method = paste(method, "test for serial dependence", sep = " "),
                 statistic = S, p.value = p.value, alternative = "two.sided",
                 null.value = zero, Z = Z)
  oldClass(retval) <- "htest"
  return(retval)
}
