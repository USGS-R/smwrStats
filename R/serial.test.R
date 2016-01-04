#' Test for Serial Correlation
#' 
#' Performs either of two nonparametric tests (Wilcoxon Rank-Sum test or runs
#' test) for serial correlation. Both tests compute Z, which is the product of
#' sequential values of \code{x} after removing any trend and normalizing so that the 
#' median of Z values is 0.
#' 
#' If the method is "wilcoxon," then the Wilcoxon Rank-Sum test (Wilcoxon,
#' 1945) is performed to compare the distribution of positive Z values to the
#' distribution of negative Z values. Excessive numbers of positive Z
#' values suggests positive serial correlation and excessive negative Z
#' values suggests negative serial correlation.
#' 
#' If the method is "runs," then the runs test (Wald and Wolfowitz, 1940) is
#' performed on sequences of positive and negative values of Z. If there is a
#' larger than expected number of runs, then that suggest negative serial correlation
#' and if fewer than the expected number of runs, then positive serial correlation.
#' 
#' @param x the numeric vector of observations. Missing values (NAs) are not
#' allowed.
#' @param method the string "wilcoxon" or "runs," depending on which test
#' should be used.  Only the first character is necessary. See \bold{Details}.
#' @return An object of class htest containing the following components:
#' \item{statistic}{ the score for either test. } \item{p.value}{ the two-sided
#' p-value of the statistic. } \item{Z}{ the raw data of the analysis.  }
#' \item{alternative}{ a string "two.sided" indicating the hypothesis. }
#' \item{method}{ a description of the method.  } \item{data.name}{ a string
#' containing the actual name of the input data.  }
#' @note The null hypothesis is that the data are uncorrelated.
#' @references Dufour, J.M., 1981, Rank test for serial dependence: Journal of
#' the Time Series Analysis, v. 2, no. 3, p. 117--128.
#' 
#' Wald, A., and Wolfowitz, J., 1940, On a test whether two samples are from
#' the same population: Annals of Mathematical Statistics, v. 11, p.
#' 147--162.
#' 
#' Wilcoxon, F., 1945, Individual comparisons by ranking methods: Biometrics,
#' v. 1, p. 80--83.
#' @keywords nonparametric htest
#' @importFrom randtests runs.test
#' @export serial.test
serial.test <- function(x, method = "wilcoxon") {
	# Coding history:
	#    2007Apr20 DLLorenz Original Coding
	#    2007May08 DLLorenz Added to USGS library
	#    2011Oct25 DLLorenz Update for package
	#    2013Mar25 DLLorenz fixed case in metthod output
	#    2014Dec29 DLLorenz Conversion to roxygen header
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
    # Use the wilcox.test
    exact <- length(Z) < 50 && !any(duplicated(Z))
    tmp <- wilcox.test(Z[Z >= 0], -Z[Z < 0], exact=exact)
    S <- tmp$statistic
    p.value <- tmp$p.value
  } # end of wilcoxon
  else {
    method <- "Runs" # Fix case for output
    ## Test the number of times sequential values of x have the same sign
    tmp <- runs.test(as.numeric(Z >= 0), threshold=0.5) # coded as 1,0
    S <- c(Runs=tmp$runs)
    p.value <- tmp$p.value
  }
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
