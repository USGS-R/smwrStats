#' Trend Test
#' 
#' Computes the regional Kendall trend test with Sen slope estimator.
#' 
#' 
#' @param series a numeric matrix with rows representing the annual observations and
#'columns representing the sites. Missing values are permitted.
#' @param correct logical, if \code{TRUE}, then apply the correction for
#'cross correlation among sites. if \code{FALSE}, then do not apply the correction.
#' @return An object of class "htest" also inhereting class "seaken" containing 
#' the following components:
#' \item{method}{ a description of the method. } 
#' \item{statistic}{ the value of Kendall's tau. } 
#' \item{p.value}{ the p-value. See \bold{Note}. }
#' \item{p.value.raw}{ the p-value computed without correction for cross
#' correlation. See \bold{Note}. } 
#' \item{p.value.corrected}{ the p-value computed with correction for cross
#' correlation. See \bold{Note}. }
#' \item{estimate}{ a named vector containing the Sen estimate of the slope in
#' units per year, the median value of the data, and the median value of time.} 
#' \item{data.name}{ a string containing the actual name of the input series
#' with the number of years and sites } 
#' \item{alternative}{ a character string describing alternative to the test 
#' ("two.sided"). }
#' \item{null.value}{ the value for the hypothesized slope (0).}
#' \item{nyears}{ the number of years.}
#' \item{nseasons}{ the number of sites}
#' \item{series}{ the data that was analyzed.}
#' \item{seasonnames}{ the names of the sites}
#' @note The value of \code{p.value} is \code{p.value.raw} if there are fewer
#' than 10 years of data and is \code{p.value.corrected} otherwise.
#' 
#' The regional Kendall is described in Douglas and others (2000) and Helsel and others (2006).
#'The adjustment for spatial correlation used in \code{regken} is based on equation 13 in 
#'Douglas and others (2000) and uses the Spearman correlation to account for monotonic correlations.
#' @seealso \code{\link{seaken}}
#' @references 
#'Douglas, E.M., Vogel, R.M., and Kroll, C.N., 2000, Trends in floods and low flows in
#'the United States: impact of spatial correlation: Journal of Hydrology, v. 240, p. 90--105.
#'
#'Helsel, D.R., Mueller, D.K., and Slack, J.R., 2006, Computer program for the Kendall 
#'family of trend tests: U.S. Geological Survey Scientific Investigations Report 2005--5275, 4 p.
#'
#' @keywords htest
#' @examples
#' # Need example?
#' 
#' @useDynLib smwrStats seakenf
#' @export regken
regken <- function(series, correct=nrow(series) > 9) {
	#
	# Do error checking before calling seaken.
	Dname <- deparse(substitute(series))
	ns <- as.integer(ncol(series))
	n <- as.integer(length(series))
	nseas <- colnames(series)
	x <- as.single(t(series)) # transform so that each site is grouped by year
	x[is.na(x)] <- as.single(-99999.0)   # Convert any NAs to -99999.0
	
	# Call the Fortran seaken DLL
	# x is the sreies grouped by site, not year
	# n is the length of the timeseries
	# ns is the negative number of sites
	# results is the statistics of the trend
	results <- .Fortran("seakenf", x , n , ns , as.single(vector("numeric", 9)))[[4]]
	## Return the statistics.
	tau <- results[1]
	names(tau) <- "tau"
	zero <- 0
	names(zero) <- "slope"
	est <- c(results[4], results[6],  results[5])
	names(est) <- c("slope", "median.data", "median.time")
	## Compute the p.value adjusted for cross correlation
	x.cor <- cor(series,  method="spearman", use="pair")
  rhoxx <- 2*mean(x.cor[lower.tri(x.cor)])
  VarS <- results[8L]*(1 + rhoxx)
  S <- results[7L] - sign(results[7L]) # continuity not applied to return value
  results[3L] <- 2*(1 - pnorm(abs(S)/sqrt(VarS)))
	## use appropriate p.value
	if(correct) {
		p.value <- results[3L]
		method <- "Regional Kendall trend test with correlation correction"
	} else {
		p.value <- results[2L]
		method <- "Regional Kendall trend test"
	}
	x <- miss2na(as.numeric(x), -99999.0)
	z <- list(method = method, data.name =
							paste(Dname, " (", n/ns, " years and ", ns, " sites)", sep=""),
						nyears = n/ns, nseasons = ns, series=x, seasonames=nseas,
						statistic=tau, p.value=p.value,
						p.value.raw = results[2], p.value.corrected = results[3],
						estimate=est, alternative = "two.sided", null.value = zero)
	oldClass(z) <- c("htest", "seaken")
	return(z)
}

