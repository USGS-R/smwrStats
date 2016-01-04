#' Trend Test
#' 
#' Computes the seasonal Kendall trend test with Sen slope estimator.
#' 
#' 
#' @param series a regularly spaced numeric vector to test for trend. Missing
#' values are permitted.
#' @param nseas the number of seasons per year. Must not exceed 52. Can also
#'be a character vector of the names of the seasons. The length of the 
#'character vector determines the number of seasons.
#' @return An object of class "htest" also inhereting class "seaken" containing 
#' the following components:
#' \item{method}{ a description of the method.  } 
#' \item{statistic}{ the value of Kendall's tau. } 
#' \item{p.value}{ the p-value. See \bold{Note}. }
#' \item{p.value.raw}{ the p-value computed without correction for serial
#' correlation. See \bold{Note}. } 
#' \item{p.value.corrected}{ the p-value computed with correction for serial
#' correlation. See \bold{Note}. }
#' \item{estimate}{ a named vector containing the Sen estimate of the slope in
#' units per year, the median value of the data, and the median value of time.} 
#' \item{data.name}{ a string containing the actual name of the input series
#' with the number of years and seasons. } 
#' \item{alternative}{ a character string describing alternative to the test 
#' ("two.sided"). }
#' \item{null.value}{ the value for the hypothesized slope (0).}
#' \item{nyears}{ the number of years.}
#' \item{nseasons}{ the number of seasons.}
#' \item{series}{ the data that was analyzed.}
#' \item{seasonnames}{ the names of the seasons.}
#' @note The value of \code{p.value} is \code{p.value.raw} if there are fewer
#' than 10 years of data and is \code{p.value.corrected} otherwise.
#' @seealso \code{\link{kensen.test}}, \code{\link{regularSeries}}
#' @references Hirsch, R.M., Alexander, R.B., and Smith, R.A., 1991, Selection
#' of methods for the detection and estimation of trends in water quality:
#' Water Resources Research, v. 27, p. 803--813.\cr
#' 
#' Hirsch, R.M., Slack, J.R., and Smith, R.A., 1982, Techniques of trend
#' analysis for monthly water quality data: Water Resources Research, v. 18, p.
#' 107--121.\cr
#' 
#' Hirsch, R.M., and Slack, J.R., 1984, A nonparametric trend test for seasonal
#' data with serial dependence: Water Resources Research, v. 20, p.
#' 727--732.\cr
#' 
#' Kendall, M.G., 1938, A new measure of rank correlation: Biometrika v. 30, p.
#' 81--89.\cr
#' 
#' Kendall, M.G., 1976, Rank correlation methods (4th ed.): London, Griffin,
#' 202 p.\cr
#' 
#' Sen, P.K., 1968, Estimates of regression coefficient based on Kendall's tau:
#' Journal of the American Statisical Association, v. 63, p. 1379--1389.\cr
#' @keywords htest
#' @examples
#' 
#' \dontrun{
#' library(smwrData)
#' library(smwrBase)
#' data(KlamathTP)
#' RegTP <- with(KlamathTP, regularSeries(TP_ss, sample_dt))
#' # The warning generated is expected and acceptable for these data
#' seaken(RegTP$Value, 12)
#' # Manaus river data is in package boot
#' library(boot)
#' data(manaus)
#' manaus.sk <- seaken(manaus, 12)
#' print(manaus.sk)
#' # Note for these data the large difference between the raw and corrected p-values.
#' #  p-value (raw) is << 0.001
#' manaus.sk$p.value.raw
#' #  p-value (with correlation correction) is = 0.10
#' manaus.sk$p.value.corrected
#' #  Hence, it may be concluded that these particular data show substantial serial correlation
#' #  as seen with see with acf(manaus).
#' }
#' 
#' @useDynLib smwrStats seakenf
#' @export seaken
seaken <- function(series, nseas=12) {
	# Coding history:
	#    2000Oct18 JRSlack  Initial coding.
	#    2000Dec22 JRSlack  Simplify the returned list and maintain mode single.
	#    2002Feb20 JRSlack  Modified for S-PLUS 6 conventions.
	#    2002Dec19 JRSlack  Remove upper limits on seasons and years.
	#    2005Jul14 DLLorenz Date fix
	#    2006Apr10 DLLorenz Modified to ouput htest and print slope and medians
	#    2006Apr11 DLLorenz Added years seasons to data.name.
	#    2006May26 DLLorenz Bug fix in retval
	#    2012Feb09 DLLorenz Conversion to R
	#    2012Aug06 DLLorenz Added series to return value and added "seaken" class
	#    2012Dec17 DLLorenz Bug fix in call to set up axis
	#    2014Dec29 DLLorenz Convert to roxygen header
	#
	# Do error checking before calling seaken.
	Dname <- deparse(substitute(series))
	x <- as.single(series)
	n <- as.integer(length(series))
	if(class(nseas) == "character") {
		ns <- as.integer(length(nseas))
	} else {
		ns <- as.integer(nseas)
		nseas <- as.character(seq(nseas))
	}
	if (n/ns < 2) stop ("seaken requires at least 2 years of data.")
	# Pad any trailing partial year to a full year. The original Fortran
	#    code truncated the series using n<-floor(n/ns)*ns
	if (n%%ns != 0) {
		nfull <- ceiling(n/ns)*ns
		npo <- n+1
		x[npo:nfull] <- as.single(-99999.0)
		warning(paste("The original series of", n,
									"values was padded with missing values to a length of", nfull,
									"so that it contains an integral number of years."))
		n <- as.integer(nfull)
	}
	
	x[is.na(x)] <- as.single(-99999.0)   # Convert any NAs to -99999.0
	
	# Call the Fortran seaken DLL
	# x is the equally spaced seasonal time series to examine for trend
	# n is the length of the timeseries
	# ns is the number of seasons in the timeseries
	# results is the statistics of the trend
	results <- .Fortran("seakenf", x , n , ns , as.single(vector("numeric", 9)))[[4]]
	## Return the statistics.
	method <- "Seasonal Kendall with correlation correction"
	tau <- results[1]
	names(tau) <- "tau"
	zero <- 0
	names(zero) <- "slope"
	est <- c(results[4], results[6],  results[5])
	names(est) <- c("slope", "median.data", "median.time")
	## use appropriate p.value
	if(n/ns < 10)
		p.value <- results[2]
	else
		p.value <- results[3]
	x <- miss2na(as.numeric(x), -99999.0)
	z <- list(method = method, data.name =
							paste(Dname, " (", n/ns, " years and ", ns, " seasons)", sep=""),
						nyears = n/ns, nseasons = ns, series=x, seasonames=nseas,
						statistic=tau, p.value=p.value,
						p.value.raw = results[2], p.value.corrected = results[3],
						estimate=est, alternative = "two.sided", null.value = zero)
	oldClass(z) <- c("htest", "seaken")
	return(z)
}

