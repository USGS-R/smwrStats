#' Test for Outlier
#' 
#' Test for a single outlier in a sample from a normal distribution.
#' 
#' 
#' @param x a vector of numeric values. Missing values are allowed, but are
#' ignored in the calculation.
#' @param alternate the alternate hypothesis; must be "two.sided" for either a
#' high or a low outlier, "high" to test only for a high outlier, or "low" to
#' test only for a lopw outlier.
#' @return An object of class "htest" having the following components:
#' \item{statistic}{ the value of the test statistic. } \item{p.value}{ the
#' attained p-value for the test. } \item{data.name}{ a character string
#' describing the name of the data used in the test. } \item{alternative}{ a
#' description of the altrnative and null hypotheses. } \item{method}{ a
#' description of the method. }
#' @seealso \code{\link{pgrubbs}}, \code{\link{qgrubbs}}
#' @references Grubbs, F., 1969, Procedures for Detecting Outlying Observations
#' in Samples, Technometrics, v. 11, no. 1, pp. 1-21.
#' @keywords htest
#' @examples
#' 
#' # A random sample with a high value
#' set.seed(100)
#' grubbstest <- rnorm(32)
#' grubbs.test(grubbstest)
#' qqnorm(grubbstest)
#' 
#' @export grubbs.test
grubbs.test <- function(x, alternate="two.sided") {
	## Arguments:
	##  x the data to test, missing values ignored
	##  alternate (character string): "two.sided", "high", or "low"
	## See http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm
	##  for a good explanaiton and details on the computation.
	##
	## Scale x
	x.name <- deparse(substitute(x))
	x <- scale(x)
	alternate <- match.arg(alternate, c("two.sided", "high", "low"))
	G <- switch(alternate,
							two.sided = max(abs(x), na.rm=TRUE),
							high = max(x, na.rm=TRUE),
							low = -min(x, na.rm=TRUE))
	names(G) <- "G"
	N <- sum(!is.na(x))
	p.val <- pgrubbs(G, N)
	## Need to protect against bizarre cases, like testing for a low outlier
	## when there is a very high outlier.
	if(alternate == "two.sided") {
		p.val <- min(p.val*2, 1)
		alt="Either a high or low outlier in the sample\nnull hypothesis: No outlier"
	}
	else {
		p.val <- min(p.val, 1)
		alt=paste("A ", alternate,
							" outlier in the sample\nnull hypothesis: No outlier", sep='')
	}
	retval <- list(statistic=G, p.value=p.val, data.name=x.name,
								 alternative=alt, method="Grubbs' test for an outlier")
	oldClass(retval) <- "htest"
	return(retval)
}
