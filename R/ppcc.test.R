#' Test for Normality
#' 
#' Computes the probability plot correlation coefficient test for departures
#' from normality.
#' 
#' @param x a vector of numeric values. Missing values are allowed, but are
#' ignored in the calculation.
#' @return An object of class "htest" having the following components:
#' \item{statistic}{ the value of the test statistic. } \item{p.value}{ the
#' attained p-value for the test. } \item{data.name}{ a character string
#' describing the name of the data used in the test. } \item{method}{ a
#' description of the method. }
#' @note The PPCC test is attractive because it has a simple, graphical
#' interpretation: it is a measure of the correlation in a Q-normal plot of the
#' data.  As such, it is related to the Shapiro-Wilk test (Shapiro and Wilk,
#' 1965) for normality.\cr
#' 
#' The distribution function of the test statistic is empirical. This
#' application uses the "pocket calculator" approximation for computing the
#' p-value of the observed statistic (Royston, 1992).\cr
#' @seealso \code{\link{shapiro.test}}
#' @references Filliben, 1975, The PPCC test for normality: Technometrics, v.
#' 17, no. 1, p. 111--117.\cr
#' 
#' Looney, S.W., and Gulledge, T.R., 1985, Use of the correlation coefficient
#' with normal probability plots: The American Statistician, v. 39, p.
#' 75--79.\cr
#' 
#' Royston, J.P., 1992, A pocket-calculator algorithm for the Shapiro-Francia
#' test of non-normality--an application to medicine: Statistics in Medicine,
#' v. 12, p. 181--184.\cr
#' 
#' Shapiro, S.S., and Wilk, M.B., 1965, An analysis of variance test for
#' normality (complete samples): Biometrika, v. 52, p. 591--611.
#' @keywords htest
#' @examples
#' 
#' ## These data should produce an attained p-value less than 0.001
#' set.seed(45)
#' ppcc.test.data <- rnorm(32)
#' qqnorm(ppcc.test.data)
#' abline(mean(ppcc.test.data), sd(ppcc.test.data))
#' ppcc.test(ppcc.test.data)
#' 
#' @export ppcc.test
ppcc.test <- function(x) {
	# Coding History
	#    2008Aug21 DLLorenz Initial dated version
	#    2008Oct20 DLLorenz USGS library version
	#    2009Sep03 DLLorenz use Royston's approximation to the p-value of W
	#    2011Oct25 DLLorenz Update for package
	#    2014Dec22 DLLorenz Roxygen header
  ##
  ## Ref is Filliben, 1975, Th PPCC test for normality, Technometrics,
  ## Vol 17, no. 1, p 111-117
  ##
  data.name <- deparse(substitute(x))
  x <- sort(x)
  pp <- qnorm(ppoints(x, a=.375))
  Stat <- cor(x, pp)
  names(Stat) <- "r"
  ## use Royston's 1992 approximation to W'
  W <- log(1 - Stat^2)
  N <- length(x)
  Adjm <- log(log(N)) - log(N)
  Adjs <- log(log(N)) + 2/log(N)
  z <- (W - (-1.2725 + 1.0521 * Adjm))/(1.0308 - 0.26758 * Adjs)
  P.val <- 1. - pnorm(z)
  retval <- list(statistic=Stat, p.value=P.val,
                 data.name=data.name,
                 alternative="Data are not from a normal distribution\nnull hypothesis: Data are from a normal distribution",
                 method="PPCC Normality Test",
  							 x=x)
  oldClass(retval) <- c("htest", "ppcc")
  return(retval)
}
