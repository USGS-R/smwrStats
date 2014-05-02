# Compute the prob. plot corr. coef. test
#
# Coding History
#    2008Aug21 DLLorenz Initial dated version
#    2008Oct20 DLLorenz USGS library version
#    2009Sep03 DLLorenz use Royston's approximation to the p-value of W
#    2011Oct25 DLLorenz Update for package
#    2012Feb03          This version.
#

ppcc.test <- function(x) {
  ## Argument:
  ##  x (numeric vector) the data to tesst for normality
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
                 alternative="Data are from a normal distribution\nnull hypothesis: Data are not from a normal distribution",
                 method="PPCC Normality Test",
  							 x=x)
  oldClass(retval) <- c("htest", "ppcc")
  return(retval)
}
