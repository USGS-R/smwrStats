# Compute skewness
#
# Coding history:
#    Unknown   DLLorenz Original Coding 
#    2011Aug24 DLLorenz Conversion to R
#    2013Feb28 DLLorenz Tweak to handling method
#    2013Feb28          This version
#

skew <- function(x, na.rm=TRUE, method="fisher") {
  ## Arguments:
  ##  x (numeric vector) 
  ##  na.rm (logical scalar) remove missings?
  ##  method (character scalar) the method to use: "fisher" or "moments"
  ##
  method <- match.arg(method, c("fisher", "moments"))
  if(na.rm)
    x <- x[!is.na(x)]
  n <- length(x)
  mn <- mean(x)
  dif.x <- x - mn
  m2 <- sum(dif.x^2)/n
  m3 <- sum(dif.x^3)/n
  b1 <- (m3/(m2^(3/2)))
  if(method == "moments")
    return(b1)
  if(n < 3)
    g1 <- NA
  else
    g1 <- (sqrt(n * (n - 1)) * b1)/(n - 2)
  return(g1)
}
