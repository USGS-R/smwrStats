# function to compute the (one-sided) Grubbs crtical value and test
#
# Coding history:
#    2008Apr24 DLLorenz Original
#    2012May11 DLLorenz Conversion to R
#    2012May21          This version.
#

qgrubbs <- function(alpha, N) {
  ## Arguments:
  ##  alpha (numeric vector) the significance level
  ##  N the number of values in the sample
  ##
  ## return the one-sided critical values for alpha
  ## to get two-sided divied alpha by 2
  if(N < 3)
    return(as.double(NA))
  G <- qt(alpha/N, N-2)
  return((N-1)/sqrt(N)*sqrt(G^2/(N-2+G^2)))
}

pgrubbs <- function(G, N) {
  ## Arguments:
  ##  G (numeric vector) the maximum or minimum scaled difference from the mean
  ##  N the number of values in the sample
  ##
  ## return the one-sided attained p-value for the observed statistic
  ## G = (Ymax - Ymean)/s or
  ## G = (Ymean - Ymin)/s
  TT <- sqrt((N-2)/((N-1)^2/(N*G^2) - 1))
  pt(-TT, N-2)*N
}

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

