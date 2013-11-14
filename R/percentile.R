# compute cumulative percent or percent exceedence at a value q
#
# Coding history:
#    2007Oct12 DLLorenz Initial Coding
#    2011Aug09 DLLorenz Conversion to R and create generic function
#    2011Oct25 DLLorenz Update for package
#    2013Apr16 DLLorenz Named percentile
#

percentile <- function(x, q, test='>=', na.rm=TRUE,
                            percent=TRUE, ...) {
  ## Arguments:
  ##  x (any object for which a method is available)
  ##  q (numeric vector) the numeric criterion
  ##  test (character scalar) the test
  ##  na.rm (logical scalar) remove missing values
  ##  percent (logical scalar) express result in percent
  ##  ... (dots) required for generic function
  ##
  UseMethod("percentile")
}

percentile.default <- function(x, q, test='>=', na.rm=TRUE,
                            percent=TRUE, ...) {
  ## Arguments:
  ##  x (numeric vector) the values to test
  ##  q (numeric vector) The numeric criterion
  ##  test (character scalar) the test
  ##  na.rm (logical scalar) remove missing values
  ##  percent (logical scalar) express result in percent
  ##  ... (dots) required for method function
  ##
  retval <- double(length(q))
  check <- test
  test <- get(test)
  if(na.rm)
    x <- x[!is.na(x)]
  N <- length(x)
  for(i in seq(along=q)) {
    Ntest <- sum(test(x, q[i]))
    retval[i] <- Ntest/N
  }
  if(percent) {
    retval <- retval *100
    names(retval) <- paste("Percent", check, q, sep=' ')
  }
  else
    names(retval) <- paste("Proportion", check, q, sep=' ')
  return(retval)
}
