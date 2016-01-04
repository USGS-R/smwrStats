#' Empirical Cumulative Percent
#' 
#' Computes the empirical cumulative percent or percent exceedance of observed data
#' for specific values.\cr
#' 
#' 
#' @aliases percentile percentile.default
#' @param x a numeric vector representing the observed data.
#' @param q a vector of quantiles for which the cumulative percent or percent
#' exceedence is desired.
#' @param test a character string indicating the test. The default value, '>=,'
#' is the percent equalling or exceeding the quantile and '<' would return the
#' cumulative percent below the quantile.
#' @param na.rm a logical value indication whether missing values (NAs) should
#' be removed or not. If na.rm is \code{FALSE} and there are missing values in
#' \code{x}, then the result will be NA.  The default value is \code{TRUE}.
#' @param percent a logical value indicating whether the result should be
#' expressed as a percent or proportion. The default value, \code{TRUE}, will
#' express the result as a percent.
#' @param \dots not used, required for method function
#' @return A named vector as long as \code{q} corresponding to the requested
#' value.
#' @note The stats package contains the \code{ecdf} function that performs a
#' similar function when \code{test} is "<=."
#' @seealso \code{\link{ecdf}}
#' @keywords univar math manip
#' @examples
#' 
#' set.seed(2342)
#' Xr <- rlnorm(24)
#' # The percentage of the observarions greater than or equal to 2
#' percentile(Xr, 1:5)
#' 
#' @export percentile
percentile <- function(x, q, test='>=', na.rm=TRUE,
                            percent=TRUE, ...) {
	# Coding history:
	#    2007Oct12 DLLorenz Initial Coding
	#    2011Aug09 DLLorenz Conversion to R and create generic function
	#    2011Oct25 DLLorenz Update for package
	#    2013Apr16 DLLorenz Named percentile
	#    2014Dec22 DLLorenz Roxygen headers
  ##
  UseMethod("percentile")
}

#' @rdname percentile
#' @export 
#' @method percentile default
percentile.default <- function(x, q, test='>=', na.rm=TRUE,
                            percent=TRUE, ...) {
  ## Arguments:
  ##  x (numeric vector) the values to test
  ##  q (numeric vector) The numeric criterion
  ##  test (character scalar) the test
  ##  na.rm (logical scalar) remove missing values
  ##  percent (logical scalar) express result in percent
  ##  ... (dots) not used, required for method function
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
